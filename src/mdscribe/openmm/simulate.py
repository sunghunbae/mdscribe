import copy
import collections
import logging
import functools
import pathlib
import sys
import shutil

import openmm
import openmm.app
import openmm.unit

import parmed

import femto.md
import femto.md.config
import femto.md.restraints



_LOGGER = logging.getLogger(__name__)


try:
    import pyspark.sql
    import pyspark.storagelevel

    def simulate_states(
        spark: pyspark.sql.session.SparkSession,
        system: openmm.System,
        topology: parmed.Structure,
        states: list[dict[str, float]],
        stages: list[femto.md.config.SimulationStage],
        platform: femto.md.constants.OpenMMPlatform,
        reporter: openmm.app.statedatareporter.StateDataReporter | None = None,
        enforce_pbc: bool = False,
        ) -> list[openmm.State]:
        """Launching Apache Spark jobs to run multi-stage simulations of each state.

        Args:
            system: The system to simulate.
            topology: The topology of the system to simulate.
            states: The states of the system to simulate.
            stages: The stages to run.
            platform: The accelerator to use.
            reporter: The reporter to use to record system statistics such as volume and
                energy.
            enforce_pbc: Whether to enforce periodic boundary conditions when retrieving
                the final coordinates.

        Returns:
            The final coordinates at each state.
        """
        state_indexes = list(range(len(states)))
        _partial_spark_simulate_state = functools.partial(
            simulate_state_index,
            system=system, # <-- serializable
            topology=topology, # <-- serializable
            states=states, # <-- list of dictionaries
            stages=stages,
            platform=platform,
            reporter=reporter,
            enforce_pbc=enforce_pbc,
            )
        rdd = spark.sparkContext.parallelize(state_indexes)
        results = rdd.map(_partial_spark_simulate_state).collect()
        return results
except:
    pass 


def _set_temperature(integrator: openmm.Integrator, temperature: openmm.unit.Quantity):
    """Sets the temperature of an OpenMM integrator.

    Args:
        integrator: The integrator to set the temperature of.
        temperature: The temperature to set the integrator to.
    """

    if hasattr(integrator, "setTemperature"):
        integrator.setTemperature(temperature.value_in_unit(openmm.unit.kelvin))
    else:
        integrator.setGlobalVariableByName(
            "kT", openmm.unit.MOLAR_GAS_CONSTANT_R * temperature
        )


def anneal_temperature(
    simulation: openmm.app.Simulation,
    temperature_initial: openmm.unit.Quantity,
    temperature_final: openmm.unit.Quantity,
    n_steps: int,
    frequency: int,
):
    """Gradually ramp the system temperature from a starting value to the final value.

    Args:
        simulation: The current simulation.
        temperature_initial: The initial temperature.
        temperature_final: The final temperature.
        n_steps: The number of steps to anneal over.
        frequency: The frequency at which to increment the temperature.
    """

    n_increments = n_steps // frequency
    increment = (temperature_final - temperature_initial) / n_increments

    temperature = temperature_initial
    _set_temperature(simulation.integrator, temperature)

    for _ in range(n_increments):
        simulation.step(frequency)

        temperature += increment
        _set_temperature(simulation.integrator, temperature)

    _set_temperature(simulation.integrator, temperature_final)



def create_position_restraints(
    topology: parmed.Structure,
    mask: str,
    config: femto.md.config.FlatBottomRestraint,
) -> openmm.CustomExternalForce:
    """Creates position restraints for all ligand atoms.

    Args:
        topology: The topology of the system being simulation.
        mask: The mask that defines which atoms to restrain.
        config: The restraint configuration.

    Returns:
        The created restraint force.
    """

    selection = parmed.amber.AmberMask(topology, mask).Selection()
    selection_idxs = tuple(i for i, matches in enumerate(selection) if matches)

    assert len(selection_idxs) > 0, "no atoms were found to restrain"

    coords = {
        i: openmm.Vec3(topology.atoms[i].xx, topology.atoms[i].xy, topology.atoms[i].xz)
        * openmm.unit.angstrom
        for i in selection_idxs
    }

    if not isinstance(config, femto.md.config.FlatBottomRestraint):
        raise NotImplementedError("only flat bottom restraints are currently supported")

    return femto.md.restraints.create_flat_bottom_restraint(config, coords)


def create_simulation(
    system: openmm.System,
    topology: parmed.Structure,
    coords: openmm.State | None,
    integrator: openmm.Integrator,
    state: dict[str, float] | None,
    platform: femto.md.constants.OpenMMPlatform,
) -> openmm.app.Simulation:
    """Creates an OpenMM simulation object.

    Args:
        system: The system to simulate
        topology: The topology being simulated.
        coords: The initial coordinates and box vectors. If ``None``, the coordinates
            and box vectors from the topology will be used.
        integrator: The integrator to evolve the system with.
        state: The state of the system to simulate.
        platform: The accelerator to run using.

    Returns:
        The created simulation.
    """

    platform_properties = (
        {"Precision": "mixed"} if platform.upper() in ["CUDA", "OPENCL"] else {}
    )
    platform = openmm.Platform.getPlatformByName(platform)

    if coords is not None:
        system.setDefaultPeriodicBoxVectors(*coords.getPeriodicBoxVectors())
    else:
        system.setDefaultPeriodicBoxVectors(*topology.box_vectors)

    simulation = openmm.app.Simulation(
        topology.topology, system, integrator, platform, platform_properties
    )

    if coords is None:
        simulation.context.setPeriodicBoxVectors(*topology.box_vectors)
        simulation.context.setPositions(topology.positions)
    else:
        simulation.context.setState(coords)

    state = femto.md.utils.openmm.evaluate_ctx_parameters(state, simulation.system)

    for k, v in state.items():
        simulation.context.setParameter(k, v)

    return simulation


def prepare_simulation(
    system: openmm.System,
    topology: parmed.Structure,
    state: dict[str, float],
    coords: openmm.State | None,
    config: femto.md.config.SimulationStage,
    platform: femto.md.constants.OpenMMPlatform,
    ) -> openmm.app.Simulation:
    """Prepare an OpenMM simulation object ready for a given stage."""
    system = copy.deepcopy(system)

    for mask, restraint in config.restraints.items():
        system.addForce(create_position_restraints(topology, mask, restraint))
    if isinstance(config, femto.md.config.Simulation) and config.pressure is not None:
        barostat = openmm.MonteCarloBarostat(
            config.pressure, config.temperature, config.barostat_frequency
        )
        system.addForce(barostat)

    if isinstance(config, femto.md.config.Anneal):
        integrator = femto.md.utils.openmm.create_integrator(
            config.integrator, config.temperature_initial
        )
    elif isinstance(config, femto.md.config.Simulation):
        integrator = femto.md.utils.openmm.create_integrator(
            config.integrator, config.temperature
        )
    else:
        integrator = openmm.VerletIntegrator(0.0001)

    femto.md.utils.openmm.assign_force_groups(system)

    simulation = create_simulation(system, topology, coords, integrator, state, platform)
    return simulation




def simulate_state_index(
    state_index: int,
    system: openmm.System,
    topology: parmed.Structure,
    states: list[dict[str, float]],
    stages: list[femto.md.config.SimulationStage],
    platform: femto.md.constants.OpenMMPlatform,
    reporter: openmm.app.statedatareporter.StateDataReporter | None = None,
    enforce_pbc: bool = False,
    ) -> openmm.State:
    """Simulate a system following the specified ``stages``, at a given 'state' (i.e.
    a set of context parameters, such as free energy lambda values)

    Args:
        system: The system to simulate.
        topology: The topology to simulate.
        state: The state to simulate at.
        stages: The stages to run.
        platform: The accelerator to use.
        reporter: The reporter to use to record system statistics such as volume and
            energy.
        enforce_pbc: Whether to enforce periodic boundary conditions when retrieving
            the final coordinates.

    Returns:
        The final coordinates and box vectors.
    """

    reporter = (
        reporter
        if reporter is not None
        else 
        openmm.app.statedatareporter.StateDataReporter(
            sys.stdout, 
            reportInterval=1000, 
            step=True, 
            time=True,
            potentialEnergy=True, 
            temperature=True, 
            volume=True, 
            density=True,
            )
        # else femto.md.reporting.openmm.OpenMMStateReporter(
        #     femto.md.reporting.NullReporter(), "", 999999999
        # )
    )

    state = states[state_index]
    stage_counter = collections.defaultdict(int)

    coords = None

    for stage in stages:
        stage_name = f"{stage.type}-{stage_counter[stage.type]}"
        stage_counter[stage.type] += 1

        reporter.tag = f"equilibration/{stage_name}"

        simulation = prepare_simulation(
            system, topology, state, coords, stage, platform
        )
        simulation.reporters.append(reporter)

        if isinstance(stage, femto.md.config.Minimization):
            reporter.report(simulation, simulation.context.getState(getEnergy=True))
            simulation.minimizeEnergy(
                stage.tolerance.value_in_unit(openmm.unit.kilojoules_per_mole / openmm.unit.nanometer), 
                stage.max_iterations
            )
            reporter.report(simulation, simulation.context.getState(getEnergy=True))
        elif isinstance(stage, femto.md.config.Anneal):
            anneal_temperature(
                simulation,
                stage.temperature_initial,
                stage.temperature_final,
                stage.n_steps,
                stage.frequency,
            )
        elif isinstance(stage, femto.md.config.Simulation):
            simulation.step(stage.n_steps)
        else:
            raise NotImplementedError(f"unknown stage type {type(stage)}")

        _LOGGER.info(
            f"after {stage_name} "
            f"{femto.md.utils.openmm.get_simulation_summary(simulation)}"
        )
        coords = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getForces=True,
            getEnergy=True,
            enforcePeriodicBox=enforce_pbc,
        )

    return coords




def pre_equilibrate(
    system: openmm.System,
    parmedstruct: parmed.Structure,
    statedict: dict[str, float],
    stages: list[femto.md.config.SimulationStage],
    platform: femto.md.constants.OpenMMPlatform = femto.md.constants.OpenMMPlatform.CUDA,
    enforce_pbc: bool = False,
    workdir: pathlib.Path | None = pathlib.Path("."),
    destdir: pathlib.Path | None = None,
    save_prefix: str = "equi",
    save_checkpoint: bool = True,
    save_state: bool = True,
    save_coord: bool = True,
    ) -> openmm.State:
    """Simulate a system following the specified ``stages``, at a given 'state' (i.e.
    a set of context parameters, such as free energy lambda values)

    Args:
        system (openmm.System): system to simulate
        parmedstruct (parmed.Structure): topology (parameters + coordinates) to simulate
        statedict (dict[str, float]): state dictionary
        stages (list[femto.md.config.SimulationStage]): stages for equilibration
        platform (femto.md.constants.OpenMMPlatform): computing platform. Defaults to CUDA.
        enforce_pbc (bool, optional): whether to enforce periodic boundary condition for final coordinates. Defaults to False.
        save_checkpoint (pathlib.Path | None, optional): save checkpoint (.cpt). Defaults to None.
        save_state (pathlib.Path | None, optional): save state (.xml). Defaults to None.
        save_coord (pathlib.Path | None, optional): save coordinates (.rst7). Defaults to None.

    Raises:
        NotImplementedError: unknown stage type

    Returns:
        openmm.State: final openmm state (position, velocity, force, energy)

    NVT annealing:
        When performing a simulated annealing process in an OpenMM molecular dynamics simulation, 
        it is generally recommended to use the NVT ensemble (constant number of particles, volume, and temperature), 
        as this allows for precise control of the system's temperature throughout the heating and cooling stages, 
        which is crucial for a successful annealing procedure.
    """

    try:
        workdir.mkdir(parents=True, exist_ok=True)
    except:
        print(f"cannot create given workdir: {workdir}")
        sys.exit(0)

    log_filename = f'{save_prefix}.log'
    checkpoint_filename = f'{save_prefix}.cpt'
    state_filename = f'{save_prefix}.xml'
    rst7_filename = f'{save_prefix}.rst7'
    prmtop_filename = f'{save_prefix}.prmtop'
    
    logging.basicConfig(filename= (workdir / log_filename).as_posix(),
                        filemode='w',
                        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
                        datefmt='%Y-%m-%d %H:%M',
                        level=logging.INFO)
    
    reporters = [
        openmm.app.statedatareporter.StateDataReporter(
            sys.stdout, 
            reportInterval=1000, 
            step=True, 
            time=True,
            potentialEnergy=True, 
            temperature=True, 
            volume=True, 
            density=True,
            ),
        openmm.app.statedatareporter.StateDataReporter(
            (workdir / log_filename).as_posix(), 
            reportInterval=1000, 
            step=True, 
            time=True,
            potentialEnergy=True, 
            temperature=True, 
            volume=True, 
            density=True,
            ),
    ]

    stage_state = None
    for stage in stages:
        simulation = prepare_simulation(system, parmedstruct, statedict, stage_state, stage, platform)
        for reporter in reporters:
            simulation.reporters.append(reporter)

        if isinstance(stage, femto.md.config.Minimization):   
            simulation.minimizeEnergy(
                stage.tolerance.value_in_unit(
                    openmm.unit.kilojoules_per_mole / openmm.unit.nanometer
                    ), 
                stage.max_iterations
            )
        elif isinstance(stage, femto.md.config.Anneal):
            anneal_temperature(
                simulation,
                stage.temperature_initial,
                stage.temperature_final,
                stage.n_steps,
                stage.frequency,
            )
        elif isinstance(stage, femto.md.config.Simulation):
            simulation.step(stage.n_steps)
        else:
            raise NotImplementedError(f"unknown stage type {type(stage)}")

        _LOGGER.info(f"{femto.md.utils.openmm.get_simulation_summary(simulation)}")

        stage_state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getForces=True,
            getEnergy=True,
            enforcePeriodicBox=enforce_pbc,
        )
    
    if save_checkpoint:
        with open(workdir / checkpoint_filename, "wb") as file:
            simulation.saveCheckpoint(file)

    if save_state:
        with open(workdir / state_filename, "w") as file:
            simulation.saveState(file)
    
    if save_coord:
        _parmedstruct = copy.deepcopy(parmedstruct)
        _parmedstruct.coordinates = stage_state.getPositions(asNumpy=True)
        _parmedstruct.save((workdir / rst7_filename).as_posix(), overwrite=True)
        _parmedstruct.save((workdir / prmtop_filename).as_posix(), overwrite=True)

    # copy output files in the `workdir` to the `destdir``
    if destdir is not None:
        destdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(workdir / log_filename, destdir)
        if save_checkpoint:
            shutil.copy(workdir / checkpoint_filename, destdir)
        if save_state:
            shutil.copy(workdir / state_filename, destdir)
        if save_coord:
            shutil.copy(workdir / rst7_filename, destdir)
            
    return stage_state