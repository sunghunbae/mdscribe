""" HREMD with PySpark """

import femto.md.config
import femto.md.constants
import femto.md.utils.openmm

import parmed
import openmm
import numpy
import pyarrow
import functools
import typing
import contextlib
import pathlib
import pickle
import itertools
import logging
import shutil
import pandas
import sys
import copy

import pyspark.sql
import pyspark.storagelevel

_LOGGER = logging.getLogger(__name__)

_HREMDStorage = typing.NamedTuple(
    "HREMDStorage",
    [
        ("file", pyarrow.OSFile),
        ("writer", pyarrow.RecordBatchStreamWriter),
        ("schema", pyarrow.Schema),
    ],
)

@contextlib.contextmanager
def _create_exchanges_storage(output_path: pathlib.Path | None, n_states: int) -> _HREMDStorage | None:
    """Open a storage ready for writing.

    Args:
        output_path (pathlib.Path | None): The path to write the output to.
        n_states (int): number of states

    Returns:
        _HREMDStorage | None: OpenMM report object or None

    Yields:
        Iterator[_HREMDStorage, None]: OpenMM report object or None
    """

    if output_path is None:
        yield None
        return

    schema = pyarrow.schema([
        ("step", pyarrow.int64()),
        ("u_kn", pyarrow.list_(pyarrow.list_(pyarrow.float64(), n_states), n_states),),
        ("replica_to_state_idx", pyarrow.list_(pyarrow.int16(), n_states),),
        ("n_proposed_swaps", pyarrow.list_(pyarrow.list_(pyarrow.int64(), n_states), n_states),),
        ("n_accepted_swaps", pyarrow.list_(pyarrow.list_(pyarrow.int64(), n_states), n_states),),
        ])

    output_path.parent.mkdir(exist_ok=True, parents=True)

    records = []

    if output_path.exists():
        # pyarrow doesn't seem to offer a clean way to stream append
        with pyarrow.OSFile(str(output_path), "rb") as file:
            with pyarrow.RecordBatchStreamReader(file) as reader:
                for record in reader:
                    records.append(record)

    with pyarrow.OSFile(str(output_path), "wb") as file:
        with pyarrow.RecordBatchStreamWriter(file, schema) as writer:
            for record in records:
                writer.write_batch(record)

            yield _HREMDStorage(file, writer, schema)


def _create_trajectory_storage(
    simulation: openmm.app.Simulation,
    n_replicas: int,
    n_steps_per_cycle: int,
    trajectory_interval: int | None,
    trajectory_dir: pathlib.Path | None,
    trajectory_prefix: str,
    exit_stack: contextlib.ExitStack,
    ) -> list[openmm.app.DCDFile] | None:
    """Open a DCD trajectory reporter per replica.

    Args:
        simulation: The simulation to report.
        n_replicas: The number of replicas being sampled on this process.
        n_steps_per_cycle: The number of steps per cycle.
        trajectory_interval: The interval with which to write the trajectory.
        output_dir: The root output directory. 
        output_prefix: Trajectories will be written to `{trajectory_dir}/{trajectory_prefix}{replica_idx}.dcd`.
        exit_stack: The exit stack to use for opening the files.

    Returns:
        The trajectory reporters if the trajectory interval is greater than zero.
    """
    if trajectory_dir is None or trajectory_interval is None or trajectory_interval <= 0:
        return

    trajectory_dir.mkdir(exist_ok=True, parents=True)
    trajectory_paths = [ trajectory_dir / f"{trajectory_prefix}{i}.dcd" for i in range(n_replicas) ]
    should_append = [path.exists() for path in trajectory_paths]

    return [
        openmm.app.DCDFile(
            exit_stack.enter_context(path.open("wb" if not append else "r+b")),
            simulation.topology,
            simulation.integrator.getStepSize(),
            n_steps_per_cycle * trajectory_interval,
            append=append,
        )
        for path, append in zip(trajectory_paths, should_append, strict=True)
    ]


def _store_trajectory(coords: list[openmm.State], storage: list[openmm.app.DCDFile] | None):
    """Store the current replica states to DCD files."""

    if storage is None:
        return

    for state, output_file in zip(coords, storage, strict=True):
        output_file.writeModel(
            positions=state.getPositions(),
            periodicBoxVectors=state.getPeriodicBoxVectors(),
        )


def _store_potentials(
    replica_to_state_idx: numpy.ndarray,
    reduced_potentials: numpy.ndarray,
    n_proposed_swaps: numpy.ndarray,
    n_accepted_swaps: numpy.ndarray,
    storage: _HREMDStorage | None,
    step: int,
    ) -> None:
    """Report the current state of the replica exchange simulation to an output file.

    Args:
        replica_to_state_idx (numpy.ndarray): _description_
        reduced_potentials (numpy.ndarray): _description_
        n_proposed_swaps (numpy.ndarray): _description_
        n_accepted_swaps (numpy.ndarray): _description_
        storage (_HREMDStorage | None): _description_
        step (int): _description_
    """

    if storage is None:
        return

    record = pyarrow.record_batch(
        [
            (step,),
            (reduced_potentials.tolist(),),
            (replica_to_state_idx,),
            (n_proposed_swaps.tolist(),),
            (n_accepted_swaps.tolist(),),
        ],
        schema=storage.schema,
    )
    storage.writer.write_batch(record)


def _store_checkpoint(
    start_cycle: int,
    coords: list[openmm.State],
    u_kn: numpy.ndarray,
    n_k: numpy.ndarray,
    has_sampled: numpy.ndarray,
    n_proposed_swaps: numpy.ndarray,
    n_accepted_swaps: numpy.ndarray,
    replica_to_state_idx: numpy.ndarray,
    path: pathlib.Path,
    ) -> None:
    """Store the state of an HREMD simulation to a pickle checkpoint.

    Args:
        start_cycle (int): _description_
        coords (list[openmm.State]): _description_
        u_kn (numpy.ndarray): _description_
        n_k (numpy.ndarray): _description_
        has_sampled (numpy.ndarray): _description_
        n_proposed_swaps (numpy.ndarray): _description_
        n_accepted_swaps (numpy.ndarray): _description_
        replica_to_state_idx (numpy.ndarray): _description_
        path (pathlib.Path): path to the checkpoint file
    """

    coords_dict = {i : coord for i, coord in enumerate(coords)}
    coords = [coords_dict[replica_to_state_idx[i]] for i in range(len(u_kn))]

    path.parent.mkdir(exist_ok=True, parents=True)

    with path.open("wb") as file:
        pickle.dump(
            (
                start_cycle,
                coords,
                u_kn,
                n_k,
                has_sampled,
                n_proposed_swaps,
                n_accepted_swaps,
                replica_to_state_idx,
            ),
            file,
        )

def _load_checkpoint(config: femto.md.config.HREMD, n_states: int, path: pathlib.Path) -> tuple[
    int,
    list[openmm.State],
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    ]:
    """Attempt to load a previous HREMD run from a pickled checkpoint.

    Returns:
        The starting cycle, initial coordinates for all states, u_kn, N_k, which states
        have been sampled, the number of proposed swaps, the number of accepted swaps,
        and the replica to state index map.
    """
    (
        start_cycle,
        initial_coords,
        u_kn,
        n_k,
        has_sampled,
        n_proposed_swaps,
        n_accepted_swaps,
        replica_to_state_idx,
    ) = pickle.loads(path.read_bytes())

    n_found_cycles = u_kn.shape[1] // n_states

    assert n_found_cycles <= config.n_cycles

    if u_kn.shape[1] != (config.n_cycles * n_states):
        n_pad = config.n_cycles - n_found_cycles

        u_kn_block = u_kn.reshape(n_states, n_states, n_found_cycles)
        u_kn_block = numpy.pad(
            u_kn_block, ((0, 0), (0, 0), (0, n_pad)), "constant", constant_values=0
        )

        u_kn = u_kn_block.reshape(n_states, n_states * config.n_cycles)

        has_sampled_block = has_sampled.reshape(n_states, n_found_cycles)
        has_sampled_block = numpy.pad(
            has_sampled_block, ((0, 0), (0, n_pad)), "constant", constant_values=False
        )

        has_sampled = has_sampled_block.reshape(n_states * config.n_cycles)

    return (
        start_cycle,
        initial_coords,
        u_kn,
        n_k,
        has_sampled,
        n_proposed_swaps,
        n_accepted_swaps,
        replica_to_state_idx,
    )


def _propose_swap(
    state_idx_i: int,
    state_idx_j: int,
    reduced_potentials: numpy.ndarray,
    n_proposed_swaps: numpy.ndarray,
    n_accepted_swaps: numpy.ndarray,
    replica_to_state_idx: numpy.ndarray,
    ):
    """Attempt to swap a pair of states between replicas.

    Notes:
        replica_to_state_idx
    """

    state_to_replica_idx = {
        state_idx: replica_idx
        for replica_idx, state_idx in enumerate(replica_to_state_idx)
    }

    replica_idx_i = state_to_replica_idx[state_idx_i]
    replica_idx_j = state_to_replica_idx[state_idx_j]

    potential_ii = reduced_potentials[state_idx_i, replica_idx_i]
    potential_ij = reduced_potentials[state_idx_i, replica_idx_j]
    potential_ji = reduced_potentials[state_idx_j, replica_idx_i]
    potential_jj = reduced_potentials[state_idx_j, replica_idx_j]

    criteria = -((potential_ij - potential_ii) - (potential_jj - potential_ji))

    if any(
        numpy.isnan(x) for x in [potential_ii, potential_ij, potential_ji, potential_jj]
    ):
        return

    n_proposed_swaps[state_idx_i, state_idx_j] += 1
    n_proposed_swaps[state_idx_j, state_idx_i] += 1

    if criteria >= 0.0 or numpy.random.rand() < numpy.exp(criteria):
        replica_to_state_idx[replica_idx_i] = state_idx_j
        replica_to_state_idx[replica_idx_j] = state_idx_i

        n_accepted_swaps[state_idx_i, state_idx_j] += 1
        n_accepted_swaps[state_idx_j, state_idx_i] += 1


def _propose_swaps(
    replica_to_state_idx: numpy.ndarray,
    reduced_potentials: numpy.ndarray,
    n_proposed_swaps: numpy.ndarray,
    n_accepted_swaps: numpy.ndarray,
    mask: set[tuple[int, int]],
    mode: femto.md.config.HREMDSwapMode | None,
    max_swaps: int | None,
    ):
    """Attempt to swap states between replicas.

    Args:
        replica_to_state_idx: The replica to state index map to modify in-place.
        reduced_potentials: The matrix of reduced potentials with
            ``reduced_potentials[state_idx][replica_idx] = value``
        n_proposed_swaps: An array tracking the number of proposed swaps to modify
            in-place.
        n_accepted_swaps: An array tracking the number of accepted swaps to modify
            in-place.
        mask: Pairs of state indices that cannot be swapped.
        mode: The swap mode. This can either be:
            * Neighbours: only try and swap adjacent states
            * All: try and swap all states stochastically
        max_swaps: The maximum number of proposals to make if running in 'all' mode.
            This variable does nothing when running in 'neighbours' mode.
    """

    if mode is None:
        return

    n_states = len(replica_to_state_idx)

    if mode == femto.md.config.HREMDSwapMode.NEIGHBOURS:
        # select whether to swap [0, 1], [2, 3], ... OR [1, 2], [3, 4], ...
        state_idx_offset = numpy.random.randint(2)

        pairs = [
            (state_idx_i, state_idx_i + 1)
            for state_idx_i in range(state_idx_offset, n_states - 1, 2)
            if (state_idx_i, state_idx_i + 1) not in mask
            and (state_idx_i + 1, state_idx_i) not in mask
        ]

    elif mode == femto.md.config.HREMDSwapMode.ALL:
        pairs = [
            (state_idx_i, state_idx_j)
            for state_idx_i, state_idx_j in itertools.combinations(range(n_states), r=2)
            if (state_idx_i, state_idx_j) not in mask
            and (state_idx_j, state_idx_i) not in mask
        ]
        max_swaps = len(pairs) if max_swaps is None else max_swaps

        pairs = (
            pairs
            if len(pairs) <= max_swaps
            else [
                pairs[i]
                for i in numpy.random.random_integers(0, len(pairs) - 1, max_swaps)
            ]
        )

    else:
        raise NotImplementedError

    for state_idx_i, state_idx_j in pairs:
        _propose_swap(
            state_idx_i,
            state_idx_j,
            reduced_potentials,
            n_proposed_swaps,
            n_accepted_swaps,
            replica_to_state_idx,
        )


def _spark_propagate_replica(
    replica_idx: int,
    structure: parmed.Structure,
    system: openmm.System,
    integrator_config: femto.md.config.LangevinIntegrator,
    temperature: openmm.unit.Quantity,
    pressure: openmm.unit.Quantity | None,
    n_steps: int,
    states: list[dict[str, float]],
    replica_to_state_idx: numpy.ndarray,
    broadcast_coords: list[openmm.State],
    force_groups: set[int] | int,
    max_retries: int,
    enforce_pbc: bool = True,
    ) -> tuple[int, openmm.State, numpy.ndarray]:
    """Propagate replica states using pyspark.

    Args:
        replica_idx (int): index of replica
        structure (parmed.Structure): topology, parmed.Structure object.
        system (openmm.System): openmm.System object.
        integrator (openmm.Integrator): openmm.Integrator object.
        temperature (openmm.unit.Quantity): simulation temperature.
        pressure (openmm.unit.Quantity | None): simulation pressure.
        states (list[dict[str, float]]): states.
        coords (list[openmm.State]): coordinates.
        n_steps (int): number of simulation steps
        replica_to_state_idx (numpy.ndarray): mapping of replica index to state index
        force_groups (set[int] | int): force groups
        max_retries (int): max. attempts to step a replica if a NaN is encountered
        enforce_pbc (bool, optional): enforce PBC. Defaults to True.

    Returns:
        tuple[list, numpy.ndarray]: (replica_coords, reduced potentials)
    """
    state_idx = replica_to_state_idx[replica_idx]
    replica_state = broadcast_coords.value[replica_idx]

    simulation = femto.md.utils.openmm.create_simulation(
        system,
        structure,
        None,  # or None to use the coordinates / box in structure
        integrator=femto.md.utils.openmm.create_integrator(integrator_config, temperature),
        state=states[state_idx],
        platform=femto.md.constants.OpenMMPlatform.CUDA,
    )
    
    for attempt in range(max_retries):
        try:
            simulation.context.setState(replica_state)

            for key, value in states[state_idx].items():
                simulation.context.setParameter(key, value)

            simulation.step(n_steps)

            replica_state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=enforce_pbc,
            )
            femto.md.utils.openmm.check_for_nans(replica_state)
        
        except openmm.OpenMMException:
            # randomize the velocities and try again
            simulation.context.setVelocitiesToTemperature(temperature)
            message = f"NaN detected for replica={replica_idx} state={state_idx}"

            if attempt == max_retries - 1:
                _LOGGER.warning(f"{message} that could not be resolved by retries.")
                raise

            _LOGGER.warning(f"{message}, retrying {attempt + 1}/{max_retries}")
            continue

        break
        
    reduced_potentials = _compute_reduced_potentials(
        simulation.context, states, temperature, pressure, force_groups)
    
    return (replica_idx, replica_state, reduced_potentials)



def _propagate_replicas(
    spark: pyspark.sql.session.SparkSession,
    simulation: openmm.app.Simulation,
    structure: parmed.Structure,
    temperature: openmm.unit.Quantity,
    pressure: openmm.unit.Quantity | None,
    states: list[dict[str, float]],
    coords: list[openmm.State],
    n_steps: int,
    replica_to_state_idx: numpy.ndarray,
    force_groups: set[int] | int,
    max_retries: int,
    enforce_pbc: bool = True,
    ) -> numpy.ndarray:
    """Propagate all replica states forward in time.

    Args:
        simulation: The main simulation object.
        states: The states being sampled.
        temperature: The temperature being sampled at.
        pressure: The pressure being sampled at.
        coords: The starting coordinates of each replica.
        n_steps: The number of steps to propagate by.
        replica_to_state_idx: A map between each replica and the state it is sampling.
        force_groups: The force groups to consider. Use -1 if all groups should be
            considered.
        max_retries: The maximum number of times to attempt to step a replica if a NaN
            is encountered
    """
    n_states = len(states)
    reduced_potentials = numpy.zeros((n_states, n_states))

    if n_steps <= 0:
        return reduced_potentials

    """
    When using PySpark's map and collect with functions that have complex arguments, you might encounter pickling errors. 
    This is because PySpark needs to serialize your function and its arguments to distribute them across the cluster, 
    and some objects aren't serializable. In our case, `openmm.app.Simulation` object itselt is NOT serializable. 
    To resolve this issue, a new simulation object is created within the worker function, `_spark_propagate_replica`, 
    from serializable `structure`, `system`, and `integrator` objects.
    """
    # integrator object will be created within the `_spark_propagaate_replica()`
    timestep = simulation.integrator.getStepSize() # openmm.unit.Quantity
    
    # Broadcast the replica states or coordinates
    broadcast_coords = spark.sparkContext.broadcast(coords)

    _partial_spark_propagate_replica = functools.partial(
        _spark_propagate_replica, 
        structure=structure, # <-- serializable
        system=simulation.system, # <-- serializable
        integrator_config=femto.md.config.LangevinIntegrator(timestep=timestep), # <-- serializable
        temperature=temperature, 
        pressure=pressure,
        states=states, # <-- list of dictionaries
        n_steps=n_steps,
        replica_to_state_idx=replica_to_state_idx,
        broadcast_coords=broadcast_coords,
        force_groups=force_groups,
        max_retries=max_retries,
        enforce_pbc=enforce_pbc,
        )

    rdd = spark.sparkContext.parallelize(list(range(len(coords))))
    # rdd.persist(pyspark.storagelevel.StorageLevel.MEMORY_AND_DISK)
    results = rdd.map(_partial_spark_propagate_replica).collect()

    for replica_idx, replica_state, red_pot in results:
        reduced_potentials[:, replica_idx] = red_pot
        coords[replica_idx] = replica_state

    return reduced_potentials



def _compute_reduced_potentials(
    context: openmm.Context,
    states: list[dict[str, float]],
    temperature: openmm.unit.Quantity,
    pressure: openmm.unit.Quantity | None,
    force_groups: set[int] | int,
    ) -> numpy.ndarray:
    """Compute the reduced potential of the given coordinates for each state being
    sampled.

    Args:
        context: The current simulation config.
        states: The states being sampled.
        temperature: The temperature being sampled at.
        force_groups: The force groups to consider. Use -1 if all groups should be
            considered.

    Returns:
        The reduced potentials with ``shape=(n_states,)``.
    """
    beta = 1.0 / (openmm.unit.BOLTZMANN_CONSTANT_kB * temperature)

    reduced_potentials = numpy.zeros(len(states))

    for state_idx, state in enumerate(states):
        for key, value in state.items():
            context.setParameter(key, value)

        reduced_potential = (
            context.getState(getEnergy=True, groups=force_groups).getPotentialEnergy()
            / openmm.unit.AVOGADRO_CONSTANT_NA
        )

        if pressure is not None:
            reduced_potential += pressure * context.getState().getPeriodicBoxVolume()

        reduced_potentials[state_idx] = beta * reduced_potential

    return reduced_potentials



def _mean_exchange_acceptance(
        n_proposed_swaps:numpy.ndarray,
        n_accepted_swaps:numpy.ndarray,
        verbose:bool = False,
        ) -> list[tuple[int, int, float]]:
    """calculate mean acceptance ratio

    Args:
        n_proposed_swaps (numpy.ndarray): number of proposed swaps for [states_i, states_j]
        n_accepted_swaps (numpy.ndarray): number of accepted swaps for [states_i, states_j]

    Returns:
        float: mean acceptance ratio [0..1]
    """

    if numpy.count_nonzero(n_proposed_swaps) == 0:
        return

    non_zero_indices = numpy.nonzero(n_proposed_swaps)
    acceptance_ratios = n_accepted_swaps[non_zero_indices]/n_proposed_swaps[non_zero_indices]
    array_i, array_j = non_zero_indices
    
    retval = []
    for i, j, acc in zip(array_i, array_j, acceptance_ratios):
        retval.append((i, j, acc))
        message = f'{i:2d} -> {j:2d} exchange acceptance= {acc:.3f}'
        _LOGGER.info(message)
        if verbose:
            print(message)
    
    message = f'mean exchange acceptance= {numpy.mean(acceptance_ratios):.3f}'
    _LOGGER.info(message)
    if verbose:
        print(message)

    return retval



def run_hremd(
    simulation: openmm.app.Simulation,
    structure: parmed.Structure,
    states: list[dict[str, float]],
    config: femto.md.config.HREMD,
    swap_mask: set[tuple[int, int]] | None = None,
    force_groups: set[int] | int = -1,
    initial_coords: list[openmm.State] | None = None,
    analysis_fn: typing.Callable[[int, numpy.ndarray, numpy.ndarray], None] | None = None,
    analysis_interval: int | None = None,
    workdir: pathlib.Path | None = None,
    destdir: pathlib.Path | None = None, # checkpoint should be stored here
    report_acceptance: bool = True,
    checkpoint_filename: str = "checkpoint.pkl",
    exchanges_filename: str = "exchanges.arrow",
    trajectory_dirname: str = "trajectories",
    trajectory_prefix: str = "replica",
    log_filename: str = "hremd.log",
    spark: pyspark.sql.SparkSession | None = None,
    ) -> tuple[numpy.ndarray, numpy.ndarray, list[openmm.State]]:
    """Run a Hamiltonian replica exchange simulation.

    Args:
        simulation: The main simulation object to sample using.
        states: The states to sample at. This should be a dictionary with keys
            corresponding to global context parameters.
        config: The sampling configuration.
        workdir: The directory to store the sampled energies and statistics to, and
            any trajectory / checkpoint files if requested in the config. If ``None``,
            no output of any kind will be written.
        swap_mask: Pairs of states that should not be swapped.
        force_groups: The force groups to consider when computing the reduced potentials
        initial_coords: The initial coordinates of each state. If not provided, the
            coordinates will be taken from the simulation object.
        analysis_fn: A function to call after every ``analysis_interval`` cycles. It
            should take as arguments the current cycle number, the reduced potentials
            with ``shape=(n_states, n_samples)`` and the number of samples of each
            state with ``shape=(n_states,)``.
        analysis_interval: The interval with which to call the analysis function.
            If ``None``, no analysis will be performed.

    Returns:
        The reduced potentials, the number of samples of each state, and the final
        coordinates of each state.
    """

    try:
        workdir.mkdir(parents=True, exist_ok=True)
    except:
        print(f"cannot create given workdir: {workdir}")
        sys.exit(0)

    logging.basicConfig(filename= workdir / log_filename,
                        filemode='w',
                        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
                        datefmt='%Y-%m-%d %H:%M',
                        level=logging.INFO)

    n_states = len(states)
    _LOGGER.info(f"number of states: {n_states}")

    # states = [
    #     femto.md.utils.openmm.evaluate_ctx_parameters(state, simulation.system)
    #     for state in states
    # ]

    swap_mask = set() if swap_mask is None else swap_mask

    n_proposed_swaps = numpy.zeros((n_states, n_states))
    n_accepted_swaps = numpy.zeros((n_states, n_states))

    replica_to_state_idx = numpy.arange(n_states) # 0, 1, 2, .., n_states-1

    u_kn, n_k = (
        numpy.empty((n_states, n_states * config.n_cycles)),
        numpy.zeros(n_states, dtype=int),
    )
    has_sampled = numpy.zeros(n_states * config.n_cycles, bool)

    pressure = femto.md.utils.openmm.get_pressure(simulation.system)
    
    exchanges_path = None if workdir is None else workdir / exchanges_filename
    checkpoint_path = (
        None
        if destdir is None
        or config.checkpoint_interval is None
        or config.checkpoint_interval <= 0
        else destdir / checkpoint_filename
    )

    start_cycle = 0

    if checkpoint_path is not None and checkpoint_path.exists():
        (
            start_cycle,
            initial_coords,
            u_kn,
            n_k,
            has_sampled,
            n_proposed_swaps,
            n_accepted_swaps,
            replica_to_state_idx,
        ) = _load_checkpoint(config, n_states, checkpoint_path)
        _LOGGER.info(f"resuming from cycle {start_cycle} samples")

    with ( 
        _create_exchanges_storage(exchanges_path, n_states) as storage, 
        contextlib.ExitStack() as exit_stack,
        ):
        # contextlib.ExitStack() is to properly close the multiple files 
        # created by `_create_exchanges_storage` in this `with` statement

        if initial_coords is None:
            coords_0 = simulation.context.getState(
                getPositions=True, 
                enforcePeriodicBox=config.trajectory_enforce_pbc
            )
            coords = [coords_0] * n_states
        else:
            coords = initial_coords

        if start_cycle == 0:
            _LOGGER.info(f"running {config.n_warmup_steps} warm-up steps")

            _propagate_replicas(
                spark,
                simulation,
                structure,
                config.temperature,
                pressure,
                states,
                coords,
                config.n_warmup_steps,
                replica_to_state_idx,
                force_groups,
                config.max_step_retries,
                config.trajectory_enforce_pbc,
            )

        trajectory_storage = _create_trajectory_storage(
            simulation,
            n_states, # <-- n_replicas
            config.n_steps_per_cycle,
            config.trajectory_interval,
            workdir / trajectory_dirname,
            trajectory_prefix, # {trajectory_prefix}{i}.dcd 
            exit_stack,
        )

        _LOGGER.info(f"running {config.n_cycles} replica exchange cycles")

        for cycle in range(start_cycle, config.n_cycles):
            reduced_potentials = _propagate_replicas(
                spark,
                simulation,
                structure,
                config.temperature,
                pressure,
                states,
                coords,
                config.n_steps_per_cycle,
                replica_to_state_idx,
                force_groups,
                config.max_step_retries,
                config.trajectory_enforce_pbc,
            )

            has_sampled[replica_to_state_idx * config.n_cycles + cycle] = True
            u_kn[:, replica_to_state_idx * config.n_cycles + cycle] = reduced_potentials

            n_k += 1

            should_save_trajectory = (
                config.trajectory_interval is not None
                and cycle % config.trajectory_interval == 0
            )

            if should_save_trajectory:
                # OSError: [Errno 95]
                # Cause - The underlying storage that is mounted to DBFS does not support append.
                # Solution - As a workaround, you should run your append on a local disk, 
                # such as /tmp, and move the entire file at the end of the operation.
                _store_trajectory(coords, trajectory_storage)

            should_analyze = (
                analysis_fn is not None
                and analysis_interval is not None
                and cycle % analysis_interval == 0
            )

            if should_analyze:
                analysis_fn(cycle, u_kn[:, has_sampled], n_k)

            _store_potentials(
                replica_to_state_idx,
                reduced_potentials,
                n_proposed_swaps,
                n_accepted_swaps,
                storage,
                cycle * config.n_steps_per_cycle,
            )

            _propose_swaps(
                replica_to_state_idx,
                reduced_potentials,
                n_proposed_swaps, # <-- to be changed
                n_accepted_swaps, # <-- to be changed
                swap_mask,
                config.swap_mode,
                config.max_swaps,
            )

            should_checkpoint = (
                checkpoint_path is not None
                and config.checkpoint_interval is not None
                and (
                    cycle % config.checkpoint_interval == 0
                    or cycle == config.n_cycles - 1
                )
            )

            if should_checkpoint:
                _store_checkpoint(
                    cycle + 1,
                    coords,
                    u_kn,
                    n_k,
                    has_sampled,
                    n_proposed_swaps,
                    n_accepted_swaps,
                    replica_to_state_idx,
                    checkpoint_path,
                )

            # report accumulated mean acceptance ratio
            mean_acc = 0.0
            if numpy.count_nonzero(n_proposed_swaps) != 0:
                non_zero_indices = numpy.nonzero(n_proposed_swaps)
                acc = n_accepted_swaps[non_zero_indices]/n_proposed_swaps[non_zero_indices]
                mean_acc = numpy.mean(acc)

            _LOGGER.info(f"done cycle {cycle+1} of {config.n_cycles}. mean acceptance ratio {mean_acc:.4f}")

        coords_dict = {i : coord for i, coord in enumerate(coords)}

        final_coords = [coords_dict[replica_to_state_idx[i]] for i in range(n_states)]
        
        # report mean acceptance ratio
        # `n_proposed_swaps` and `n_accepted_swaps` have accumulated counts
        # of proposed and accepted exchanges
        if report_acceptance:
            acc = _mean_exchange_acceptance(n_proposed_swaps, n_accepted_swaps, verbose=True)

        # copy output files in the `workdir` to the `destdir``
        if destdir is not None:
            
            destdir.mkdir(parents=True, exist_ok=True)

            # OSError: [Errno 95]
            # Cause - The underlying storage that is mounted to DBFS does not support append.
            # Solution - As a workaround, you should run your append on a local disk, 
            # such as /tmp, and move the entire file at the end of the operation.

            if (workdir / log_filename).is_file():
                shutil.copy(workdir / log_filename, destdir)
            
            if (workdir / exchanges_filename).is_file():
                shutil.copy(workdir / exchanges_filename, destdir)
            
            for file_path in (workdir / trajectory_dirname).iterdir():
                if file_path.is_file():
                    shutil.copy(file_path, destdir)

    return u_kn, n_k, final_coords


def plot_acceptance(exchanges_filename:str):
    """_summary_

    Args:
        exchanges_filename (str): exchanges data filename

    Returns:
        _type_: _description_
    """
    with pyarrow.OSFile(exchanges_filename, "rb") as file:
        with pyarrow.RecordBatchStreamReader(file) as reader:
            table = reader.read_all()
            S = table['step'].to_pylist()
            P = table['n_proposed_swaps'].to_pylist()
            A = table['n_accepted_swaps'].to_pylist()
            try:
                assert len(P) == len(A)
            except:
                print("something is wrong")
                return
            
            step = []
            mean_acc = []
            for s,p,a in zip(S,P,A):
                if numpy.count_nonzero(p) == 0:
                    continue
                non_zero_indices = numpy.nonzero(p)
                acc = a[non_zero_indices]/p[non_zero_indices]
                step.append(s)
                mean_acc.append(numpy.mean(acc))
                
            df = pandas.DataFrame({"step":step, "mean_acceptance":mean_acc})
            g = df.plot(x="step", 
                        y="mean_acceptance", 
                        xlabel='Step', 
                        ylabel='Mean Exchange Acceptance Ratio', 
                        xlim=0,
                        ylim=0,
                       )
            return g