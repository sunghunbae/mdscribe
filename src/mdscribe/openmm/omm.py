__all__ = [ 'OpenMMSetup', 'create_openmm_md_jobs' ]

from openmm.app import ForceField, PDBFile, Modeller, HBonds, PME, Simulation
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm.app import DCDReporter, CheckpointReporter, StateDataReporter
from openmm import LangevinIntegrator, CustomExternalForce, MonteCarloBarostat
from openmm.unit import kelvin, bar, nanosecond, picosecond, femtosecond
from openmm.unit import nanometer, angstrom
from openmm.unit import kilocalories_per_mole

import argparse
import pathlib
import sys
import logging
from operator import itemgetter


class OpenMMSetup(object):
    def __init__(self, args:argparse.Namespace,
                 prmtop_path:pathlib.Path, 
                 inpcrd_path:pathlib.Path,
                 outfile_prefix:str, 
                 logger:logging.Logger) -> None:

        # output files
        self.log_path = args.workdir / f'{outfile_prefix}.log'
        self.state_data_path = args.workdir / f'{outfile_prefix}.ene'
        self.checkpoint_path = args.workdir / f'{outfile_prefix}.cpt'
        self.trajectory_path = args.workdir / f'{outfile_prefix}.dcd'

        self.logger = logger
        self.logger_file_h = logging.FileHandler(self.log_path, mode="a", encoding='utf-8')
        self.logger_file_h.setFormatter(logging.Formatter(
            fmt='%(asctime)s %(message)s', datefmt='%Y/%m/%d %H:%M:%S'))
        self.logger.addHandler(self.logger_file_h)

        self.temperature = args.temperature * kelvin
        self.pressure = args.pressure * bar
        self.prod_time = args.prod * nanosecond
        self.equi_time = args.equi * picosecond
        self.time_step = args.time_step * femtosecond
        self.trajectory_interval = args.trajectory_interval * picosecond
        self.checkpoint_interval = args.checkpoint_interval * picosecond
        self.posres_k = args.posres # kilocalories_per_mole/angstrom**2
        self.equi_steps = int(self.equi_time / self.time_step + 0.5)
        self.prod_steps = int(self.prod_time / self.time_step + 0.5)
        self.trajectory_steps = int(self.trajectory_interval / self.time_step + 0.5)
        self.checkpoint_steps = int(self.checkpoint_interval / self.time_step + 0.5)
        self.trajectory_frames = self.prod_steps // self.trajectory_steps

        self.logger.info(f"{'time step':<20}= {str(self.time_step):>12}")
        self.logger.info(f"{'equilibration':<20}= {str(self.equi_time):>12} | {self.equi_steps:>9} steps")
        self.logger.info(f"{'production':<20}= {str(self.prod_time):>12} | {self.prod_steps:>9} steps")
        self.logger.info(f"{'checkpoint interval':<20}= {str(self.checkpoint_interval):>12} | " +
                    f"{self.checkpoint_steps:>9} steps")
        self.logger.info(f"{'trajectory interval':<20}= " +
                    f"{str(self.trajectory_interval.in_units_of(picosecond)):>12} | " +
                    f"{self.trajectory_steps:>9} steps")
        self.logger.info(f"{'trajectory frames':<20}= {str(self.trajectory_frames):>9}")
        self.logger.info(f"topol= {prmtop_path}")
        self.logger.info(f"coord= {inpcrd_path}")

        self.prmtop = AmberPrmtopFile(prmtop_path.as_posix())
        self.inpcrd = AmberInpcrdFile(inpcrd_path.as_posix())
        self.system = None
        self.integrator = None
        self.simulation = None
        self.barostat = None


    def check_solute(self) -> None:
        """check solute and terminate"""
        for atom in self.prmtop.topology.atoms():
            if not (atom.residue.name in ('HOH', 'Na+', 'Cl-') or atom.name.startswith('H')):
                self.logger(f"solute {atom}")
        sys.exit(0)


    def create_system(self) -> None:
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, 
                                               nonbondedCutoff=1.0*nanometer, 
                                               constraints=HBonds)


    def set_posres(self) -> None:
        self.logger.info(f"posres (k={self.posres_k} kcal/(A**2 mol)) for solute heavy atoms")
        force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addGlobalParameter("k", self.posres_k * kilocalories_per_mole/angstrom**2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for i, (atom, coor) in enumerate(zip(self.prmtop.topology.atoms(), self.inpcrd.positions)):
            if not (atom.residue.name in ('HOH', 'Na+', 'Cl-') or atom.name.startswith('H')): 
                # solute heavy atoms
                force.addParticle(i, coor.value_in_unit(nanometer))
        self.system.addForce(force)


    def set_integrator_and_simulation(self) -> None:
        # Using LangevinIntegator means NVT simulation
        self.integrator = LangevinIntegrator(self.temperature, 1 / picosecond, 2 * femtosecond)
        self.integrator.setConstraintTolerance(0.00001)
        self.simulation = Simulation(self.prmtop.topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.inpcrd.positions)
        if self.inpcrd.boxVectors is not None:
            self.simulation.context.setPeriodicBoxVectors(*self.inpcrd.boxVectors)

    def energy_min(self) -> None:
        self.logger.info('energy minimization')
        self.simulation.reporters.append(StateDataReporter(self.state_data_path.as_posix(),
                                                           10,step=True,
                                                           potentialEnergy=True,
                                                           temperature=True,
                                                           volume=True,
                                                           density=True))
        self.simulation.minimizeEnergy()
        self.simulation.reporters.pop()


    def nvt_warmup(self) -> None:
        self.logger.info(f'NVT equi. ({6000*self.time_step.in_units_of(picosecond)}) / warming up')
        dT = 5.0 * kelvin
        self.simulation.context.setVelocitiesToTemperature(dT)
        self.simulation.reporters.append(StateDataReporter(self.state_data_path.as_posix(), 
                                                           100,
                                                           step=True,
                                                           potentialEnergy=True,
                                                           temperature=True,
                                                           volume=True,
                                                           density=True))
        # 100 steps for every 5 degree increase
        for i in range(60): # 0..59
            T = min((i+1)*dT, self.temperature)
            self.integrator.setTemperature(T)
            self.simulation.step(100) 
        self.simulation.reporters.pop()

    def npt_releasing_posres(self) -> None:
        # Using Barostat means NPT simulation
        self.logger.info(f'NPT equi. ({self.equi_time.in_units_of(picosecond)}) / releasing posres')
        self.barostat = self.system.addForce(MonteCarloBarostat(self.pressure, self.temperature))
        self.simulation.context.reinitialize(True)
        self.simulation.reporters.append(StateDataReporter(self.state_data_path.as_posix(), 
                                                           1000,
                                                           step=True,
                                                           potentialEnergy=True,
                                                           temperature=True,
                                                           volume=True,
                                                           density=True))
        # 20 x 5 % decrease of POSRES force constant
        # end point should be zero
        dk = self.posres_k / 20.0
        for i in range(20): # 0..19
            k = min(self.posres_k - dk * (i+1), self.posres_k) * kilocalories_per_mole/angstrom**2
            self.simulation.context.setParameter('k', k)
            self.simulation.step(int(self.equi_steps/20))
        self.simulation.reporters.pop()

    def run(self) -> None:
        # Production
        self.logger.info(f"NPT prod. ({self.prod_time.in_units_of(nanosecond)})")
        self.simulation.reporters.append(CheckpointReporter(self.checkpoint_path.as_posix(), 
                                                            self.checkpoint_steps))
        self.simulation.reporters.append(DCDReporter(self.trajectory_path.as_posix(), 
                                                     self.trajectory_steps))
        self.simulation.reporters.append(StateDataReporter(self.state_data_path.as_posix(), 
                                                           self.trajectory_steps,
                                                           step=True,
                                                           potentialEnergy=True,
                                                           temperature=True,
                                                           volume=True,
                                                           density=True))
        self.simulation.step(self.prod_steps)
        self.logger.info(f"trajectory= {self.trajectory_path}")

        # create an empty file to mark completion
        self.logger.info("done")
        self.logger.removeHandler(self.logger_file_h)



def create_openmm_md_jobs(args: argparse.Namespace, logger: logging.Logger) -> list:
    """return list of md simulation jobs after checking input and output files"""
    priority_group = {
        'rec_solvated' : 1, 
        'complex_solvated' : 2,
        }

    jobs = []
    for filename in args.prmtop:
        prmtop_path = pathlib.Path(filename)
        # determine inpcrd filename based on prmtop filename
        if prmtop_path.name.endswith(".prmtop"):
            # jobname from 'mol-0167_solvated.prmtop' will be 'mol-0167_solvated'
            jobname = prmtop_path.name.replace(".prmtop", "")
            inpcrd_path = prmtop_path.parent / f"{jobname}.inpcrd"
        elif prmtop_path.name.endswith(".top"):
            jobname = prmtop_path.name.replace(".top", "")
            inpcrd_path = prmtop_path.parent / f"{jobname}.crd"
        
        priority = 3
        for k, v in priority_group.items():
            if k in jobname:
                priority = v

        # test whether input .prmtop and .inpcrd files exist
        if not prmtop_path.exists():
            logger.error(f"{prmtop_path} does not exist")
            sys.exit(1)
        if not inpcrd_path.exists():
            logger.error(f"{inpcrd_path} does not exist")
            sys.exit(1)
        
        outfile_prefix = f"openmm_md_job_{jobname}-out"
        # end_path is created when simulation is done
        trajectory_path = args.workdir / f'{outfile_prefix}.dcd'
        end_path        = args.workdir / f'{outfile_prefix}.end'

        # adding to jobs considering --overwrite option
        if trajectory_path.exists() and end_path.exists():
            if args.overwrite:
                logger.info(f"{outfile_prefix} exists / overwriting")
                jobs.append((prmtop_path, inpcrd_path, outfile_prefix, priority, jobname))
            else:
                logger.info(f"{outfile_prefix} exists / skipping")
                pass
        else:
            jobs.append((prmtop_path, inpcrd_path, outfile_prefix, priority, jobname))

    return sorted(jobs, key=itemgetter(3,4)) # sort by priority and jobname in ascending order
