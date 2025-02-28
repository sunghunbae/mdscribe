#!/usr/bin/env python3

import os
import sys
import pathlib
import logging
import argparse

from mdscribe.omm import create_openmm_md_jobs, OpenMMSetup
from mdscribe.amber import run_mmpbsa_py

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description='run OpenMM & MMPBSA.py sequentially', 
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("--gpu", dest="gpu", default=0, type=int, choices=[0,1], help="GPU device to run md")
    argparser.add_argument("-t", "--prod", dest="prod", default=1.0, type=float, 
                           help="production simulation time (ns)")
    argparser.add_argument("-i", "--trajectory-interval", dest="trajectory_interval", default=10.0, type=float, 
                           help="trajectory interval (ps)")
    argparser.add_argument("--equi", dest="equi", default=100.0, type=float, help="NPT equilibration time (ps)")
    argparser.add_argument("--time-step", dest="time_step", default=4.0, type=float, help="time step (fs)")
    argparser.add_argument("--checkpoint-interval", dest="checkpoint_interval", default=200.0, type=float, 
                           help="checkpoint interval (ps)")
    argparser.add_argument("--posres", dest="posres", default=50.0, type=float, 
                           help="force constant for positional restraints (kcal/mol/A^2)")
    argparser.add_argument("--temperature", dest="temperature", default=300.0, type=float, help="temperature (K)")
    argparser.add_argument("--pressure", dest="pressure", default=1.0, type=float, help="pressure (Bar)")
    argparser.add_argument("--workdir", dest="workdir", default=".", help="working directory")
    argparser.add_argument("--dcd-prefix", dest="dcd_prefix", default="openmm_md_job", help="prefix for .dcd files")
    argparser.add_argument("--dcd-suffix", dest="dcd_suffix", default="solvated-out", help="suffix for .dcd files")
    argparser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true", 
                           help="overwrite or skip if output exists")
    argparser.add_argument("--method", dest="method", default=1, type=int, choices=[1, 3], 
                           help="GBSA 1 or 3 trajectories?")
    argparser.add_argument("--start-frame", dest="start_frame", default=1, type=int, help="GBSA start frame")
    argparser.add_argument("--end-frame", dest="end_frame", default=100, type=int, help="GBSA end frame")
    argparser.add_argument("--interval", dest="interval", default=1, type=int, help="GBSA interval")
    argparser.add_argument("--igb", dest="igb", default=2, type=int, help="GBSA igb")
    argparser.add_argument("--salt", dest="salt", default=0.15, type=float, help="GBSA salt conc. (M)")
    argparser.add_argument("--receptor-name", dest="receptor_name", default="rec", help="GBSA receptor name")
    argparser.add_argument("--prmtop-dir", dest="prmtop_dir", help="GBSA directory for the Amber .prmtop files" 
                           "(to be set by the .prmtop files automatically)")
    
    argparser.add_argument("prmtop", nargs="+", help="Amber topology(.prmtop) file(s)")

    # argparse run into recursion error when a number of arguments exceed some limit
    if len(sys.argv) > 20:
        args = argparser.parse_args(sys.argv[1:20])
        args.prmtop += sys.argv[20:]
    else:
        args = argparser.parse_args()

    args.workdir = pathlib.Path(args.workdir)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO) # level: DEBUG < INFO < WARNING < ERROR < CRITICAL
    logger_formatter = logging.Formatter(fmt='%(asctime)s %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    logger_stream_h = logging.StreamHandler(stream=sys.stdout)
    logger_stream_h.setFormatter(logger_formatter)
    logger.addHandler(logger_stream_h)
    logger.info(f'GPU device=', args.gpu)

    os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu)

    jobs = create_openmm_md_jobs(args, logger)
    for (prmtop_path, inpcrd_path, outfile_prefix, priority, jobname) in jobs:
        omm = OpenMMSetup(args, prmtop_path, inpcrd_path, outfile_prefix, logger)
        omm.create_system()
        omm.set_posres()
        omm.set_integrator_and_simulation()
        omm.energy_min()
        omm.nvt_warmup()
        omm.npt_releasing_posres()
        omm.run()
        if "complex" in jobname: # ex. mol-167-complex
            args.prmtop_dir = prmtop_path.parent.as_posix()
            result = run_mmpbsa_py(args, 0, omm.trajectory_path.resolve())
            logger.info(f"{result.output_path.name.as_posix():<20} "
                        f"dG= {result.dG:9.4f} +/- {result.stdev:9.4f} "
                        f"(SEM= {result.sem:9.4f})")
