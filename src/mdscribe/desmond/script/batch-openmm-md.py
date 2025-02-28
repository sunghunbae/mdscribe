#!/usr/bin/env python3

import os
import sys
import pathlib
import logging
import argparse

from mdscribe.omm import create_openmm_md_jobs, OpenMMSetup


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description='Run OpenMM Molecular Dynamics Simulation', 
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
    argparser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                           help="overwrite if output exists otherwise skip")
    argparser.add_argument("prmtop", nargs="+", help="Amber topology file(s) (.prmtop)")

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
