#!/usr/bin/env python3

import pathlib
import argparse
import sys
import concurrent.futures
import logging

from mdscribe.amber import check_ambertools, run_amber_setup

"""
Protein Force Fields
For a full description of the force fields and how to load them, 
please read Chapter 3 of the Amber Reference Manual.

ff19SB (recommended)
The ff19SB force field is the most current protein force field.
The new ff19SB forcefield has shown to improve amino-acid dependent properties 
such as helical propensities and reproduces the differences in amino acid-specific Ramachandram Map 
using amino acid specific CMAPS.

The ff19SB force field contains amino-acid specific backbone parameters. 
ff19SB pairs best with the more accurate OPC water model. 
Note: This adds up to 33% computational time over the ff14SBonlysc/OPC3 and ff14SB/TIP3P combinations 
but is expected to be more accurate.
Ref: C. Tian; K. Kasavajhala; K. Belfon; L. Raguette; H. Huang; A. Migues; J. Bickel; 
Y. Wang; J. Pincay; Q. Wu; C. Simmerling. ff19SB: Amino-Acid-Specific Protein Backbone 
Parameters Trained against Quantum Mechanics Energy Surfaces in Solution. 
J. Chem. Theory Comput., 2020, 16, 528–552.

ff14SBonlysc
In the ff14SBonlysc force field, the side-chain dihedral parameters were 
fit to quantum mechanical data for each amino acid. 
ff14SBonlysc is the same model as ff14SB but without the empirical backbone corrections for TIP3P.
For simulations in implicit solvent (igb=8), ff14SBonlysc is best.

If you want a 3-point water model to minimize computational cost, this is a good choice

ff14SB
The ff14SB force field is intended for use with the TIP3P water model. 
Backbone parameters were based upon alanine and glycine, 
including a TIP3P-specific correction to the backbone parameters.

referece - https://ambermd.org/AmberModels_proteins.php
"""

    
argparser = argparse.ArgumentParser(description='Setup Amber prmtop and inpcrd files in batch mode',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
argparser.add_argument("--ff-protein", dest="ff_protein", default="ff14SB", help="force field for protein")
argparser.add_argument("--ff-water", dest="ff_water", default="tip3p", help="force field for water")
argparser.add_argument("--ff-ligand", dest="ff_ligand", default="gaff2", help="force field for ligand")
argparser.add_argument("--charge", dest="charge", default="mol2",
                        help="charge method - mol2[from .mol2 file] gas[Gasteiger] bcc[AM1-BCC] mul[Mulliken;AM1]")
argparser.add_argument("--water", dest="water", default=15.0, help="water box buffer (A)")
argparser.add_argument("--salt", dest="salt", default=0.15, help="salt concentration (M)")
argparser.add_argument("--gbsa", dest="gbsa", default=False, action="store_true", help="setup for GBSA")
argparser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true", help="overwrite or skip if output exists")
argparser.add_argument("--max-workers", dest="max_workers", default=10, type=int, help="max concurrent workers")
argparser.add_argument("--workdir", dest="workdir", default=".", help="working directory")
argparser.add_argument("--receptor-name", dest="receptor_name", default="rec", help="receptor name")
argparser.add_argument("--receptor-pdb", "-r", dest="receptor_pdb", help="receptor PDB or PDBQT")
argparser.add_argument("--no-leap", dest="leap", default=True, action="store_false", help="do not run leap")
argparser.add_argument("ligand_mol2", nargs="+", help="ligand MOL2")

# argparse run into recursion error when a number of arguments exceed some limit
if len(sys.argv) > 20:
    args = argparser.parse_args(sys.argv[1:20])
    args.ligand_mol2 += sys.argv[20:]
else:
    args = argparser.parse_args()

ligand_mol2_paths = []

if args.receptor_pdb:
    args.receptor_pdb = pathlib.Path(args.receptor_pdb).resolve()
    ligand_mol2_paths.append(None) # receptor

args.workdir = pathlib.Path(args.workdir).resolve()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # level: DEBUG < INFO < WARNING < ERROR < CRITICAL
logger_formatter = logging.Formatter(fmt='%(asctime)s %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
logger_file_h = logging.FileHandler(f"amber-setup.log", mode="a", encoding='utf-8')
logger_file_h.setFormatter(logger_formatter)
logger.addHandler(logger_file_h)

check_ambertools(logger)

logger.info(f"{'protein force field':<20}: {args.ff_protein}")
logger.info(f"{'water force field':<20}: {args.ff_water}")
logger.info(f"{'ligand force field':<20}: {args.ff_ligand}")
logger.info(f"{'charge method':<20}: {args.charge}")
logger.info(f"{'water box buffer (A)':<20}: {args.water}")
logger.info(f"{'salt conc. (M)':<20}: {args.salt}")
logger.info(f"{'GBSA radii':<20}: {args.gbsa}")
logger.info(f"{'receptor PDB':<20}: {args.receptor_pdb}")
logger.info(f"{'receptor name':<20}: {args.receptor_name}")
logger.info(f"{'leap':<20}: {args.leap}")
logger.info(f"generating Amber files (max_workers={args.max_workers})")

# get the full path of the ligand .mol2

ligand_mol2_paths += [ pathlib.Path(p).resolve() for p in args.ligand_mol2 ]

n = min(args.max_workers, len(ligand_mol2_paths))
with concurrent.futures.ThreadPoolExecutor(max_workers=n) as executor:
    futures = {
        executor.submit(run_amber_setup, args, i, p): p for i, p in enumerate(ligand_mol2_paths)
    }
    for future in concurrent.futures.as_completed(futures):
        ligand_mol2_path = futures[future]
        output_paths = future.result()
        logger.info(f"{ligand_mol2_path}:")
        for p in sorted(output_paths):
            logger.info(f"    {p}")