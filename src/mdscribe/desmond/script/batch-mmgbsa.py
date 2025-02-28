#!/usr/bin/env python3

import pathlib
import argparse
import sys
import logging
import concurrent.futures
import pandas as pd

from mdscribe.amber import check_ambertools, run_mmpbsa_py


argparser = argparse.ArgumentParser(description='run MMPBSA.py in parallel',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    epilog='''B. R. Miller 3rd, et al., 
                                    MMPBSA.py: An Efficient Program for End-State Free Energy Calculations.
                                    J. Chem. Theory Comput. 8, 3314–3321 (2012).''')
argparser.add_argument("-m", "--method", dest="method", default=1, type=int, choices=[1,3], help="1 or 3 trajectories?")
argparser.add_argument("-s", "--start-frame", dest="start_frame", default=1, type=int, help="start frame")
argparser.add_argument("-e", "--end-frame", dest="end_frame", default=100, type=int, help="end frame")
argparser.add_argument("-i", "--interval", dest="interval", default=1, type=int, help="interval")
argparser.add_argument("--igb", dest="igb", default=2, type=int, help="igb")
argparser.add_argument("--salt", dest="salt", default=0.15, type=float, help="salt conc. (M)")
argparser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true", help="overwrite or skip if output exists")
argparser.add_argument("--receptor-name", dest="receptor_name", default="rec", help="receptor name")
argparser.add_argument("--workdir", dest="workdir", default=".", help="working directory")
argparser.add_argument("--max-workers", dest="max_workers", default=10, type=int, help="max concurrent workers")
argparser.add_argument("--dcd-prefix", dest="dcd_prefix", default="openmm_md_job", help="prefix for .dcd files")
argparser.add_argument("--dcd-suffix", dest="dcd_suffix", default="solvated-out", help="suffix for .dcd files")
argparser.add_argument("-p", "--prmtop-dir", dest="prmtop_dir", required=True, help="directory for the Amber .prmtop files")
argparser.add_argument("complex_dcd", nargs="+", help="complex md trajectory file(s) (.dcd)")

# argparse run into recursion error when a number of arguments exceed some limit
if len(sys.argv) > 20:
    args = argparser.parse_args(sys.argv[1:20])
    args.complex_dcd += sys.argv[20:]
else:
    args = argparser.parse_args()

args.workdir = pathlib.Path(args.workdir).resolve()
args.prmtop_dir = pathlib.Path(args.prmtop_dir).resolve()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # level: DEBUG < INFO < WARNING < ERROR < CRITICAL
logger_formatter = logging.Formatter(fmt='%(asctime)s %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
logger_file_h = logging.FileHandler(f"amber-mmgbsa.log", mode="a", encoding='utf-8')
logger_file_h.setFormatter(logger_formatter)
logger.addHandler(logger_file_h)

check_ambertools(logger)

logger.info(f"{'method':<20} {args.method}")
logger.info(f"{'start frame':<20} {args.start_frame}")
logger.info(f"{'end frame':<20} {args.end_frame}")
logger.info(f"{'interval':<20} {args.interval}")
logger.info(f"{'igb':<20} {args.igb}")
logger.info(f"{'salt':<20} {args.salt}")

data = {"ligand":[], "dG":[], "stdev":[], "sem":[], "filename":[]}

complex_dcd_paths = [ pathlib.Path(p).resolve() for p in args.complex_dcd ]
with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
    futures = {
        executor.submit(run_mmpbsa_py, args, i, dcd): dcd for i, dcd in enumerate(complex_dcd_paths)
    }
    for future in concurrent.futures.as_completed(futures):
        result = future.result()
        if result is None:
            continue
        data["ligand"].append(result.ligand)
        data["dG"].append(result.dG)
        data["stdev"].append(result.stdev)
        data["sem"].append(result.sem)
        data["filename"].append(result.output_path.name)
        logger.info(f"{result.output_path.name:<20} dG= {result.dG:9.4f} +/- {result.stdev:9.4f} (SEM= {result.sem:9.4f})")

df = pd.DataFrame(data).sort_values(by="dG", ascending=True)
df.to_csv("amber-mmgbsa.csv", index=False, float_format='%.4f')