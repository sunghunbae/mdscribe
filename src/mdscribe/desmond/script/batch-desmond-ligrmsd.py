""" Analyze ligand pose stability """

try:
    from schrodinger.application.desmond.packages import analysis, topo, traj_util
except ImportError:
    print("schrodinger python API is required")
    sys.exit(0)

import sys
import os
import argparse
import pandas as pd
import numpy as np

# ASL ex. not (chain. D) and (not (atom.ele H)) 

parser = argparse.ArgumentParser(description="trajectory dihedral angles",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--receptor', dest="receptor", default="mol. 1", 
    help="receptor ASL")
parser.add_argument('--ligand', dest="ligand", default="mol. 2", 
    help="ligand ASL")
parser.add_argument('--threshold', dest="threshold", type=float, default=2.0, 
    help="rmsd threshold")
parser.add_argument('--out', dest="out", default="ligand-rmsd.csv.gz", 
    help="output csv or csv.gz filename")
parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
args = parser.parse_args()

if __name__ == "__main__" :

    data = {
        "name" : [], 
        f"rmsd_percent({args.threshold:.1f}A)" : [], 
        "rmsd_median" : [], 
        "rmsd_mean" : [], 
        "rmsd_max" : [],
        "dist_init" : [], 
        "dist_median" : [], 
        "dist_mean" : [], 
        "dist_max" : [],
    }
    total_cms_files = len(args.cms)
    count_cms_files = 0

    # by default, heavy atom and no water
    receptor_asl = f"({args.receptor}) and (not atom.ele H) and (not water)"
    ligand_asl = f"({args.ligand}) and (not atom.ele H) and (not water)"

    logfile = args.out
    if logfile.endswith(".gz"):
        logfile = logfile[:-3]
    if logfile.endswith(".csv"):
        logfile = logfile[:-4]
    logfile = logfile + ".log"
    logging = open(logfile,"w")

    for input_cms in args.cms:

        (msys_model, cms_model, traj) = traj_util.read_cms_and_traj(input_cms)
        num_frames = len(traj)
        count_cms_files += 1

        name, ext_ = os.path.splitext(os.path.basename(input_cms))
        name = name.replace("desmond_md_job_","").replace("-out","")

        print(f"[{count_cms_files}/{total_cms_files}] {input_cms} ({num_frames} frames)")
        logging.write(f"[{count_cms_files}/{total_cms_files}] {input_cms} ({num_frames} frames)\n")

        receptor_aids = cms_model.select_atom(receptor_asl)
        print(f"\treceptor {len(receptor_aids)} atoms : {receptor_asl}")
        logging.write(f"\treceptor {len(receptor_aids)} atoms : {receptor_asl}\n")

        receptor_gids = topo.aids2gids(cms_model, receptor_aids, include_pseudoatoms=False)
        receptor_ref_pos = traj[0].pos(receptor_gids)

        ligand_aids = cms_model.select_atom(ligand_asl)
        print(f"\tligand {len(ligand_aids)} atoms : {ligand_asl}")
        logging.write(f"\tligand {len(ligand_aids)} atoms : {ligand_asl}\n")

        ligand_gids = topo.aids2gids(cms_model, ligand_aids, include_pseudoatoms=False)
        ligand_ref_pos = traj[0].pos(ligand_gids)
        
        # Ligand Root Mean Square Deviation from reference positions, 
        # with optional alignment fitting. 
        # Taking conformational symmetry into account.
        ligand_rmsd_analyzer = analysis.LigandRMSD(msys_model, cms_model, 
            aids = ligand_aids, 
            ref_pos = ligand_ref_pos,
            fit_aids = receptor_aids,
            fit_ref_pos = receptor_ref_pos,
            )

        receptor_centroid = analysis.Centroid(msys_model, cms_model, asl=receptor_asl)
        ligand_centroid = analysis.Centroid(msys_model, cms_model, asl=ligand_asl)
        ligand_dist_analyzer = analysis.Distance(msys_model, cms_model, 
            xid0 = receptor_centroid, 
            xid1 = ligand_centroid,
            )

        results = analysis.analyze(traj, ligand_rmsd_analyzer, ligand_dist_analyzer)
        # results[0] : list of rmsd
        # results[1] : list of centroid-to-centroid distance 
        
        rmsd_within_threshold = np.sum([1 for v in results[0] if v < args.threshold])
        rmsd_percent = 100.0 * rmsd_within_threshold / num_frames
        rmsd_mean_ = np.mean(results[0])
        rmsd_median_ = np.median(results[0])
        rmsd_max_ = np.max(results[0])
        dist_init_ = results[1][0]
        dist_mean_ = np.mean(results[1])
        dist_median_ = np.median(results[1])
        dist_max_ = np.max(results[1])
        print(f"\tLigandRMSD {rmsd_within_threshold} out of {num_frames} within {args.threshold:.1f}A ({rmsd_percent:.1f} %) ",end="")
        print(f"median {rmsd_median_:.2f} mean {rmsd_mean_:.2f} max {rmsd_max_:.2f}")
        print(f"\tCentroid Distance init {dist_init_:.2f} median {dist_median_:.2f} mean {dist_mean_:.2f} max {dist_max_:.2f}")
        print()
        logging.write(f"\tLigandRMSD {rmsd_within_threshold} out of {num_frames} within {args.threshold:.1f}A ({rmsd_percent:.1f} %) ")
        logging.write(f"median {rmsd_median_:.2f} mean {rmsd_mean_:.2f} max {rmsd_max_:.2f}\n")
        logging.write(f"\tCentroid Distance init {dist_init_:.2f} median {dist_median_:.2f} mean {dist_mean_:.2f} max {dist_max_:.2f}\n\n")
        
        data["name"].append(name)
        data[f"rmsd_percent({args.threshold:.1f}A)"].append(rmsd_percent)
        data["rmsd_median"].append(rmsd_median_)
        data["rmsd_mean"].append(rmsd_mean_)
        data["rmsd_max"].append(rmsd_max_)
        data["dist_init"].append(dist_init_)
        data["dist_median"].append(dist_median_)
        data["dist_mean"].append(dist_mean_)
        data["dist_max"].append(dist_max_)
    
    df = pd.DataFrame(data)
    df.to_csv(args.out, sep=",", index=False, header=True, float_format="%.3f")
    print(f"dataframe written to {args.out}")
    print()
    logging.write(f"dataframe written to {args.out}\n\n")
