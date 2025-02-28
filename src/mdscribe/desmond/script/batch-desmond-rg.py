import sys
import csv
import os
import argparse
import pandas as pd
import numpy as np

""" Analyze radius of gyration """


try:
    from schrodinger.application.desmond.packages import analysis, topo, traj, traj_util
    from schrodinger.structutils import analyze
except ImportError:
    print("schrodinger python API is required")
    sys.exit(0)


def check_asl_atoms(st, asl):
    # check ASL
    aids = analyze.evaluate_asl(st, asl)
    # collect chain, resSeq, resName, atomnames matching the ASL
    chain_ = set()
    resSeq_ = set()
    resName_ = set()
    # atomName_ = set()
    for a in st.atom:
        if a.index in aids:
            chain_.add(a.property['s_m_chain_name'])
            resSeq_.add(str(a.property['i_m_residue_number']))
            resName_.add(a.property['s_m_pdb_residue_name'].strip())
            # atomName_.add(a.property['s_m_pdb_atom_name'].strip())
    print("  asl: ", asl, "(number of atoms= ", len(aids),")")
    print("  chain(s): " + " ".join(chain_))
    print("  residue number(s): " + " ".join(resSeq_))
    print("  residue name(s): " + " ".join(resName_))
    # print("atom name(s): " + " ".join(atomName_))


parser = argparse.ArgumentParser(description="trajectory Rg",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--asl', dest="asl", action="append", help="atom selection language(s)")
parser.add_argument('--start', dest="start", type=float, default=0.0, help="start fraction for analysis")
parser.add_argument('--end', dest="end", type=float, default=1.0, help="end fraction for analysis")
parser.add_argument('--n', dest="n", type=int, default=0, help="target sample points")
parser.add_argument('--csv', dest="csv", default="rg.csv.gz", help="output csv or csv.gz filename")
parser.add_argument('--check', dest="check", default=False, action="store_true", help="just check atom index")
parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
args = parser.parse_args()


if __name__ == "__main__" :
    data = {"name":[], "asl":[], "Rg":[]}
    total_cms_files = len(args.cms)
    count_cms_files = 0

    if not args.check:
        print("trajectory fractions : ", args.start, "-", args.end) 
    
    for input_cms in args.cms:
        # id
        name = os.path.basename(input_cms).replace("-out.cms","").replace("desmond_md_job_","")
        count_cms_files += 1

        # read .cms and its associated trajectory file
        # <jobname>-out.cms
        # <jobname>_trj

        if args.check:
            # only use cms
            (msys_model, cms_model) = topo.read_cms(input_cms)
        else:
            try:
                (msys_model, cms_model, trj) = traj_util.read_cms_and_traj(input_cms)
            except:
                (msys_model, cms_model) = topo.read_cms(input_cms)
                try:
                    trj_path = topo.find_traj_path(cms_model)  # or
                except:
                    try:
                        trj_path = topo.find_traj_path_from_cms_path(input_cms)
                    except:
                        trj_path = input_cms.replace("-out.cms","_trj")
                trj = traj.read_traj(trj_path)
        
            # number of frames
            num_frames = len(trj)
            startframe = max(int(num_frames * args.start), 0)
            endframe   = min(int(num_frames * args.end  ), num_frames)
            if args.n > 1:
                step = (endframe - startframe + 1) // args.n
            else:
                step = 1
            sys.stdout.write(f"[{count_cms_files}/{total_cms_files}] ({startframe}-{endframe})/{num_frames} {input_cms:60s}\r")
            sys.stdout.flush()

        list_of_obj = []
        list_of_tag = []
        for asl in args.asl:
            if args.check:
                print(input_cms)
                check_asl_atoms(cms_model, asl)
                continue
            obj = analysis.Gyradius(msys_model, cms_model, asl=asl)
            list_of_obj.append(obj)
            list_of_tag.append(asl)

        if not args.check:
            list_of_results = analysis.analyze(trj[startframe:endframe:step], *list_of_obj)
            # analysis.analyze would return a list when one distance is calculated
            # analysis.analyze would return a list of list for more than one distances
            # for uniform process, let's make a list of list for one distance
            if len(list_of_tag) == 1:
                list_of_results = [ list_of_results ]

            for tag, results in zip(list_of_tag, list_of_results):
                num_data = len(results)
                data["name"]    += num_data * [name,]
                data["asl"]     += num_data * [tag,]
                data["Rg"]      += results


    if not args.check:
        df = pd.DataFrame(data)
        df.to_csv(args.csv, sep=",", index=False, header=True, float_format="%.1f")
        print()
        print(f"dataframe written to {args.csv}")
        print()
