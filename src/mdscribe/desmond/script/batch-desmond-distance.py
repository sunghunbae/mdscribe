import sys
import csv
import os
import argparse
import pandas as pd

try:
    from schrodinger.application.desmond.packages import analysis, topo, traj, traj_util
    from schrodinger.structutils import analyze
    from schrodinger import structure
except ImportError:
    print("Schrodinger Python API is required")
    sys.exit(0)



""" Analyze distance 

    The .cms files only contain physical atoms whereas the trajectory files may contain additional pseudo atoms. 
    This leads to two atom indices:
        Atom ID (AID) for structural files such as cms and mae files
            physical atoms only
            1-indexed
            results of ASL/SMARTS evaluation

    Global ID (GID) for trajectory files and Desmond MSYS files
            all atoms, both physical and pseudo
            0-indexed

Usage:

$SCHRODINGER/run ~/bucket/schrodinger/batch-desmond-distance.py \
    --receptor "C-555-SG" \
    --ligand "smarts. [CH3]"  \
    r0?/desmond_md_job_*-out.cms
"""

def get_atom_id(st:structure.Structure, chain: str, resSeq: int, name: str) -> int:
    """Get an aid with chain/resSeq/name.

    Args:
        st (structure.Structure): structure object
        chain (str): chainId
        resSeq (int): residue number
        name (str): atom name

    Returns:
        int: atom index
    """
    for a in st.atom:
        if (a.property['s_m_chain_name'] == chain) and \
            (a.property['i_m_residue_number'] == resSeq) and \
            (a.property['s_m_pdb_atom_name'].strip() == name) :
            return a.index
    return None



parser = argparse.ArgumentParser(
    description="trajectory distance",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--atoms', dest="atoms", action="append", help="pair(s) of chain-resSeq-name separated by '-'")
parser.add_argument('--receptor', dest="receptor", default="", help="receptor chain-resSeq-name separated by '-'")
parser.add_argument('--ligand', dest="ligand", default="", help="ligand ASL")
parser.add_argument('--start', dest="start", type=float, default=0.0, help="start fraction for analysis")
parser.add_argument('--end', dest="end", type=float, default=1.0, help="end fraction for analysis")
parser.add_argument('--n', dest="n", type=int, default=0, help="target sample points")
parser.add_argument('--csv', dest="csv", default="distance.csv.gz", help="output csv or csv.gz filename")
parser.add_argument('--check', dest="check", default=False, action="store_true", help="just check atom index")
parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
args = parser.parse_args()


if __name__ == "__main__" :

    # check give atom representations
    distance_atom_pairs = []
    distance_names = []
    if not (args.receptor and args.ligand):
        if len(args.atoms) >= 2 and len(args.atoms) % 2 == 0:
            for i in range(0, len(args.atoms)//2):
                try:
                    (chain_1, resSeq_1, atomname_1) =  args.atoms[2*i+0].strip().split('-')
                    resSeq_1 = int(resSeq_1)
                except:
                    print("Error:", args.atoms[2*i+0])
                    print("each atom must be defined by {chain}-{reqSeq}-{atom_name} separated by '-'")
                    sys.exit(0)
                try:
                    (chain_2, resSeq_2, atomname_2) =  args.atoms[2*i+1].strip().split('-')
                    resSeq_2 = int(resSeq_2)
                except:
                    print("Error:", args.atoms[2*i+1])
                    print("each atom must be defined by {chain}-{reqSeq}-{atom_name} separated by '-'")
                    sys.exit(0)
                distance_name = args.atoms[2*i+0] + ':' + args.atoms[2*i+1]
                distance_names.append(distance_name)
                distance_atom_pairs.append((chain_1, resSeq_1, atomname_1, chain_2, resSeq_2, atomname_2, distance_name))

        else:
            print("Error: N x 2 atoms (or residues) are required to define distance(s)")
            sys.exit(0)

    else:
        # define distance_names
        distance_names = [ args.receptor + ':' + args.ligand ]

    data = {"name":[], "time":[], "atom_pair":[], "distance":[]}
    total_cms_files = len(args.cms)
    count_cms_files = 0
    
    print("distance name(s): ")
    for i, distance_name in enumerate(distance_names, start=1):
        print(f"    [{i}] {distance_name}")

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

        
        list_of_obj = [] # distances to analyze
        list_of_tag = []

        if distance_atom_pairs:
            for (chain_1, resSeq_1, atomname_1, chain_2, resSeq_2, atomname_2, distance_name) in distance_atom_pairs:
                    aid_1 = get_atom_id(cms_model, chain_1, resSeq_1, atomname_1)
                    aid_2 = get_atom_id(cms_model, chain_2, resSeq_2, atomname_2)
                    # returns None unless four atom indexes are determined
                    if not (aid_1 and aid_2):
                        continue
                    if args.check:
                        print(input_cms)
                        print(f"    distance= {distance_name} aids= {aid_1} {aid_2}")
                        continue
                    obj = analysis.Distance(msys_model, cms_model, aid_1, aid_2)
                    list_of_obj.append(obj)
                    list_of_tag.append(distance_name)
        else:
            (rec_chain, rec_resSeq, rec_atomname) =  args.receptor.strip().split('-')
            rec_resSeq = int(rec_resSeq)
            aid_1 = get_atom_id(cms_model, rec_chain, rec_resSeq, rec_atomname)
            ligands = analyze.find_ligands(cms_model)
            aid_2 = analyze.evaluate_asl(cms_model, args.ligand)
            # specify substructure of the ligand
            aid_2 = list(set(aid_2) & set(ligands[0].atom_indexes))[0]
            # define list_of_obj, list_of_tag
            if args.check:
                print(input_cms)
                print(f"    distance= {distance_name} aids= {aid_1} {aid_2}")
                continue
            obj = analysis.Distance(msys_model, cms_model, aid_1, aid_2)
            list_of_obj.append(obj)
            list_of_tag.append(distance_name)
                
        if not args.check:
            
            list_of_results = analysis.analyze(trj[startframe:endframe:step], *list_of_obj)
            #list_of_results = analysis.analyze(trj, *list_of_obj)

            # analysis.analyze would return a list when one distance is calculated
            # analysis.analyze would return a list of list for more than one distances
            # for uniform process, let's make a list of list for one distance

            if len(list_of_tag) == 1:
                list_of_results = [ list_of_results ]

            for distance_name, results in zip(list_of_tag, list_of_results):
                num_data = len(results)
                data["name"] += num_data * [name,]
                data["time"] += [ fr.time for fr in trj[startframe:endframe:step] ] # (ps)
                data["atom_pair"] += num_data * [distance_name,]
                data["distance"] += results


    if not args.check:
        df = pd.DataFrame(data)
        df.to_csv(args.csv, sep=",", index=False, header=True, float_format="%.1f")
        print()
        print(f"dataframe written to {args.csv}")
        print()
