import sys
import csv
import os
import argparse
import pandas as pd

""" Analyze dihedral and pseudo-dihedrals """


try:
    from schrodinger.application.desmond.packages import analysis, topo, traj, traj_util
except ImportError:
    print("schrodinger python API is required")
    sys.exit(0)

# definition of dihedral angles

dihedral_def = {
    # standard protein backbone dihedral angles: phi, psi
    "phi"   : [(-1, "C"), (0, "N" ), (0, "CA"), (0, "C")],
    "psi"   : [( 0, "N"), (0, "CA"), (0, "C" ), (1, "N")],
    # standard nucleic acids backbone dihedral angles: alpha, beta, gamma, delta, epsilon, zeta
    "a"     : [(-1, "O3'"), (0, "P"  ), (0, "O5'"), (0, "C5'")],
    "b"     : [( 0, "P"  ), (0, "O5'"), (0, "C5'"), (0, "C4'")],
    "g"     : [( 0, "O5'"), (0, "C5'"), (0, "C4'"), (0, "C3'")],
    "d"     : [( 0, "C5'"), (0, "C4'"), (0, "C3'"), (0, "O3'")],
    "e"     : [( 0, "C4'"), (0, "C3'"), (0, "O3'"), (1, "P"  )],
    "z"     : [( 0, "C3'"), (0, "O3'"), (1, "P"  ), (1, "O5'")],
    # pseudo-dihedral angles
    "eta"   : [(-1, "C4'"), (0, "P"  ), (0, "C4'"), (1, "P"  )],
    "theta" : [( 0, "P"  ), (0, "C4'"), (1, "P"  ), (1, "C4'")],
    }



def get_dihedral_atom_ids(st, chain, resSeq, angle_name):
    if not angle_name in dihedral_def :
        return None

    residue_numbers = set()
    for (offset, atom_name) in dihedral_def[angle_name]:
        residue_numbers.add(offset+resSeq)
    atoms_ = []
    for a in st.atom:
        if (a.property['s_m_chain_name'] == chain) and \
            (int(a.property['i_m_residue_number']) in residue_numbers) :
            atoms_.append(a)

    atom_ids=[]
    for (offset, atom_name) in dihedral_def[angle_name]:
        for a in atoms_:
            if (int(a.property['i_m_residue_number']) == (resSeq+offset)) and \
                (a.property['s_m_pdb_atom_name'].strip() == atom_name) :
                atom_ids.append(a.index)
    
    if len(atom_ids) == 4 :
        return atom_ids
    
    return None


parser = argparse.ArgumentParser(description="trajectory dihedral angles",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--chain', dest="chain", default="A", help="chain to be calculated")
parser.add_argument('--pseudo', dest="pseudo", default=False, action="store_true" , help="pseudo dihedral angles")
parser.add_argument('--eta', dest="eta", default=False, action="store_true" , help="pseudo dihedral angle eta")
parser.add_argument('--theta', dest="theta", default=False, action="store_true" , help="pseudo dihedral angle theta")
parser.add_argument('--nucleic', dest="nucleic", default=False, action="store_true" , help="nucleic acid dihedral angles")
parser.add_argument('--alpha', dest="alpha", default=False, action="store_true" , help="nucleic acid dihedral alpha")
parser.add_argument('--beta', dest="beta", default=False, action="store_true" , help="nucleic acid dihedral beta")
parser.add_argument('--gamma', dest="gamma", default=False, action="store_true" , help="nucleic acid dihedral gamma")
parser.add_argument('--delta', dest="delta", default=False, action="store_true" , help="nucleic acid dihedral delta")
parser.add_argument('--epsilon', dest="epsilon", default=False, action="store_true" , help="nucleic acid dihedral epsilon")
parser.add_argument('--zeta', dest="zeta", default=False, action="store_true" , help="nucleic acid dihedral zeta")
parser.add_argument('--protein', dest="protein", default=False, action="store_true" , help="protein backbone dihedral angles")
parser.add_argument('--phi', dest="phi", default=False, action="store_true" , help="protein backbone dihedral phi")
parser.add_argument('--psi', dest="psi", default=False, action="store_true" , help="protein backbone dihedral psi")
parser.add_argument('--resSeq_start', dest="resSeq_start", type=int, default=0, help="first resSeq for analysis")
parser.add_argument('--resSeq_end', dest="resSeq_end", type=int, default=0, help="last resSeq for analysis")
parser.add_argument('--traj_start', dest="traj_start", type=float, default=0.0, help="traj_start fraction for analysis")
parser.add_argument('--traj_end', dest="traj_end", type=float, default=1.0, help="end fraction for analysis")
parser.add_argument('--csv', dest="csv", default="dihedrals.csv.gz", help="output csv or csv.gz filename")
parser.add_argument('--check', dest="check", default=False, action="store_true", help="just check atom index")
parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
args = parser.parse_args()


if __name__ == "__main__" :
    # dihedral angles to be calculated
    names = []
    if args.nucleic :
        names.extend(["a", "b", "g", "d", "e", "z"])
    if args.alpha:
        names.append("a")
    if args.beta:
        names.append("b")
    if args.gamma:
        names.append("g")
    if args.delta:
        names.append("d")
    if args.epsilon:
        names.append("e")
    if args.zeta:
        names.append("z")
    if args.pseudo:
        names.extend(["eta", "theta"])
    if args.eta:
        names.extend("eta")
    if args.theta:
        names.append("theta")
    if args.protein:
        names.extend(["phi", "psi"])
    if args.phi:
        names.append("phi")
    if args.psi:
        names.append("psi")

    data = {"name":[], "resid":[], "angle_name":[], "angle":[]}
    total_cms_files = len(args.cms)
    count_cms_files = 0
    
    print("dihedral angles      : ", ",".join(names))
    if not args.check:
        print("trajectory fractions : ", args.traj_start, "-", args.traj_end) 
    
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
            startframe = max(int(num_frames * args.traj_start), 0)
            endframe   = min(int(num_frames * args.traj_end  ), num_frames)
        
            sys.stdout.write(f"[{count_cms_files}/{total_cms_files}] ({startframe}-{endframe})/{num_frames} {input_cms:60s}\r")
            sys.stdout.flush()

        # residue numbers
        list_resSeq = []
        list_resName = []
        for c in cms_model.chain:
            if c.name == args.chain:
                for r in c.residue:
                    first_atom = list(r.atom)[0]
                    list_resSeq.append(int(first_atom.property['i_m_residue_number']))
                    list_resName.append(first_atom.property['s_m_pdb_residue_name'])
                # resNames = [ r.pdbres.strip() for r in c.residue ] # residue name
        
        list_resSeq = sorted(list_resSeq)
        resSeq_first = max(args.resSeq_start, list_resSeq[0])
        resSeq_last = min(args.resSeq_end, list_resSeq[-1])
        list_resSeq = [ x for x in list_resSeq if x >= resSeq_first and x <= resSeq_last ]

        list_of_obj = [] # dihedrals to analyze
        list_of_tag = []

        for angle_name in names:
            for resSeq in list_resSeq:
                dihed_aids = get_dihedral_atom_ids(cms_model, args.chain, resSeq, angle_name)
                # returns None unless four atom indexes are determined

                if not dihed_aids:
                    continue

                if args.check:
                    print(f"angle= {angle_name:8s} res= {resSeq:3d} aids= {dihed_aids}")
                    continue
                
                # convert atom index to global index
                #dihed_gids = topo.aids2gids(cms_model, dihed_aids, include_pseudoatoms=False)

                """
                print("angle_name=", angle_name, "resSeq=", resSeq)
                print("dihed_aids=", dihed_aids)
                print("    atom=", [cms_model.atom[x].property['s_m_pdb_atom_name'].strip() for x in dihed_aids])
                print("    res=", [cms_model.atom[x].property['i_m_residue_number'] for x in dihed_aids])
                print("dihed_gids=", dihed_gids)
                print("    atom=", [msys_model.atom(x).mass for x in dihed_gids])
                """
                
                # dihedral angle
                # * = unpacking argument list
                obj = analysis.Torsion(msys_model, cms_model, *dihed_aids)
                #print("frame_torsions=")
                #print([float(f"{v:.1f}") for v in analysis.analyze(trj[startframe:endframe], obj)])

                list_of_obj.append(obj)
                list_of_tag.append((angle_name, resSeq))
                #print()
                
        if not args.check:

            list_of_results = analysis.analyze(trj[startframe:endframe], *list_of_obj) 
            #print("list_of_results=", list_of_results[:20])
            # analysis.analyze returns a list of list

            for (angle_name, resSeq), results in zip(list_of_tag, list_of_results):
                #angle_in_degree = [v+360. if v < 0 else v for v in results] # 0...360
                angle_in_degree = [(v+360.0) % 360 for v in results] # 0...360
                num_data = len(angle_in_degree)
                data["name"]        += num_data * [name,]
                data["resid"]       += num_data * [resSeq,]
                data["angle_name"]  += num_data * [angle_name,]
                data["angle"]       += angle_in_degree


    if not args.check:
        df = pd.DataFrame(data)
        df.to_csv(args.csv, sep=",", index=False, header=True, float_format="%.1f")
        print()
        print(f"dataframe written to {args.csv}")
        print()
