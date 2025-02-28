import os
import sys
import argparse
import subprocess
import re
import csv
import os
import argparse
import random
import glob
import shutil

import numpy as np
import pandas as pd

from typing import List

run = "{}/run".format(os.environ["SCHRODINGER"])
multisim = "{}/utilities/multisim".format(os.environ["SCHRODINGER"])

try:
    from schrodinger import structure
    from schrodinger.structutils import analyze
    from schrodinger.application.desmond.packages import analysis, topo, traj, traj_util
except ImportError:
    print("Schrodinger Python API is required")
    sys.exit(0)


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


def get_dihedral_atom_indices(st, chain, resSeq, angle_name) -> List[int] | None:
    """Get 4 atom indices for a dihedral angle.
    """
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


def check_asl_atoms(st, asl:str) -> None:
    """Check ASL matching atoms.

    Args:
        st (_type_): _description_
        asl (_type_): _description_
    """
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



def batch_aslcheck() -> None:
    """Check ASL expression.

    entry > molecule > chain > residue > atom

    Examples:
        >>> 'entry. e1'
        >>> 'entry.name recep, lig*'
        
        >>> 'mol. 1-4'
        >>> 'mol. 1,2,3,4'
        >>> 'mol.atoms 200-300'
        >>> 'mol.a > 200'
        >>> 'mol.weight <=300.0'
        >>> 'mol.weight > 200.0'
        
        >>> 'chain.name A'
        >>> 'chain. A'
        >>> 'c. A'
        
        >>> 'res. ala val leu'
        >>> 'res. 1 2 3'
        >>> 'res.ptype gly,val,ala' # three-letter PDB code
        >>> 'res. gly val ala'
        >>> 'res.mtype g,v,a' # one-letter residue code in Maestro
        >>> 'res.m g,v,a'
        >>> 'residue.polarity hydrophobic'
        >>> 'residue.pol pos,neg'
        >>> 'res.pol h pos neg'
        >>> 'residue.sec helix # secondary structure
        >>> 'residue.sec hel, str'
        >>> 'res.sec l, s'
        >>> 'residue.pos 0.0 0.1' # faction position of the residue e.g. from the N-ter to 10% of the total residues
        
        >>> 'atom. 1,2,3,CA' # atom ptypes and numbers may be mixed
        >>> 'atom.ptype N,CA,C,O'
        >>> 'atom. N,CA,C,O'
        >>> 'a. n,ca,c,o'
        >>> 'atom.name the_36th_carbon'
        >>> atom.na C15, O:66, H-77
        >>> atom.na C* # returns atoms with name starting with C
        >>> atom.nam ??0* # returns atoms whose name's 3rd character is '0')
        >>> atom.num 1,2,3,4
        >>> a. 1 2 3 4
        >>> atom. 1-4
        >>> atom.molnum 1 # returns the first atom in each molecule
        >>> atom.molnum 1-10 # returns the first 10 atoms in each molecule
        >>> atom.entrynum 1 # returns the first atom in each entry
        >>> atom.entrynum 1-10 # returns the first 10 atoms in each entry
        >>> atom.mtype C2,O2 # Maestro atom type
        >>> atom.ele C,O # element
        >>> atom.att 1 # returns all terminal atoms. 
        >>> atom.att <=2 # returns all terminal atoms and other atoms with 2 or fewer bonds
        >>> atom.att 0 # returns all isolated atoms
        >>> atom.atomicnum 1 # returns all hydrogen atoms
        >>> atom.ato 1-6 # returns all atoms in the range H to C
        >>> atom.charge 0.400 # returns atoms with partial charges of 0.400. 
        >>> atom.charge -0.6--0.4 # returns atoms with partial charges -0.6 to -0.4,
        >>> atom.charge <0.0 # returns atoms with negative partial charges
        >>> atom.charge >=0.5 # returns atoms with charges of 0.5 or greater
        >>> atom.formal 0 # returns atoms with formal charges of 0
        >>> atom.formal -2 -- 1 # returns atoms with formal charges -2 to -1
        >>> atom.formal <0 # returns atoms with negative formal charges
        >>> atom.formal >=1 # returns atoms with formal charges of 1 or greater

        >>> 'smarts. [N]C(=O)[OH1]'
        >>> '(chain.name "A" AND res.num 833 AND atom.name "SG")'

    """
    parser = argparse.ArgumentParser(description="batch ASL expression check",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--asl', dest="asl", default="", help="ASL expression")
    parser.add_argument('cms', nargs="+", help="...-out.cms file")
    args = parser.parse_args()

    for cms_file in args.cms :
        
        st = structure.StructureReader.read(cms_file)
        
        if args.asl :
            for i, mol in enumerate(st, start=1):
                asl_aids = analyze.evaluate_asl(mol, args.asl)
                for j in asl_aids:
                    a = st.atom[j]
                    chain_.add(a.property['s_m_chain_name'])
                    resseq_.add(str(a.property['i_m_residue_number']))
                    resname_.add(a.property['s_m_pdb_residue_name'].strip())
                print("Molecule", i)
                print("    number of atoms= ", len(asl_aids))
                print("    chain(s): " + " ".join(chain_))
                print("    residue number(s): " + " ".join(resseq_))
                print("    residue name(s): " + " ".join(resname_))
        
        else: # detect ligand
            ligand_aids_list = analyze.find_ligands(st) # list of list
            for i, ligand_aids in enumerate(ligand_aids_list, start=1):
                chain_ = set()
                resseq_ = set()
                resname_ = set()
                for j in ligand_aids:
                    a = st.atom[j]
                    chain_.add(a.property['s_m_chain_name'])
                    resseq_.add(str(a.property['i_m_residue_number']))
                    resname_.add(a.property['s_m_pdb_residue_name'].strip())
                print("Ligand", i)
                print("    number of atoms= ", len(ligand_aids))
                print("    chain(s): " + " ".join(chain_))
                print("    residue number(s): " + " ".join(resseq_))
                print("    residue name(s): " + " ".join(resname_))



def batch_report():
    parser = argparse.ArgumentParser(description="batch desmond report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--protein', dest="protein", default="", help="protein ASL")
    parser.add_argument('-l', '--ligand', dest="ligand", default="none", help="ligand ASL")
    parser.add_argument('-s', '--slice-trj', dest="slice", default="", help="slice trajectory (start:end:step)")
    parser.add_argument('-w', '--overwrite', dest="overwrite", default=False, action="store_true", help="overwrite output(s)")
    parser.add_argument('-j', '--jobfile', dest="job_file", default="desmond_report_1.sh", help="job filename")
    parser.add_argument('cms', nargs="+", help="...-out.cms file")
    args = parser.parse_args()

    excluded_residues = ["POPC", "SPC", "T3P", "TIP3P", "T4P", "TIP4P", "NA", "K", "CL"]

    job_file = args.job_file
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"

    if args.slice:
        slice_trj = f'-slice-trj {args.slice}'
    else:
        slice_trj = ''

    # ex. md_1/r00/desmond_md_job_00_md2_m_naloxone-out.cms
    # mandatory input
    if not args.cms :
        parser.print_help()
        sys.exit(0)

    # provide options if protein and ligand ASL are not given
    if not (args.protein and args.ligand):
        ASL_choice = []
        protein_asl_index = None
        ligand_asl_index = None

        # read .cms and its associated trajectory file: <jobname>-out.cms, <jobname>_trj
        input_cms = args.cms[0]

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

        print()
        print("Title:", cms_model.title)
        print("Excluded residues:", ",".join(excluded_residues))
        print()

        # chain based ASL choices
        for chain in cms_model.chain:
            # r and s for the first and last residues
            for r in chain.residue:
                break
            res = [ s for s in chain.residue if s.molecule_number == r.molecule_number ]
            s = res[-1]
            if (r.pdbres.strip() in excluded_residues) or (s.pdbres.strip() in excluded_residues):
                continue
            if chain.name.strip() :
                ASL_choice.append((
                    f"mol. {r.molecule_number} and chain. {chain.name}",
                    f"{len(res):6d} residues ({r.pdbres}{r.resnum}...{s.pdbres}{s.resnum})",
                    ))
            else: # if chain.name is blank
                ASL_choice.append((
                    f"mol. {r.molecule_number}",
                    f"{len(res):6d} residues ({r.pdbres}{r.resnum}...{s.pdbres}{s.resnum})",
                    ))

        # molecule based ASL choices
        for molecule in cms_model.molecule:
            r = None # first residue
            s = None # last residue
            for s in molecule.residue:
                if not r:
                    r = s
                pass
            if r.pdbres.strip() in excluded_residues :
                continue
            ASL_choice.append((
                f"mol. {molecule.number}",
                f"{len(molecule.residue):6d} residues ({r.pdbres}{r.resnum}...{s.pdbres}{s.resnum})",
                ))

        # sort by ASL expression text
        ASL_choice = sorted(ASL_choice)
        ASL_choice.append(("none", ""))

        # show choices
        # default
        protein_asl_index = 0
        ligand_asl_index = len(ASL_choice)-1
        for idx, (asl, info) in enumerate(ASL_choice):
            print(f"[{idx}] {asl:20s} {info}")
        print()
        
        # set protein ASL
        args.protein = ASL_choice[protein_asl_index][0] # default
        ret = input(f'Enter protein of interest (number) [{protein_asl_index}] or ASL: ')
        ret = ret.strip()
        if ret:
            if re.match('^[0-9]', ret):
                protein_asl_index = int(ret)
                args.protein = ASL_choice[protein_asl_index][0]
            else:
                args.protein = ret
        # check protein ASL
        protein_atoms = analyze.evaluate_asl(cms_model, args.protein)
        print(f"  -prot {args.protein}  [{len(protein_atoms)} atoms]")

        # set ligand ASL
        args.ligand = ASL_choice[ligand_asl_index][0] # default
        ret = input(f'Enter ligand  of interest (number) [{ligand_asl_index}] or ASL: ')
        ret = ret.strip()
        if ret:
            if re.match('^[0-9]', ret):
                ligand_asl_index = int(ret)
                args.ligand = ASL_choice[ligand_asl_index][0]
            else:
                args.ligand = ret   
        # check ligand ASL
        if not args.ligand.startswith("none"):  
            ligand_atoms = analyze.evaluate_asl(cms_model, args.ligand)
            print(f"  -lig {args.ligand}  [{len(ligand_atoms)} atoms]")

    total_cms_files = len(args.cms)
    count_cms_files = 0
    with open(job_file,"w") as job:
        for cms_file in args.cms :
            cms_base = os.path.basename(cms_file).replace("-out.cms","")
            cms_prefix = cms_file[:-8]
            trj_dir  = f"{cms_prefix}_trj"
            name = cms_base.replace("desmond_md_job_", "")
            count_cms_files += 1

            if (not os.path.exists(f'{name}-in.eaf')) or args.overwrite:
                job.write(f'{run} event_analysis.py analyze \\\n')
                job.write(f'\t{cms_file} \\\n')
                job.write(f'\t-prot "{args.protein}" \\\n')
                job.write(f'\t-lig "{args.ligand}" \\\n')
                job.write(f'\t-out {name}\n\n')

            if (not os.path.exists(f'{name}-out.eaf')) or args.overwrite:
                """
                -s START:END:STEP, -slice-trj START:END:STEP
                                Use the sliced trajectory. We use Python's slice
                                notation. START, END, and STEP should be integer
                                numbers.
                """
                job.write(f'{run} analyze_simulation.py \\\n')
                job.write(f'\t{slice_trj} {cms_file} \\\n')
                job.write(f'\t{trj_dir} \\\n')
                job.write(f'\t{name}-out.eaf \\\n')
                job.write(f'\t{name}-in.eaf\n\n')

            if (not os.path.exists(f'report_{name}.pdf')) or args.overwrite:
                job.write(f'{run} event_analysis.py report -pdf \\\n')
                job.write(f'\treport_{name}.pdf \\\n')
                job.write(f'\t{name}-out.eaf\n\n')

    os.chmod(job_file, 0o777)


def batch_dihedral() -> None:
    """Analyze dihedral and pseudo-dihedrals
    """
    parser = argparse.ArgumentParser(description="batch trajectory dihedral angles",
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
                dihed_aids = get_dihedral_atom_indices(cms_model, args.chain, resSeq, angle_name)
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


def batch_ligrmsd() -> None:
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


def batch_rg() -> None:
    """Calculate Radius of Gyration (Rg).
    """
    parser = argparse.ArgumentParser(description="batch trajectory radius of gyration (Rg)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--asl', dest="asl", action="append", help="atom selection language(s)")
    parser.add_argument('--start', dest="start", type=float, default=0.0, help="start fraction for analysis")
    parser.add_argument('--end', dest="end", type=float, default=1.0, help="end fraction for analysis")
    parser.add_argument('--n', dest="n", type=int, default=0, help="target sample points")
    parser.add_argument('--csv', dest="csv", default="rg.csv.gz", help="output csv or csv.gz filename")
    parser.add_argument('--check', dest="check", default=False, action="store_true", help="just check atom index")
    parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
    args = parser.parse_args()

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
