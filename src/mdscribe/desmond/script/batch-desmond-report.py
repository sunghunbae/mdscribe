import os
import sys
import argparse
import subprocess
import re


from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.structutils import analyze

run = os.environ["SCHRODINGER"]+"/run"

excluded_residues = ["POPC", "SPC", "T3P", "TIP3P", "T4P", "TIP4P", "NA", "K", "CL"]

parser = argparse.ArgumentParser(description="batch gdesmond md jobs",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p', '--protein', dest="protein", default="", help="protein ASL")
parser.add_argument('-l', '--ligand', dest="ligand", default="none", help="ligand ASL")
parser.add_argument('-s', '--slice-trj', dest="slice", default="", help="slice trajectory (start:end:step)")
parser.add_argument('-w', '--overwrite', dest="overwrite", default=False, action="store_true", help="overwrite output(s)")
parser.add_argument('-j', '--jobfile', dest="job_file", default="desmond_report_1.sh", help="job filename")
parser.add_argument('cms', nargs="+", help="...-out.cms file")
args = parser.parse_args()

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