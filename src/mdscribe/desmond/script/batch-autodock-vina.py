#!/usr/bin/env python3

import pathlib
import pandas as pd
import argparse
import subprocess
import sys
import numpy as np
import io
import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import namedtuple

MolLines = namedtuple("MolLines", ["start", "end"])


def setup_docking_box(args):
    if args.glidegrid:
        """
        Binding pocket information from GLIDE receptor grid
        (ex) ../../../share/glide-grid-7EPZ-B/glide-grid-7EPZ-B.in

        FORCEFIELD   OPLS_2005
        GRID_CENTER   148.44219276056342, 147.82381628169014, 119.45476397183101
        GRIDFILE   glide-grid-7EPZ-B.zip
        INNERBOX   20, 15, 15
        OUTERBOX   38.23687762887996, 33.23687762887996, 33.23687762887996
        RECEP_FILE   glide-grid-7EPZ-B.maegz
        """
        with open(args.glidegrid, "r") as f:
            for line in f:
                if line.startswith("GRID_CENTER"):
                    xyz = line.strip().replace(",", "").split()[1:]
                    [args.center_x, args.center_y, args.center_z] = map(float, xyz)
                if line.startswith("INNERBOX"):
                    xyz = line.strip().replace(",", "").split()[1:]
                    [ix, iy, iz] = map(float, xyz)
                if line.startswith("OUTERBOX"):
                    xyz = line.strip().replace(",", "").split()[1:]
                    [ox, oy, oz] = map(float, xyz)
        [args.size_x, args.size_y, args.size_z] = 0.5*np.array([ix+ox, iy+oy, iz+oz])
    
    print(f"Center = {args.center_x:8.3f} {args.center_y:8.3f} {args.center_z:8.3f}")
    print(f"Size   = {args.size_x:8.3f} {args.size_y:8.3f} {args.size_z:8.3f}")


def openbabel_version(args):
    """returns openbabel version"""
    p = subprocess.run((args.obabel, '-V'), capture_output=True, text=True)
    if p.stderr:
        print(f"openbabel error: {p.stderr.rstrip()}")
        return "0.0.0"
    # ex. p.stdout = Open Babel 3.1.0 -- Sep 24 2023 -- 08:44:32
    return p.stdout.split()[2]


def convert_ligand_mol2_to_pdbqt(args):
    p = subprocess.run([args.obabel, args.lig, '-opdbqt'], capture_output=True, text=True)
    lines = p.stdout.splitlines()

    # split converted pdbqt
    molecules = []
    for i, line in enumerate(lines):
        if line.startswith("MODEL"):
            mol = MolLines(i+1, 0) # exclude this MODEL line
        elif line.startswith("ENDMDL"):
            mol.end = i # this ENDMDL line will be excluded by python slice rule
            molecules.append(mol)

    pathlib.Path(args.workdir / f"{args.jobname}-ligands").mkdir(parents=True, exist_ok=True)

    for i in range(len(molecules)):
        pdbqt_outfile = (args.workdir / f"{args.jobname}-ligands" / f"mol-{i+1:04d}.pdbqt").as_posix()
        with open(pdbqt_outfile, "w") as f:
            # both MODEL and ENDMDL lines are exlcuded in the final .pdbqt
            f.write('\n'.join(lines[molecules[i].start:molecules[i].end]))


def run_vina_batch(args):
    """
    #################################################################
    # If you used AutoDock Vina in your work, please cite:          #
    #                                                               #
    # J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli  #
    # AutoDock Vina 1.2.0: New Docking Methods, Expanded Force      #
    # Field, and Python Bindings, J. Chem. Inf. Model. (2021)       #
    # DOI 10.1021/acs.jcim.1c00203                                  #
    #                                                               #
    # O. Trott, A. J. Olson,                                        #
    # AutoDock Vina: improving the speed and accuracy of docking    #
    # with a new scoring function, efficient optimization and       #
    # multithreading, J. Comp. Chem. (2010)                         #
    # DOI 10.1002/jcc.21334                                         #
    #                                                               #
    # Please see https://github.com/ccsb-scripps/AutoDock-Vina for  #
    # more information.                                             #
    #################################################################

    Input:
    --receptor arg             rigid part of the receptor (PDBQT)
    --flex arg                 flexible side chains, if any (PDBQT)
    --ligand arg               ligand (PDBQT)
    --batch arg                batch ligand (PDBQT)
    --scoring arg (=vina)      scoring function (ad4, vina or vinardo)

    Search space (required):
    --maps arg                 affinity maps for the autodock4.2 (ad4) or vina
                                scoring function
    --center_x arg             X coordinate of the center (Angstrom)
    --center_y arg             Y coordinate of the center (Angstrom)
    --center_z arg             Z coordinate of the center (Angstrom)
    --size_x arg               size in the X dimension (Angstrom)
    --size_y arg               size in the Y dimension (Angstrom)
    --size_z arg               size in the Z dimension (Angstrom)
    --autobox                  set maps dimensions based on input ligand(s) (for
                                --score_only and --local_only)

    Output (optional):
    --out arg                  output models (PDBQT), the default is chosen based
                                on the ligand file name
    --dir arg                  output directory for batch mode
    --write_maps arg           output filename (directory + prefix name) for
                                maps. Option --force_even_voxels may be needed to
                                comply with .map format

    Misc (optional):
    --cpu arg (=0)             the number of CPUs to use (the default is to try
                                to detect the number of CPUs or, failing that, use
                                1)
    --seed arg (=0)            explicit random seed
    --exhaustiveness arg (=8)  exhaustiveness of the global search (roughly
                                proportional to time): 1+
    --max_evals arg (=0)       number of evaluations in each MC run (if zero,
                                which is the default, the number of MC steps is
                                based on heuristics)
    --num_modes arg (=9)       maximum number of binding modes to generate
    --min_rmsd arg (=1)        minimum RMSD between output poses
    --energy_range arg (=3)    maximum energy difference between the best binding
                                mode and the worst one displayed (kcal/mol)
    --spacing arg (=0.375)     grid spacing (Angstrom)
    --verbosity arg (=1)       verbosity (0=no output, 1=normal, 2=verbose)

    Configuration file (optional):
    --config arg               the above options can be put here

    Information (optional):
    --help                     display usage summary
    --help_advanced            display usage summary with advanced options
    --version                  display program version

    R. Quiroga, M. A. Villarreal, Vinardo: A Scoring Function Based on 
    Autodock Vina Improves Scoring, Docking, and Virtual Screening. 
    PLoS One 11, e0155183 (2016).

    
    example
    ${vina} --receptor 7EPZ-B_prepared.pdbqt \
      --batch ligands_2/ligprep-???.pdbqt \
      --scoring vina \
      --center_x 148.44219276056342 \
      --center_y 147.82381628169014 \
      --center_z 119.45476397183101 \
      --size_x 30.0 \
      --size_y 25.0 \
      --size_z 25.0 \
      --exhaustiveness 32 \
      --dir output
    """

    try:
        assert args.rec
    except:
        print("receptor .pdbqt is required")
        sys.exit(0)
    try:
        assert args.center_x
        assert args.center_y
        assert args.center_z
    except:
        print("docking box should be defined")
        sys.exit(0)

    subprocess.run([
        args.vina,
        '--verbosity', '0',
        '--receptor', args.rec,
        '--scoring', args.scoring,
        '--center_x', str(args.center_x),
        '--center_y', str(args.center_y),
        '--center_z', str(args.center_z),
        '--size_x', str(args.size_x),
        '--size_y', str(args.size_y),
        '--size_z', str(args.size_z),
        '--exhaustiveness', str(args.exhaustiveness),
        '--dir', args.dir,
        '--batch'] + 
        [v.as_posix() for v in sorted(args.workdir.glob(f"{args.jobname}-ligands/mol-*.pdbqt"))]
        )


def collect_vina_scores(args):
    """returns a pandas DataFrame of all vina scores
    Example:
    REMARK VINA RESULT:    -8.761      0.000      0.000
    REMARK VINA RESULT:    -8.750      5.164     11.063
    REMARK VINA RESULT:    -8.413      6.516      8.780
    REMARK VINA RESULT:    -8.255      1.546      3.964
    REMARK VINA RESULT:    -8.016      6.838     10.588
    REMARK VINA RESULT:    -8.007      7.359     11.015
    REMARK VINA RESULT:    -8.001      2.366      4.707
    REMARK VINA RESULT:    -7.938      4.287      7.180
    REMARK VINA RESULT:    -7.923      7.978     10.234

    1st column - binding free energy
    2nd column - RMSD Lower bound
    3rd column - RMSD upper bound
    Both RMSD values are calculated relative to the best model
    """
    data = {'name':[], 'mol':[], 'model':[], 'vina':[], 'path':[] }
    for filename in tqdm.tqdm(sorted(pathlib.Path(args.dir).glob("*_out.pdbqt"))):
        mol = int(pathlib.Path(filename).name.replace("mol-","").replace("_out.pdbqt","")) # 1-based
        with open(filename, "r") as f:
            model = None
            name = None
            dG = None
            for line in f:
                if line.startswith("MODEL"):
                    cols = line.strip().split()
                    model = int(cols[1])
                    continue
                if line.startswith("REMARK VINA RESULT:"):
                    cols = line.strip().split()
                    [dG, rmsd_lb, rmsd_ub] = cols[3:6]
                    continue
                if line.startswith("REMARK  Name ="):
                    cols = line.strip().split()
                    name = cols[3]
                    data['name'].append(name)
                    data['mol'].append(mol)
                    data['vina'].append(float(dG))
                    data['model'].append(model)
                    data['path'].append(pathlib.Path(filename).resolve())
                    model = None
                    name = None
                    dG = None
                    continue

    return pd.DataFrame(data)


def replace_mol2_coord_dict(mol2block:list, coord:dict):
    with io.StringIO() as f:
        atom_block = False
        for line in mol2block:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_block = True
                f.write(line)
                continue
            elif line.startswith("@<TRIPOS>BOND"):
                atom_block = False
            if not atom_block:
                f.write(line)
                continue
            # replace coordinates in the atom block
            cols= line.strip().split()
            if len(cols) == 9:
                (atom_id,atom_name,x,y,z,atom_type,subst_id,subst_name,charge) = cols
                x = float(x)
                y = float(y)
                z = float(z)
                charge = float(charge)
                if atom_name in coord:
                    (x, y, z) = coord[atom_name]
                f.write(f"{atom_id:>7} {atom_name:<8} {x:9.4f} {y:9.4f} {z:9.4f} " \
                        f"{atom_type:<8} {subst_id:>2} {subst_name:<9} {charge:8.4f}\n")
            else:
                f.write(line)
        new_mol2block = f.getvalue()
    
    return new_mol2block


def replace_mol2_coord_list(mol2block:list, coord:list):
    with io.StringIO() as f:
        atom_block = False
        for line in mol2block:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_block = True
                f.write(line)
                continue
            elif line.startswith("@<TRIPOS>BOND"):
                atom_block = False
            if not atom_block:
                f.write(line)
                continue
            # replace coordinates in the atom block
            cols= line.strip().split()
            if len(cols) == 9:
                (atom_id,atom_name,x,y,z,atom_type,subst_id,subst_name,charge) = cols
                idx = int(atom_id) - 1
                charge = float(charge)
                (x, y, z) = coord[idx]
                f.write(f"{atom_id:>7} {atom_name:<8} {x:9.4f} {y:9.4f} {z:9.4f} " \
                        f"{atom_type:<8} {subst_id:>2} {subst_name:<9} {charge:8.4f}\n")
            else:
                f.write(line)
        new_mol2block = f.getvalue()
    
    return new_mol2block


def export_docked_pose_to_mol2(args, df:pd.DataFrame, outdir):
    """overwriting coordinates of original mol2 with those of pdbqt

    In this workflow, .pdbqt files were originally converted from .mol2.

    We'd like to convert vina output .pdbqt to .mol2 files because 
    .mol2 files are most convenient for preparing MD simulations following 
    vina docking in MM/GBSA or binding free energy calculations.
    
    PDBQT files do not contain hydrogens and causes trouble when 
    converting PDBQT back to MOL2 by using obabel. So, we can use 
    the original .mol2 and import the XYZ coordinates from .pdbqt
    and make the final .mol2 keeping all hydrogens and bond connectivities.
    """

    # read the original .mol2 file
    with open(args.lig, "r") as f:
        orig_mol2_lines = f.readlines()
        orig_mol2 = [i for i, line in enumerate(orig_mol2_lines) if line.strip().startswith("@<TRIPOS>MOLECULE")]

    (pathlib.Path(args.workdir) / outdir).mkdir(parents=True, exist_ok=True)

    for pdbqt_outfile, model, mol in tqdm.tqdm(zip(df_repr.path, df_repr.model, df_repr.mol)):
        mol2_outfile = (pathlib.Path(args.workdir) / outdir / f"mol-{mol:04d}_out.mol2").as_posix()
        
        # get coordinates from vina .pdbqt output
        coord = {}
        with open(pdbqt_outfile, "r") as f:
            export = False
            for line in f:
                if line.startswith('MODEL'):
                    if model == int(line.split()[1]):
                        export = True
                    else:
                        export = False
                if export:
                    if line.startswith("ATOM"):
                        cols = line.split()
                        atom_name = cols[2]
                        [x, y, z] = map(float, cols[6:9])
                        coord[atom_name] = (x, y, z)

        start = orig_mol2[mol-1]
        end = len(orig_mol2_lines) if (mol == len(orig_mol2)) else orig_mol2[mol]
    
        new_mol2block = replace_mol2_coord_dict(orig_mol2_lines[start:end], coord)
        
        # fix hydrogen atom positions using energy minimization
        rdmol = Chem.MolFromMol2Block(new_mol2block, removeHs=False)

        #  rdkit: Can't kekulize mol errors
        if rdmol is None:
            continue

        rdmp = AllChem.MMFFGetMoleculeProperties(rdmol)
        rdff = AllChem.MMFFGetMoleculeForceField(rdmol, rdmp)
        for a in rdmol.GetAtoms():
            if a.GetAtomicNum() > 1: 
                # fix heavy atoms
                # force constant should be too high
                rdff.MMFFAddPositionConstraint(a.GetIdx(), 0.0, 1.0e+4)
        rdff.Minimize(maxIts=2000)
        
        for conf  in rdmol.GetConformers():
            coord_list = conf.GetPositions()
        
        new_mol2block = replace_mol2_coord_list(orig_mol2_lines[start:end], coord_list)
        
        with open(mol2_outfile, "w") as f:
            f.write(new_mol2block)


def export_vina_pose(pdbqt_infile:str, model:int, output_dir:str, 
                     name:str, mol2:bool=True):
    """extract a model from a vina output .pdbqt"""
    try:
        pathlib.Path(output_dir).mkdir()
    except FileExistsError:
        pass

    with open(pdbqt_infile, "r") as f, io.StringIO() as pdbqt_model:
        export = False
        for line in f:
            if line.startswith('MODEL'):
                if model == int(line.split()[1]):
                    export = True
                else:
                    export = False
            if export:
                pdbqt_model.write(line)
        pdbqt_output = pdbqt_model.getvalue()

    if pdbqt_output:
        if mol2:
            outfile = (pathlib.Path(output_dir) / f"{name}.mol2").as_posix()
            # convert .pdbqt to .mol2 using openbabel
            p = subprocess.Popen(['obabel', '-i','pdbqt', '-o','mol2'], 
                            stdin=subprocess.PIPE, 
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)
            mol2_output, obabel_err = p.communicate(pdbqt_output)
            with open(outfile, "w") as f:
                f.write(mol2_output)
        else:
            outfile = (pathlib.Path(output_dir) / f"{name}.pdbqt").as_posix()
            with open(outfile, "w") as f:
                f.write(pdbqt_output)



if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description='Run Autodock Vina in Batch Mode', 
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("--vina", dest="vina", help="vina executable path", 
                           default="/home2/shbae/local/bin/vina")
    argparser.add_argument("--obabel", dest="obabel", help="obabel executable path", 
                           default="/home2/shbae/miniconda3/envs/openmm/bin/obabel")
    argparser.add_argument("--scoring", dest="scoring", 
                           help="scoring function (ad4, vina or vinardo)", default="vina")
    argparser.add_argument("--exhaustiveness", dest="exhaustiveness", 
                           help="exhaustiveness of the global search", default=32)
    argparser.add_argument("--glidegrid", dest="glidegrid", 
                           help="Glide grid .in file to setup docking box")
    argparser.add_argument("--center_x", dest="center_x", help="X coordinate of the center (Angstrom)")
    argparser.add_argument("--center_y", dest="center_y", help="Y coordinate of the center (Angstrom)")
    argparser.add_argument("--center_z", dest="center_z", help="Z coordinate of the center (Angstrom)")
    argparser.add_argument("--size_x", dest="size_x", help="size in the X dimension (Angstrom)", default=25.0)
    argparser.add_argument("--size_y", dest="size_y", help="size in the Y dimension (Angstrom)", default=25.0)
    argparser.add_argument("--size_z", dest="size_z", help="size in the Z dimension (Angstrom)", default=25.0)
    argparser.add_argument("--dir", dest="dir", help="vina output directory", default="vina-output")
    argparser.add_argument("--workdir", dest="workdir", help="working directory", default=pathlib.Path("."))
    argparser.add_argument("--jobname", dest="jobname", help="jobname", default="vina")
    argparser.add_argument("--rec", dest="rec", help="receptor PDBQT")
    argparser.add_argument("--lig", dest="lig", help="ligands MOL2")
    args = argparser.parse_args()

    if not (openbabel_version(args) > "3.1"):
        print("openbabel version 3.1 or later is required")
        sys.exit(0)


    setup_docking_box(args)

    convert_ligand_mol2_to_pdbqt(args)

    run_vina_batch(args)

    df = collect_vina_scores(args)
    df.to_csv(args.workdir / f"{args.jobname}-out.csv", index=False, float_format="%7.3f")
    
    df_repr = df.iloc[df.groupby("name").vina.idxmin()]
    df_repr.to_csv(args.workdir / f"{args.jobname}-repr-out.csv", index=False, float_format="%7.3f")
    
    export_docked_pose_to_mol2(args, df_repr, outdir=f"{args.jobname}-pose")
