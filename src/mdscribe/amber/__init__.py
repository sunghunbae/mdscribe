__all__ = [ 'check_ambertools', 'AmberSetup', 'run_amber_setup', 'run_mmpbsa_py' ]

import argparse
import pathlib
import subprocess
import logging
from collections import namedtuple
import shutil

MMPBSA_py_result = namedtuple('MMPBSA_py_result', 
                              ['ligand', 
                               'dG', 
                               'stdev', 
                               'sem', 
                               'output_path',
                               ])



def check_ambertools(logger:logging.Logger) -> None:
    """check ambertools exist"""
    if shutil.which('antechamber') is None:
        logger.error(f"{'antechamber':<20} not available")
    else:
        logger.info(f"{'antechamber':<20} {shutil.which('antechamber')}")
    
    if shutil.which('parmchk2') is None:
        logger.error(f"{'parmchk2':<20} not available")
    else:
        logger.info(f"{'parmchk2':<20} {shutil.which('parmchk2')}")
    
    if shutil.which('tleap') is None:
        logger.error(f"{'tleap':<20} not available")
    else:
        logger.info(f"{'tleap':<20} {shutil.which('tleap')}")

    if shutil.which('MMPBSA.py') is None:
        logger.error(f"{'MMPBSA.py':<20} not available")
    else:
        logger.info(f"{'MMPBSA.py':<20} {shutil.which('MMPBSA.py')}")
    


class AmberSetup(object):
    def __init__(self, 
                 args: argparse.Namespace, 
                 ligand_mol2: pathlib.Path | None = None, 
                 tempdir: pathlib.Path | None = None) -> None:
        self.ff_protein = args.ff_protein
        self.ff_water = args.ff_water
        self.ff_ligand = args.ff_ligand
        self.charge = args.charge
        self.salt = args.salt
        self.water = args.water
        self.gbsa = args.gbsa
        self.gbsa_line = "set default PBRadii mbondi2" if self.gbsa else ""
        self.receptor_pdb = args.receptor_pdb
        self.receptor_name = args.receptor_name
        if ligand_mol2 is None:
            self.ligand_mol2_orig = ''
            self.ligand_mol2 = ''
            self.ligand_frcmod = ''
            self.ligand_name = ''
            self.complex_name = ''
        else:
            # input .mol2 filename
            self.ligand_mol2_orig = ligand_mol2
            name =self.ligand_mol2_orig.name.replace("_out.mol2", "").replace(".mol2", "")
            # output filenames
            self.ligand_mol2 = f"{name}_ac.mol2"
            self.ligand_frcmod = f"{name}_ac.frcmod"
            self.ligand_name = f"{name}"
            self.complex_name = f"{name}-complex"
        
        # temporary working directory
        self.tempdir = args.workdir if tempdir is None else tempdir

        self.tleap_input = {
            "complex-waterbox":
                (f"""source leaprc.protein.{self.ff_protein}
                    source leaprc.water.{self.ff_water}
                    source leaprc.{self.ff_ligand}
                    loadAmberParams {self.ligand_frcmod}
                    ligand = loadMol2 {self.ligand_mol2}
                    receptor = loadPDB {self.receptor_pdb}
                    complex = combine {{receptor ligand}}
                    {self.gbsa_line}
                    charge complex
                    saveAmberParm ligand {self.ligand_name}.prmtop {self.ligand_name}.inpcrd
                    saveAmberParm receptor {self.receptor_name}.prmtop {self.receptor_name}.inpcrd
                    saveAmberParm complex {self.complex_name}.prmtop {self.complex_name}.inpcrd
                    solvateOct complex TIP3PBOX {self.water}
                    saveAmberParm complex {self.complex_name}_waterbox.prmtop {self.complex_name}_waterbox.inpcrd
                    quit""", "leap_com_solvated.in"),
            "complex-solvated":
                (f"""source leaprc.protein.{self.ff_protein}
                    source leaprc.water.{self.ff_water}
                    source leaprc.{self.ff_ligand}
                    loadAmberParams {self.ligand_frcmod}
                    ligand = loadMol2 {self.ligand_mol2}
                    receptor = loadPDB {self.receptor_pdb}
                    complex = combine {{receptor ligand}}
                    {self.gbsa_line}
                    solvateBox complex TIP3PBOX {self.water}
                    addIons complex Cl- __Cl__
                    addIons complex Na+ __Na__
                    charge complex
                    saveAmberParm complex {self.complex_name}_solvated.prmtop {self.complex_name}_solvated.inpcrd
                    quit""", "leap_com_solvated_ions.in"),
            "receptor-waterbox":
                (f"""source leaprc.protein.{self.ff_protein}
                    source leaprc.water.{self.ff_water}
                    receptor = loadPDB {self.receptor_pdb}
                    {self.gbsa_line}
                    charge receptor
                    solvateOct receptor TIP3PBOX {self.water}
                    saveAmberParm receptor {self.receptor_name}_waterbox.prmtop {self.receptor_name}_waterbox.inpcrd
                    quit""", "leap_rec_solvated.in"),
            "receptor-solvated":
                (f"""source leaprc.protein.{self.ff_protein}
                    source leaprc.water.{self.ff_water}
                    receptor = loadPDB {self.receptor_pdb}
                    {self.gbsa_line}
                    solvateOct receptor TIP3PBOX {self.water}
                    addIons receptor Cl- __Cl__
                    addIons receptor Na+ __Na__
                    charge receptor
                    saveAmberParm receptor {self.receptor_name}_solvated.prmtop {self.receptor_name}_solvated.inpcrd
                    quit""", "leap_rec_solvated_ions.in"),
            "ligand-waterbox":
                (f"""source leaprc.water.{self.ff_water}
                    source leaprc.{self.ff_ligand}
                    loadAmberParams {self.ligand_frcmod}
                    ligand = loadMol2 {self.ligand_mol2}
                    {self.gbsa_line}
                    charge ligand
                    solvateOct ligand TIP3PBOX {self.water}
                    saveAmberParm ligand {self.ligand_name}_waterbox.prmtop {self.ligand_name}_waterbox.inpcrd
                    quit""", "leap_lig_solvated.in"),
            "ligand-solvated":
                (f"""source leaprc.water.{self.ff_water}
                    source leaprc.{self.ff_ligand}
                    loadAmberParams {self.ligand_frcmod}
                    ligand = loadMol2 {self.ligand_mol2}
                    {self.gbsa_line}
                    solvateOct ligand TIP3PBOX {self.water}
                    addIons ligand Cl- __Cl__
                    addIons ligand Na+ __Na__
                    charge ligand
                    saveAmberParm ligand {self.ligand_name}_solvated.prmtop {self.ligand_name}_solvated.inpcrd
                    quit""", "leap_lig_solvated_ions.in")
        }

    def get_net_charge(self) -> tuple[int, list]:
        with open(self.ligand_mol2_orig, "r") as f:
            mol2_charges = []
            atom_block = False
            for line in f:
                if line.startswith("@<TRIPOS>ATOM"):
                    atom_block = True
                    continue
                elif line.startswith("@<TRIPOS>BOND"):
                    atom_block = False
                if atom_block:
                    mol2_charges.append(float(line.strip().split()[8])) # 9th column
            q = sum(mol2_charges)
            net_charge = int(q+0.5) if q > 0 else int(q-0.5)
            return (net_charge, mol2_charges)


    def antechamber(self) -> None:
        """run antechamber and parmchk2
        
        Externally generated charges can be used by a charge file
        $ antechamber -fi mol2 -i in.mol2 -fo mol2 -o out.mol2 -c rc -cf input.chg -at gaff2 -nc 0
        
        input.chg can be generated by antechamber:
        $ antechamber -fi mol2 -i in.mol2 -c wc -cf input.chg

        [input.chg]
        -0.177500  0.468500 -0.177500 -0.177500 -0.143200 -0.046700 -0.186600  0.186500
        -0.186600 -0.046700 -0.336000  0.118600  0.479900 -0.505100 -0.367000  0.029400
         0.074700 -0.702800  0.074700  0.029400  0.224100 -0.171100 -0.116200  0.011300
        -0.068100 -0.083400 -0.196400  0.186000 -0.348300  0.222100 -0.227100 -0.180200
        -0.180200 -0.104500 -0.139200 -0.115300 -0.131500  0.147200  0.146600  0.146600
         0.147200  0.054600  0.054600  0.090700  0.090700  0.070200  0.070200  0.070200
         0.070200  0.090700  0.090700  0.046400  0.068700  0.068700  0.068700  0.134300
         0.135100  0.054600  0.054600  0.114400  0.099900  0.099900  0.099900  0.099900
         0.131200  0.131000  0.131300  0.130600

        note that charges are listed as the same order as in the in.mol2 file
        """

        (net_charge, mol2_charges) = self.get_net_charge()

        if self.charge == 'mol2':
            # write charges taken from input .mol2 file to a file
            # and use for antechamber
            input_charge_path = self.tempdir / "input.charges"
            with input_charge_path.open("w") as f:
                for c in mol2_charges:
                    f.write(f"{c:9.6f} ")
            p = subprocess.run(['antechamber', 
                            '-i', self.ligand_mol2_orig.as_posix(), 
                            '-fi', 'mol2', 
                            '-o', self.ligand_mol2, 
                            '-fo', 'mol2',
                            '-c', 'rc', # read charge
                            '-cf',  input_charge_path.as_posix(), # charge file
                            '-s', '2', # status information 0:brief, 1:default, 2:verbose
                            '-at', self.ff_ligand, 
                            '-nc', str(net_charge)],
                            cwd=self.tempdir.as_posix(),
                            capture_output=True) # redirect stdout and stderr
        
        elif self.charge in ['gas', 'bcc', 'mul']:
            p = subprocess.run(['antechamber', 
                            '-i', self.ligand_mol2_orig.as_posix(), 
                            '-fi', 'mol2', 
                            '-o', self.ligand_mol2, 
                            '-fo', 'mol2',
                            '-c', self.charge, 
                            '-s', '2', # status information 0:brief, 1:default, 2:verbose
                            '-at', self.ff_ligand, 
                            '-nc', str(net_charge)],
                            cwd=self.tempdir.as_posix(),
                            capture_output=True) # redirect stdout and stderr
        
        p = subprocess.run(['parmchk2', 
                        '-i', self.ligand_mol2, 
                        '-f', 'mol2', 
                        '-o', self.ligand_frcmod], 
                        cwd=self.tempdir.as_posix(),
                        capture_output=True) # redirect stdout and stderr


    def antechamber_outfiles_exist(self) -> bool:
        """return True if all required files exist"""
        name = self.ligand_name
        return all([ p.exists() for p in [
                self.tempdir / f"{name}_ac.mol2", 
                self.tempdir / f"{name}_ac.frcmod",
                ]
            ])
    

    def num_pos_neg_ions(self, leap_output:str) -> tuple[int, int]:
        """return number of postive and negative ions to make the salt concentration
        How many ions to add for 0.15 M ?
        
        Matias Machado's split method
        
        https://doi.org/10.1021/acs.jctc.9b00953
        http://archive.ambermd.org/202002/0194.html
        
        1) Calculate the expected number of ions (No) for a macroscopic
           salt concentration (Co, in molar units) in a simulation box
           of Nw water molecules:
        
           No = Nw*Co/56 # Best proved No estimator
        
        2) Get the actual number of positive (N+) and negative (N-) monovalent
           ions to add in the box due to the the solute's charge (Q):
        
           N+ = No - Q/2 # round up in case of odd Q (e.g. N+ = 2.5 ~ 3)
           N- = No + Q/2 # round up in case of odd Q (e.g. N- = 5.5 ~ 6)
        
        example:
           Lysozyme (PDB: 2VB1), Q = +8e, Water box 12Å => Nw ~ 7000 molecules, [NaCl]=0.15M:
           No = 7000*0.15/56 ~19 => No/Q > 1
           Then,
           Na+ = 19 -(+8)/2 = 15
           Cl- = 19 +(+8)/2 = 23
        """

        for line in leap_output.splitlines():
            line = line.strip()
            if line.startswith("Total unperturbed charge:"):
                Q = float(line.split()[3])
            elif line.startswith("WAT"):
                Nw = int(line.split()[1])
        
        No = Nw * self.salt / 56.
        
        pos = int(No - Q/2.0 +0.5)
        neg = int(No + Q/2.0 +0.5)
        
        system_charge = Q + pos - neg 
        # logger.info(f"system net charge= {system_charge:6.3f} positive= {pos} negative= {neg}")
        # system_charge should be close to zero
        # if not, adjust number of positive or negative ions
        while system_charge < -0.5 or system_charge > 0.5:
            if system_charge < -0.5:
                pos += 1
            elif system_charge > 0.5:
                neg += 1
            system_charge = Q + pos - neg
            # logger.info(f"system net charge= {system_charge:6.3f} positive= {pos} negative={neg}")

        return (pos, neg)
    

    def leap(self, kind:str) -> None:
        """run tLeap"""
        assert kind in ['receptor','ligand', 'complex']
        if kind == "receptor":
            name = self.receptor_name
        elif kind == "complex":
            name = self.complex_name
        else:
            name = self.ligand_name
            
        # solvate system (+ tip3box)
        (script, filename) = self.tleap_input[f"{kind}-waterbox"]
        with open((self.tempdir / filename).as_posix(), "w") as f:
            f.write('\n'.join([l.strip() for l in script.splitlines()]))
        p = subprocess.run(['tleap', '-f', filename], capture_output=True, text=True, cwd=self.tempdir.as_posix())
        (Na, Cl) = self.num_pos_neg_ions(leap_output=p.stdout)
        
        # remove intermediate output files
        (self.tempdir / f"{name}_waterbox.prmtop").unlink()
        (self.tempdir / f"{name}_waterbox.inpcrd").unlink()

        # solvated system (+ ions)
        (script, filename) = self.tleap_input[f"{kind}-solvated"]
        script = script.replace("__Cl__", str(Cl)).replace("__Na__", str(Na))
        with open((self.tempdir / filename).as_posix(), "w") as f:
            f.write('\n'.join([l.strip() for l in script.splitlines()]))
        p = subprocess.run(['tleap', '-f', filename], capture_output=True, text=True, cwd=self.tempdir.as_posix())


    def leap_outfiles_exist(self, kind) -> bool:
        """return True if output files exist"""
        assert kind in ['receptor','ligand', 'complex']
        
        if kind == "receptor":
            name = self.receptor_name
        elif kind == "complex":
            name = self.complex_name
        else:
            name = self.ligand_name
        
        return all([ p.exists() for p in [
            self.tempdir / f"{name}.prmtop",
            self.tempdir / f"{name}.inpcrd",
            self.tempdir / f"{name}_solvated.prmtop",
            self.tempdir / f"{name}_solvated.inpcrd",
            ]])



def run_amber_setup(args: argparse.Namespace, 
                    idx: int, 
                    ligand_mol2: pathlib.Path | None) -> list[ pathlib.Path ]:
    
    tempdir = args.workdir / f"temp-{idx:04d}"
    tempdir.mkdir()

    abt = AmberSetup(args, ligand_mol2=ligand_mol2, tempdir=tempdir)

    if ligand_mol2 is None: # receptor
        if args.leap:
            if args.overwrite or (not abt.leap_outfiles_exist("receptor")):
                abt.leap("receptor")
    
    else: # ligand
        if args.overwrite or (not abt.antechamber_outfiles_exist()):
            abt.antechamber()
        
        if args.leap:
            if args.overwrite or (not abt.leap_outfiles_exist("ligand")):
                abt.leap("ligand")
                
        if args.leap and args.receptor_pdb:
            if args.overwrite or (not abt.leap_outfiles_exist("complex")):
                abt.leap("complex")

    generated_output_paths = []
    for p in tempdir.glob("*.*"):
        if (args.leap and p.name.endswith(".inpcrd") or p.name.endswith(".prmtop")) or \
            ((not args.leap) and p.name.endswith("_ac.mol2") or p.name.endswith("_ac.frcmod")):
            target = args.workdir / p.name
            if not target.exists():
                p.rename(target)
                generated_output_paths.append(target)
            else:
                p.unlink()
        else:
            p.unlink()
    
    tempdir.rmdir()

    return generated_output_paths


"""
    Amber GB/SA input parameters

    igb - there are several "flavors" of GB available, depending upon the value of igb. The version
    that has been most extensively tested corresponds to igb=1; the "OBC" models (igb=2 and 5) are newer, but appear
    to give significant improvements and are recommended for most projects (certainly for peptides or proteins).
    The newest, most advanced, and least extensively tested model, GBn (igb=7), yields results in considerably better
    agreement with molecular surface Poisson-Boltzmann and explicit solvent results than the "OBC" models under
    many circumstances.

    PBRadii - choose various sets of atomic radii for generalized Born or Poisson-Boltzmann calculations.
        bondi       igb = 7             ref. 333    
        mbondi      igb = 1             ref. 198    
        mbondi2     igb = 2 or 5        ref. 182    
        mbondi3     igb = 8             ref. 25     
        amber6      igb = 1             ref. 190    only to be used for reproducing very early calc.
        
    reference - Amber20 manual, see Leap

    [333] A. Bondi. van der Waals volumes and radii. J. Phys. Chem., 1964, 68, 441–451.
    [198] V. Tsui; D. A. Case. Theory and applications of the generalized Born solvation model in macromolecular
            simulations. Biopolymers (Nucl. Acid. Sci.), 2001, 56, 275–291.
    [182] A. Onufriev; D. Bashford; D. A. Case. Exploring protein native states and large-scale conformational
            changes with a modified generalized Born model. Proteins, 2004, 55, 383–394.        
    [25] H. Nguyen; D. R. Roe; C. Simmerling. Improved Generalized Born Solvent Model Parameters for Protein
            Simulations. J. Chem. Theory Comput., 2013, 9, 2020–2034.
    [190] A. Okur; L. Wickstrom; C. Simmerling. Evaluation of salt bridge structure and energetics in peptides using
            explicit, implicit and hybrid solvation models. J. Chem. Theory Comput., 2008, 4, 488–498.

    In standard md in explicit solvent, the radii were not used
    The radii is for post-processing ONLY for both GB and PB
    You may change the radii after MD simulations by ParmEd

    parmed.py -p [old_radii.prmtop] -i input_file

    input_file:
        changeRadii mbondi3
        parmout new_radii.prmtop

    reference - http://archive.ambermd.org/201604/0297.html
"""

def run_mmpbsa_py(args:argparse.Namespace, 
                  idx:int, 
                  complex_dcd:pathlib.Path) -> MMPBSA_py_result | None:
    # dcd example: ../openmm/openmm_md_job_mol-0376-complex_solvated-out.dcd
    ligand = complex_dcd.name.replace(f"{args.dcd_prefix}_", "").replace(f"-complex_{args.dcd_suffix}.dcd", "")
    # ligand example: mol-0376
    receptor = args.receptor_name
    
    mdtraj_dir = complex_dcd.parent.resolve()
    outfile = f"mmgbsa-{ligand}-{args.method}-out.dat"

    tempdir = args.workdir / f"temp-{idx:04d}"
    tempdir.mkdir()
    
    if (args.workdir / outfile).exists() and (not args.overwrite):
        return None

    # check required files
    required_files = [
            args.prmtop_dir / f"{ligand}.prmtop",
            args.prmtop_dir / f"{receptor}.prmtop",
            args.prmtop_dir / f"{ligand}-complex.prmtop",
            args.prmtop_dir / f"{ligand}-complex_solvated.prmtop",
            complex_dcd,
        ]
    if args.method == 3: # additional files
        required_files += [
            args.prmtop_dir / f"{receptor}_solvated.prmtop",
            args.prmtop_dir / f"{ligand}_solvated.prmtop",
            mdtraj_dir / f"{args.dcd_prefix}_{receptor}_{args.dcd_suffix}.dcd",
            mdtraj_dir / f"{args.dcd_prefix}_{ligand}_{args.dcd_suffix}.dcd",
        ]

    if not all([ p.exists() for p in required_files ]):
        return None
    
    with open((tempdir / "mmgbsa.in").as_posix(), "w") as f:
        f.write(f"&general\n")
        f.write(f"  startframe={args.start_frame}, endframe={args.end_frame}, interval={args.interval},\n")
        f.write(f"  verbose=2, keep_files=0,\n")
        f.write(f"/\n")
        f.write(f"&gb\n")
        f.write(f"  igb={args.igb},\n")
        f.write(f"  saltcon={args.salt},\n")
        f.write("/\n")
    
    cmd = [ "MMPBSA.py", 
            "-O", 
            "-i", "mmgbsa.in", 
            "-o", outfile,
            "-cp", (args.prmtop_dir / f"{ligand}-complex.prmtop").as_posix(), 
            "-lp", (args.prmtop_dir / f"{ligand}.prmtop").as_posix(), 
            "-rp", (args.prmtop_dir / f"{receptor}.prmtop").as_posix(),
            "-sp", (args.prmtop_dir / f"{ligand}-complex_solvated.prmtop").as_posix(),
            "-y",  complex_dcd.as_posix(), 
            ]
        
    if args.method == 3: # additional arguments
        cmd += [
            "-srp", (args.prmtop_dir / f"{receptor}_solvated.prmtop").as_posix(),
            "-slp", (args.prmtop_dir / f"{ligand}_solvated.prmtop").as_posix(),
            "-yr", (mdtraj_dir / f"{args.dcd_prefix}_{receptor}_{args.dcd_suffix}.dcd").as_posix(),
            "-yl", (mdtraj_dir / f"{args.dcd_prefix}_{ligand}_{args.dcd_suffix}.dcd").as_posix(),
            ]

    p = subprocess.run(cmd, 
                       cwd=tempdir.as_posix(), 
                       capture_output=True) # redirect stdout and stderr

    generated_output_path = None
    for p in tempdir.glob("*.*"):
        if p.name.endswith(".dat"):
            target = args.workdir / p.name
            if not target.exists():
                p.rename(target)
                generated_output_path = target
            else:
                p.unlink()
        else:
            p.unlink()
    tempdir.rmdir()

    if generated_output_path is None:
        return None
    
    with generated_output_path.open() as f:
        for line in f:
            if line.startswith("DELTA TOTAL"):
                (dG, stdev, sem) = line.strip().split()[2:]

    return MMPBSA_py_result(ligand, 
                     float(dG), 
                     float(stdev), 
                     float(sem), 
                     generated_output_path)