import io
import os
import sys
import argparse
import numpy as np
import pandas as pd
import warnings

from pathlib import Path
from typing import List, Optional


try:
    from Bio import BiopythonWarning
    from Bio.PDB import PDBParser, Superimposer
    from Bio.PDB.PDBIO import PDBIO
    from Bio.PDB.Structure import Structure
    from Bio.PDB.StructureBuilder import StructureBuilder
    warnings.simplefilter('ignore', BiopythonWarning)
except:
    print("requires biopython")
    print()
    print("    $ pip install biopython")
    print()
    sys.exit(0)


try:
    import pypdb
    # A Python 3 toolkit for performing searches with the RCSB Protein Data Bank (PDB)
except:
    print("requires pypdb")
    print()
    print("    $ pip install pypdb")
    print()
    sys.exit(0)

protein_3_1 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
    "GLU":"E", "GLN":"Q", "GLY":"G", "HIS":"H", "ILE":"I",
    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
    "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
    }

protein_1_3 = { 
    "G":"GLY", "A":"ALA", "V":"VAL", "L":"LEU", "I":"ILE",
    "M":"MET", "P":"PRO", "F":"PHE", "W":"TRP", "S":"SER",
    "T":"THR", "N":"ASN", "Q":"GLN", "Y":"TYR", "C":"CYS",
    "K":"LYS", "R":"ARG", "H":"HIS", "D":"ASP", "E":"GLU",
    }

nucleic_1_3 = {
    'A':'ADE', 'G':'GUA', 'C':'CYT', 'U':'URA', 'T':'THY',
    }

nucleic_3_1 = {
    'ADE':'A', 'GUA':'G', 'CYT':'C', 'URA':'U', 'THY':'T',
    }

GlcNAC = {
    "HN2" : "H2N",
    "C7"  : "C2N",
    "O7"  : "O2N",
    "C8"  : "CME",
    "H81" : "H1M",
    "H82" : "H3M",
    "H83" : "H2M",
    "HO3" : "H3O",
    "HO4" : "H4O",
    "H61" : "H62",
    "H62" : "H61",
    "HO6" : "H6O",
    }


class PDBfile:
    def __init__(self, pdbfile:str | Path) -> None:
        self.parser = PDBParser()
        self.st = self.parser.get_structure("USER_PDB", pdbfile)
        # atomic coordinates
        coords = []
        for model in self.st:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        coords.append(atom.get_coord())
        self.coords = np.array(coords)
        self.shape = self.coords.shape
        

    def center_of_geometry(self) -> np.ndarray:
        return np.mean(self.coords, axis=0)
    
    
    def update_coord(self, coords:np.ndarray) -> None:
        assert self.coords.shape == coords.shape, "Matrix dimensions must match"
        # update atomic coordinates
        atom_idx = 0
        for model in self.st:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom.coord = coords[atom_idx,:]
                        atom_idx += 1


    def transform(self,
                  centroid: np.ndarray = np.array([0., 0., 0.]),
                  rot:np.ndarray = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]),
                  trans:np.ndarray = np.array([0., 0., 0.])) -> None:
        try:
            assert rot.shape == (3,3)
        except:
            raise ValueError('shape of rotation matrix should be (3,3)')
        try:
            assert trans.shape == (3,)
        except:
            raise ValueError('shape of translation matrix should be (3,)')
        
        # rotation
        self.coords = np.dot(self.coords - centroid, rot.T) + centroid
        
        # translation
        self.coords = self.coords + trans
        
        # update atomic coordinates
        atom_idx = 0
        for model in self.st:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom.coord = self.coords[atom_idx,:]
                        atom_idx += 1


    def get_seqeunce(self) -> None:
        wrap= 50
        for model in self.st:
            prev_residue_number = None
            header = "Model "+ str(model.id)
            sequence = ""
            for chain in model:
                print("chain=", chain)
                for residue in chain:
                    resname = residue.get_resname()
                    if not resname in aminoacids:
                        continue
                    hetflag,resseq,iCode = residue.get_id()
                    resSeq = int(resseq)
                    if prev_residue_number is None:
                        print("(first residue) %4d" % resSeq)
                    else:
                        num_missing_residues = resSeq-prev_residue_number-1
                        if num_missing_residues > 0:
                            sequence += "-" * num_missing_residues
                            print("(missing residues) %4d - %4d (%d)" %
                                (prev_residue_number+1,resSeq-1,num_missing_residues))
                    sequence += aminoacids[resname]
                    prev_residue_number = resSeq
                print("(last residue) %4d" % resSeq)
                # output
                print(header)
                l= len(sequence)//wrap
                for i in range(0,l):
                    start = wrap*i
                    end = max(wrap*(i+1),l)
                    print(sequence[start:end])


    @staticmethod
    def readlines(pdbfile:str):
        """(Not Used) Read PDB Lines.

        Args:
            pdbfile (str): filename.
        """
        pdblines = {1: []}
        with open(pdbfile, "r") as f:
            for line in f:
                if line.startswith("MODEL"):
                    model_number = int(line.strip().split()[1])
                    pdblines[model_number] = []

                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                    # PDB format version 2.3
                    serial = line[6:12]
                    name = line[12:16].strip()
                    altLoc = line[16:17]
                    if altLoc=="":
                        altLoc=" "
                    resName = line[17:21].strip()
                    chainId = line[21:22]
                    if chainId=="":
                        chainId=" "
                    resSeq = int(line[22:26])
                    iCode = line[26:27]
                    if iCode=="":
                        iCode=" "
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    occupancy = line[54:60]
                    tempFactor = float(line[60:66])
                    segId = line[72:76]
                    element = line[76:78]
                    charge = line[78:80]


    def write(self, filename: str | Path, 
              model:Optional[List[int]] = None,
              chain:Optional[List[str]] = None, 
              resname:Optional[List[str]] = None,
              resnumber:Optional[str] = None,
              ) -> None:
        """Write to PDB with selections.

        Examples:
            >>> pdb.write(model=[0], chain=['A'])
            write chain A and residues first-10,22-50,130-200,550,600-last
            >>> pdb.write(chain=['A'], resnumber="-10,22-50,130-200,550,600-")

        Args:
            filename (str | Path): output filename or path
            model (Optional[List[int]], optional): list of model numbers. Defaults to None.
            chain (Optional[List[str]], optional): list of chain ids. Defaults to None.
            resname (Optional[List[str]], optional): list of residue names. Defaults to None.
            resnumber (Optional[str], optional): residue number ranges. Defaults to None.
        """
        
        io = PDBIO()

        if (model is None) and (chain is None) and (resname is None) and (resnumber is None):
            # write structure as it as
            io.set_structure(self.st)
        else: 
            # write only select model(s), chain(s) or residue(s)
            # select residues by numbers
            if resnumber is not None:
                select_resseqs = [tuple(map(lambda s: int(s) if s else -1, r.split("-"))) for r in resnumber.split(",")]
                # for resnumber="-10,22-50,130-200,550,600-"
                # [(-1, 10), (22, 50), (130, 200), (550,), (600, -1)]
            builder = StructureBuilder()
            builder.init_structure(self.st.id)
            for _model in self.st:
                if (model is not None) and (_model.id not in model):
                    continue # reject model
                builder.init_model(_model.id, serial_num=None)
                for _chain in _model:
                    if (chain is not None) and (_chain.id not in chain):
                        continue # reject chain
                    builder.init_chain(_chain.id)
                    for _residue in _chain:
                        hetflag, resseq, iCode = _residue.get_id()
                        if (resname is not None) and (_residue.resname not in resname):
                            continue # reject residue by name
                        # select residue numbers
                        include_flags = []
                        for lu in select_resseqs:
                            if len(lu) == 1:
                                if lu == resseq:
                                    include_flags.append(True)
                                else:
                                    include_flags.append(False)
                            else:
                                (l,u) = lu
                                if (l == -1 or l <= resseq) and (u ==-1 or u >= resseq):
                                    include_flags.append(True)
                                else:
                                    include_flags.append(False)
                        if not any(include_flags):
                            continue # reject residue by number
                        builder.init_residue(_residue.resname, hetflag, resseq, iCode)
                        for atom in _residue:
                            builder.init_atom(atom.name, 
                                            atom.coord, 
                                            atom.bfactor, 
                                            atom.occupancy,
                                            atom.altloc,
                                            atom.fullname,
                                            None, # serial_number
                                            atom.element,
                                            )
            io.set_structure(builder.get_structure())
        with open(filename, "w") as f:
            io.save(f)
        


    def reorder(self) -> None:
        builder = StructureBuilder()
        builder.init_structure(self.st.id)
        for model in self.st:
            builder.init_model(model.id, serial_num=None)
            for chain in sorted(model, key=lambda x:x.id):
                builder.init_chain(chain.id)
                for residue in sorted(chain, key=lambda x: x.get_id()[1]):
                    hetflag, resseq, iCode = residue.get_id()
                    builder.init_residue(residue.resname, hetflag, resseq, iCode)
                    for atom in residue:
                        builder.init_atom(atom.name, 
                                        atom.coord, 
                                        atom.bfactor, 
                                        atom.occupancy,
                                        atom.altloc,
                                        atom.fullname,
                                        None, # serial_number
                                        atom.element,
                                        )
        self.st = builder.get_structure()



    def rename(self, 
               chain:Optional[dict] = None, 
               resname:Optional[dict] = None, 
               atom:Optional[dict] = None) -> None:
        """Rename PDB chains/residues/atoms.

        Examples:
            Rename chain 'C' to 'A' and 'D' to 'B'.
            >>> rename(chain={'C':'A', 'D':'B'})
            
            Rename residue 'UNL' to 'LIG' for all chains.
            >>> rename(resname={'UNL' : 'LIG'})

            Rename residue 'UNL' to 'LIG' only for chain B.
            >>> rename(chain={"B":"B"}, resname={'UNL' : 'LIG'})

            Rename atom '2H2' to 'H22' for all residues and chains.
            >>> rename(atom={"2H2": "H22"})

            Rename atoms 'H1' and 'H2' to 'H11' and 'H12' for chain C and residue UNL.
            >>> rename(chain={"C:C"}, resname={"UNL":"UNL"}, atom={"H1" : "H11", "H2": "H12"})

        Args:
            chain (Optional[dict], optional): map of chain ids {old:new}. Defaults to None.
            resname (Optional[dict], optional): map of residue names {old:new}. Defaults to None.
            atom (Optional[dict], optional): map of atom names {old:new}. Defaults to None.
        """
        for model in self.st:
            for _chain in model:
                if chain is None:
                    subject_chains = [c.id for c in model.get_chains()]
                else:
                    subject_chains = [k for k in chain]
                
                if _chain.id not in subject_chains:
                    continue
                
                if chain:
                    old_chain_id = _chain.id
                    new_chain_id = chain[old_chain_id]
                    _chain.id = new_chain_id
                    if old_chain_id != new_chain_id:
                        print(f"renamed chain id : ", end=" ")
                        print(f"{old_chain_id}", end=" -> ")
                        print(f"{new_chain_id}")
                        
                for _residue in _chain:
                    if resname is None:
                        subject_resnames = [r.resname for r in _chain.get_residues()]
                    else:
                        subject_resnames = [k for k in resname]
                    
                    if _residue.resname not in subject_resnames:
                        continue

                    hetflag, resseq, iCode = _residue.get_id()
                    
                    if resname:
                        old_resname = _residue.resname
                        new_resname = resname[old_resname]
                        _residue.resname = new_resname
                        if old_resname != new_resname:
                            print(f"renamed residue {_chain.id}", end=" : ")
                            print(f"{old_resname}({resseq})", end=" ->")
                            print(f"{new_resname}({resseq})")
                    else:
                        old_resname = None
                        new_resname = None
                        
                    for _atom in _residue:
                        if atom is None:
                            subject_atoms = [a.name for a in _residue.get_atoms()]
                        else:
                            subject_atoms = [k for k in atom]
                        
                        if _atom.name not in subject_atoms:
                            continue

                        if atom:
                            old_atom_name = _atom.name
                            new_atom_name = atom[old_atom_name]
                            _atom.name = new_atom_name
                            if old_atom_name != new_atom_name:
                                print(f"renamed atom {_chain.id}.{_residue.resname}({resseq})", end=" : ")
                                print(f"{old_atom_name}", end=" -> ")
                                print(f"{new_atom_name}")


    def reorient(self, 
                 residue_selection:str="", 
                 invert_Z:bool=False, 
                 offset:np.ndarray=np.array([0.,0.,0.])) -> None:
        """Reorient coordinates according to the principal axis.

        Examples:
            Chain A residues 12-50 and chain B residues 100-200 will be used for principal axis calculation.
            >>> pdb.reorient("A:12-50,B:100-200") 

        Args:
            residue_selection (str, optional): residues for principal axis calculation. Defaults to "".
            invert_Z (bool, optional): whether to invert Z axis. Defaults to False.
            offset (np.ndarray, optional): translate coordinates. Defaults to np.array([0.,0.,0.]).
        """
        coords = []
        if residue_selection:
            subset_residues = {}
            # subset of atoms for Eigenvector/Eigenvalue calculations
            for chain in residue_selection.split(","):
                (chainId, resSeq_range) = chain.split(":")
                resSeq_range_tuple = tuple(sorted(map(int, resSeq_range.split("-"))))
                if chainId in subset_residues:
                    subset_residues[chainId].append(resSeq_range_tuple)
                else:
                    subset_residues[chainId] = [resSeq_range_tuple]
            for model in self.st:
                for chain in model:
                    if chain.id not in subset_residues:
                        continue
                    for residue in chain:
                        hetflag, resseq, iCode = residue.get_id()
                        flag = []
                        for (lb, ub) in subset_residues[chain.id]:
                            if resseq >= lb and resseq <= ub:
                                flag.append(True)
                            else:
                                flag.append(False)
                        if any(flag):
                            coords.extend([atom.coord for atom in residue])
        else:
            for model in self.st:
                for chain in model:
                    for residue in chain:
                        coords.extend([atom.coord for atom in residue])

        coords = np.array(coords, float)
        centroid = np.mean(coords, axis=0)
        box_min = np.min(coords, axis=0)
        box_max = np.max(coords, axis=0)
        print("REMARK geometric center:", centroid)
        print("REMARK coordinate range:", box_min, box_max)
        print("REMARK box size:", box_max-box_min)

        coords = coords - centroid
        
        # import numpy as np
        # data = np.array([[1, 2], [3, 4], [5, 6]])  # Example data
        # cov_matrix = np.cov(data, rowvar=False)
        # eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
        # rotation_matrix = eigenvectors
        # rotated_data = np.dot(data, rotation_matrix.T)

        inertia = np.dot(coords.transpose(), coords)
        w,v = np.linalg.eig(inertia)

        # warning eigen values are not necessary ordered!
        # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html
        #--------------------------------------------------------------------------
        # order eigen values (and eigen vectors)
        #
        # axis1 is the principal axis with the biggest eigen value (eval1)
        # axis2 is the principal axis with the second biggest eigen value (eval2)
        # axis3 is the principal axis with the smallest eigen value (eval3)
        #--------------------------------------------------------------------------
        # axis1 --> Z
        # axis2 --> Y
        # axis3 --> X
        order = np.argsort(w)
        eval3, eval2, eval1 = w[order]
        axis3, axis2, axis1 = v[:, order]

        print("REMARK x-axis",axis3,"eval=",eval3)
        print("REMARK y-axis",axis2,"eval=",eval2)
        print("REMARK z-axis",axis1,"eval=",eval1)

        R = np.array([axis3, axis2, axis1]) # decreasing order
        # R_inv = np.linalg.inv(R)

        self.transform(rot=R, trans=-centroid)
        # xyz = np.array([x,y,z])
        # xyz = xyz - center
        # xyz = np.dot(R_inv, xyz)
        # x, y, z = xyz
        # if args.inverse_z_axis:
        #     z = -z
        # x += args.offset_x
        # y += args.offset_y
        # z += args.offset_z


    @staticmethod
    def expand_residues(txt):
        residues = []
        step1 = txt.split(",")
        for item in step1:
            step2 = item.split("-")
            if len(step2) == 2:
                residues += range(int(step2[0]), int(step2[1])+1)
            elif len(step2) == 1:
                residues += [ int(step2[0]) ]
        return sorted(residues)


    def contacts(self, 
                       residue_selection_1:str, 
                       residue_selection_2:str,
                       cutoff:float=4.0) -> None:
        residues_1 = PDBfile.expand_residues(residue_selection_1)
        residues_2 = PDBfile.expand_residues(residue_selection_2)

        source_atoms = []
        target_atoms = []
        for model in self.st:
            for chain in model:
                for residue in chain:
                    hetflag, resSeq, iCode = residue.get_id()
                    resName = residue.get_resname()
                    if resSeq in residues_1:
                        source_atoms += [ atom for atom in residue ]
                    elif resSeq in residues_2:
                        target_atoms += [ atom for atom in residue ]

        close_contacts = []
        for ai in source_atoms:
            for aj in target_atoms:
                d = ai-aj
                if d < cutoff:
                    close_contacts.append((ai, aj, d))
        close_contacts = sorted(close_contacts, key=lambda x: (x[0], x[1], x[2]))
        for (ai, aj, d) in close_contacts:
            i_hetflag, i_resSeq, i_iCode = ai.get_parent().get_id()
            j_hetflag, j_resSeq, j_iCode = aj.get_parent().get_id()
            print("%4d %4s %4d %4s %3.1f" % (i_resSeq, ai.id, j_resSeq, aj.id, d))



def renumber(residue:Optional[dict]):
    pass


def download(PDBID:str, verbose:bool=True) -> None:
    pdb = pypdb.get_pdb_file(PDBID, filetype="pdb")
    filename = PDBID+".pdb"
    if verbose:
        print("downloading ",PDBID,"as", filename, end="")
    with open(filename,"w") as f:
        f.write(pdb)
        if verbose:
            print("done")


def get_ligands(PDBID:str) -> None:
    """
    Examples:

        @structureId : 4KEB
        @chemicalID : 1QZ
        @type : non-polymer
        @molecularWeight : 423.51
        chemicalName : 6-ethyl-5-{(3S)-3-[3-(isoquinolin-5-yl)-5-methoxyphenyl]but-1-yn-1-yl}pyrimidine-2,4-diamine
        formula : C26 H25 N5 O
        InChI : InChI=1S/C26H25N5O/c1-4-24-23(25(27)31-26(28)30-24)9-8-16(2)18-12-19(14-20(13-18)32-3)21-7-5-6-17-15-29-11-10-22(17)21/h5-7,10-16H,4H2,1-3H3,(H4,27,28,30,31)/t16-/m1/s1
        InChIKey : MGLLCDAARSVGLO-MRXNPFEDSA-N
        smiles : CCc1c(c(nc(n1)N)N)C#C[C@@H](C)c2cc(cc(c2)OC)c3cccc4c3ccnc4

        @structureId : 4KEB
        @chemicalID : FOL
        @type : non-polymer
        @molecularWeight : 441.397
        chemicalName : FOLIC ACID
        formula : C19 H19 N7 O6
        InChI : InChI=1S/C19H19N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)7-21-10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,8,12,21H,5-7H2,(H,24,29)(H,27,28)(H,31,32)(H3,20,22,25,26,30)/t12-/m0/s1
        InChIKey : OVBPIULPVIDEAO-LBPRGKRZSA-N
        smiles : c1cc(ccc1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCc2cnc3c(n2)C(=O)N=C(N3)N
              
        @structureId : 4KEB
        @chemicalID : NDP
        @type : non-polymer
        @molecularWeight : 745.421
        chemicalName : NADPH DIHYDRO-NICOTINAMIDE-ADENINE-DINUCLEOTIDE PHOSPHATE
        formula : C21 H30 N7 O17 P3
        InChIKey : ACFIXJIJDZMPPO-NNYOXOHSSA-N
        InChI : InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1
        smiles : c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@@](=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)N5C=CCC(=C5)C(=O)N)O)O)O)OP(=O)(O)O)N
    """

    ligandInfo = pypdb.get_ligands(PDBID)
    
    try:    
        ligands = ligandInfo["ligandInfo"]["ligand"]
    except:
        ligands = []
    
    for ligand_dict in ligands:
        for k in ligand_dict:
            print("%20s : %s" % (k, ligand_dict[k]))

    # chem_desc= pypdb.describe_chemical(PDBID)
    # chem_desc["describeHet"]["ligandInfo"]["ligand"]["smiles"]
    # pdb_desc= pypdb.describe_pdb("5cgc")
    # pdb content as list of strings
    data = {
        "name": [],
        "smiles": [],
        "molwt": [],
        "chemical name": [],
        "InChiKey": [],
        "InChi": [],
    }
    if "," in PDBID:
        PDBID_list = PDBID.split(",")
        for PDBID in PDBID_list:
            for d in ligands:
                try:
                    data["name"].append(d['@chemicalID'])
                    data["smiles"].append(d['smiles'])
                    data["molwt"].append(float(d['@molecularWeight']))
                    data["chemical name"].append(d['chemicalName'])
                    data["InChiKey"].append(d['InChIKey'])
                    data["InChi"].append(d['InChi'])
                except:
                    continue
    else:
        for d in ligands:
            try:
                data["name"].append(d['@chemicalID'])
                data["smiles"].append(d['smiles'])
                data["molwt"].append(float(d['@molecularWeight']))
                data["chemical name"].append(d['chemicalName'])
                data["InChiKey"].append(d['InChIKey'])
                data["InChi"].append(d['InChi'])
            except:
                continue





def reorient():
    argparser = argparse.ArgumentParser(description='Reorient PDB')
    argparser.add_argument("filename", help="pdb filename")
    argparser.add_argument("--inverse-z", dest="inverse_z_axis", help="inverse z axis", default=False, action="store_true")
    argparser.add_argument("--residue", help="residues for principal axes calc. ex. A:23-50,B:1-90", default="")
    argparser.add_argument("--offset-x", dest="offset_x", help="translate x coordinates", type=float, default=0.0)
    argparser.add_argument("--offset-y", dest="offset_y", help="translate y coordinates", type=float, default=0.0)
    argparser.add_argument("--offset-z", dest="offset_z", help="translate z coordinates", type=float, default=0.0)
    argparser.add_argument("--segid", dest="segId", help="override segment id", default="")
    argparser.add_argument("--chainid", dest="chainId", help="override chain id", default="")
    args = argparser.parse_args()

    # subset of atoms for Eigenvector/Eigenvalue calculations
    subset_atoms = []
    if args.residue:
        for chain in args.residue.split(","):
            (chainId, resSeq_range) = chain.split(":")
            (resSeq_begin, resSeq_end) = resSeq_range.split("-")
            subset_atoms.append((chainId, int(resSeq_begin), int(resSeq_end)))

    with open(args.filename, "r") as pdb_file:
        xyz = []
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chainId = line[21:22]
                resSeq = int(line[22:26])
                if subset_atoms:
                    for chainId_, resSeq_begin, resSeq_end in subset_atoms:
                        if not (chainId == chainId_ and resSeq >= resSeq_begin and resSeq <= resSeq_end):
                            continue
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                xyz.append([x,y,z])

        print("REMARK reoriented according to its principal axes")
        if args.residue:
            print(f"REMARK principal axes are calculated using --residue {args.residue}")
        print(f"REMARK coordinates from {args.filename} ({len(xyz)} atoms)")
        if args.inverse_z_axis:
            print("REMARK Z axis inverted")
        print(f"REMARK coordinates offset X {args.offset_x:8.3f} Y {args.offset_x:8.3f} Z {args.offset_x:8.3f}")

        coord = np.array(xyz, float)

        # geometric center
        center = np.mean(coord, 0)
        box_min = np.min(coord, 0)
        box_max = np.max(coord, 0)
        print("REMARK geometric center:", center)
        print("REMARK coordinate range:", box_min, box_max)
        print("REMARK box size:", box_max-box_min)

        coord = coord - center

        # compute principal axis matrix
        inertia = np.dot(coord.transpose(), coord)
        w,v = np.linalg.eig(inertia)

        # warning eigen values are not necessary ordered!
        # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html
        #--------------------------------------------------------------------------
        # order eigen values (and eigen vectors)
        #
        # axis1 is the principal axis with the biggest eigen value (eval1)
        # axis2 is the principal axis with the second biggest eigen value (eval2)
        # axis3 is the principal axis with the smallest eigen value (eval3)
        #--------------------------------------------------------------------------
        # axis1 --> Z
        # axis2 --> Y
        # axis3 --> X
        order = np.argsort(w)
        eval3, eval2, eval1 = w[order]
        axis3, axis2, axis1 = v[:, order]

        print("REMARK x-axis",axis3,"eval=",eval3)
        print("REMARK y-axis",axis2,"eval=",eval2)
        print("REMARK z-axis",axis1,"eval=",eval1)

        R = np.array([axis3, axis2, axis1]) # decreasing order
        R_inv = np.linalg.inv(R)

        pdb_file.seek(0)
        for line in pdb_file:
            if not line: continue
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                xyz = np.array([x,y,z])
                xyz = xyz - center
                xyz = np.dot(R_inv, xyz)
                x, y, z = xyz
                if args.inverse_z_axis:
                    z = -z
                x += args.offset_x
                y += args.offset_y
                z += args.offset_z

                # keep the rest of information as they are
                line = line[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + line[54:]

                # override chain id
                if args.chainId:
                    line = line[:21] + "%1s" % args.chainId + line[22:]

                # override segment id
                if args.segId:
                    line = line[:72] + "%4s" % args.segId + line[76:]

                print(line, end='')
            else:
                print(line, end='')


