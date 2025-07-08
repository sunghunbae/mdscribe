import sys
import numpy as np
import pandas as pd
import yaml
import parmed

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from pathlib import Path

from mdscribe.helper import PDBfile, kabsch_algorithm, get_closest_distance, set_apart, view_svg

from rdkit import Chem
from rdkit.Chem import rdFMCS


@dataclass
class Molecule:
    """Class for protein information."""
    st: parmed.structure.Structure
    idx: List[int] = field(default_factory=list)
    
    def coor(self) -> np.ndarray:
        return self.st.coordinates[self.idx,:]

    def centroid(self) -> np.ndarray:
        return np.mean(self.coor(), axis=0)


@dataclass
class Pocket:
    """Class for pocket information."""
    st: parmed.structure.Structure
    idx: List[int] = field(default_factory=list)
    residues: List = field(default_factory=list)
    imap: Dict[int,int] = field(default_factory=dict) # pocket atom idx -> st atom idx
    
    def coor(self) -> np.ndarray:
        return self.st.coordinates[self.idx,:]
    
    def centroid(self) -> np.ndarray:
        return np.mean(self.coor(), axis=0)
    
    def name(self) -> List[str]:
        return [self.st.atoms[i].name for i in self.idx ]
    
    def residue_number(self) -> List[int]:
        return [self.st.atoms[i].residue.number for i in self.idx ]
        # return [ r.number for r in self.residues ]
    
    def residue_name(self) -> List[str]:
        return [self.st.atoms[i].residue.name for i in self.idx ]
        # return [ r.name for r in self.residues ]


@dataclass
class Ligand:
    """Class for ligand information."""
    st: parmed.structure.Structure | None = None
    idx0: List[int] = field(default_factory=list) # all atoms
    idx: List[int] = field(default_factory=list) # non-hydrogens
    imap: Dict[int,int] = field(default_factory=dict) # pdb atom idx -> st atom idx
    mol: Chem.Mol | None = None
    highlight: List[int] = field(default_factory=list)
    
    def coor(self) -> np.ndarray:
        return self.st.coordinates[self.idx,:]
    
    def centroid(self) -> np.ndarray:
        return np.mean(self.coor(), axis=0)
    
    def name(self) -> List[str]:
        return [self.st.atoms[i].name for i in self.idx ]
    
    def atomic_number(self) -> List[int]:
        return [self.st.atoms[i].atomic_number for i in self.idx ]
    
    def residue_name(self) -> List[str]:
        return [self.st.atoms[i].residue.name for i in self.idx ]
    
    def residue_number(self) -> List[int]:
        return [self.st.atoms[i].residue.number for i in self.idx ]


@dataclass
class Reference:
    """Class for reference information."""
    st: parmed.structure.Structure
    pocket: Pocket
    ligand: Ligand
    
    def coor(self) -> np.ndarray:
        return self.st.coordinates
    
    def centroid(self) -> np.ndarray:
        return np.mean(self.coor(), axis=0)


class PreTernary:
    """Staging proteins and ligand for ternary complex modeling."""

    def __init__(self,
                 config: Dict | str | Path,
                 workdir: str | Path = '.',
                ):
       
        if isinstance(config, str):
            self.config = yaml.safe_load(Path(config).read_text())
        elif isinstance(config, Path):
            self.config = yaml.safe_load(config.read_text())
        elif isinstance(config, dict):
            self.config = config

        self.references = self.config['references']

        if not isinstance(workdir, Path):
            workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)        
        self.workdir = workdir

        self.molecule = {}
        self.pocket = {}
        self.ligand = {}
        self.centroid = {}

        
        try:
            mol = Chem.MolFromSmiles(self.config['ligand']['smiles'])
            assert mol
            self.ligand['0'] = Ligand(mol=mol)
        except:
            print("Error: ligand cannot be defined")
            sys.exit(1)

        for k in self.references:
            
            try:
                self.ligand[k] = Ligand(mol=Chem.MolFromSmiles(self.references[k]['ligand']['smiles']))
                assert self.ligand[k]
            except:
                print(f"Error: invalid ligand SMILES for the reference {k}")
                sys.exit(1)
            
            try:
                mcs_res = rdFMCS.FindMCS([self.ligand['0'].mol, self.ligand[k].mol])
                mcs_mol = Chem.MolFromSmarts(mcs_res.smartsString)
                assert self.ligand['0'].mol.HasSubstructMatch(mcs_mol)
                self.ligand[k].highlight = self.ligand[k].mol.GetSubstructMatch(mcs_mol) # tuple
                self.ligand['0'].highlight += self.ligand['0'].mol.GetSubstructMatch(mcs_mol)
            except:
                print(f"Error: no common substructure between the ligand and the reference {k}")
                sys.exit(1)
            
            self.references[k]["smarts"] = mcs_res.smartsString
            
            if isinstance(self.references[k]["pdb"], str):
                pdb = parmed.load_file(self.references[k]["pdb"])
            elif isinstance(self.references[k]["pdb"], Path):
                pdb = parmed.load_file(self.references[k]["pdb"].as_posix())
            self.molecule[k] = Molecule(st = pdb)
            
            self.pocket[k] = self.create_pocket(self.molecule[k].st, 
                                                self.references[k]["ligand"]["resname"], 
                                                self.references[k]["pocket"]["distance"],
                                                self.references[k]["pocket"]["include"],
                                                self.references[k]["pocket"]["exclude"],
                                                )


    def molview(self, 
                ref_id:str | None = None, 
                width:int=600, 
                height:int=400) -> None:
        """Depict a ligand with substructure highlights.

        Args:
            ref_id (str | None, optional): id of reference structure. Defaults to None.
            width (int, optional): width. Defaults to 600.
            height (int, optional): height. Defaults to 400.
        """
        if ref_id is None:
            view_svg(self.ligand['0'].mol, 
                     highlight=self.ligand['0'].highlight, 
                     width=width, height=height)
        else:
            view_svg(self.ligand[ref_id].mol, 
                     highlight=self.ligand[ref_id].highlight,
                     width=width, height=height)


    def create_pocket(self, 
                      st : parmed.structure.Structure, 
                      resname : str, 
                      dist : float, 
                      include : str, 
                      exclude : str) -> Pocket:
        """Create a Pocket class object.

        Args:
            st (parmed.structure.Structure): input structure.
            resname (str): residue name for ligand.
            dist (float): distance cutoff from ligand which is regarded as pocket.
            include (str): residue or atom selection to include as pocket.
            exclude (str): residue or atom selection to exlude from pocket.

        Raises:
            ValueError: when pocket cannot be defined.

        Returns:
            Pocket: object containing pocket information.
        """
        try:
            # select all `pocket_atom` atoms within `pocket_dist` from the `ref_resname` residues and save to a pdb.
            # the whole residue is selected if any atom satisfies the distance criteria
            pocket = st[f"(:{resname}<:{dist})&(!:{resname})&({include})&!({exclude})"]
            imap = {(a.residue.number, a.name) : a.idx for a in st.atoms}
            indices = [imap[(b.residue.number, b.name)] for b in pocket.atoms]
            assert len(pocket.residues) > 0
            assert len(pocket.atoms) > 0
        except:
            raise ValueError("cannot define reference ligand binding pocket")

        return Pocket(st=st, idx=indices, residues=pocket.residues, imap=imap)
        

    def distance(self) -> None:
        """Summary of centroid distances.
        """
        for k in self.references:
            self.centroid[k] = {
                "pocket": self.pocket[k].centroid(),
            }
            print(f"Centroid({k}):")
            for kk, v in self.centroid[k].items():
                print(f"    {kk}: {v}")   
        [k1, k2] = [ k for k in self.references ]
        dpp = np.linalg.norm(self.centroid[k1]["pocket"] - self.centroid[k2]["pocket"])
        print(f"Distance Centroid({k1} pocket)-Centroid({k2} pocket) {dpp:8.3f}")
        cc = get_closest_distance(self.molecule[k1].st.coordinates, 
                                  self.molecule[k2].st.coordinates)
        print(f"Closest Distance: {cc:8.3f}")

        
    def stage(self, 
              outfile: str | Path,
              overwrite: bool = False,
              r: float=20.0,
              steric_clash_cutoff: float=3.0):
        """Place molecule for MD simulation.

        Args:
            outfile (str | Path): _description_
            overwrite (bool, optional): _description_. Defaults to False.
            r (float, optional): _description_. Defaults to 20.0.
            steric_clash_cutoff (float, optional): _description_. Defaults to 3.0.
        """
        
        (fixed, movable) = [k for k in self.molecule ]
        (trans, c, n) = set_apart(self.molecule[fixed].st.coordinates, 
                                  self.molecule[movable].st.coordinates,
                                  r, steric_clash_cutoff)
        print(f"Setting apart {movable} coordinates (r= {r:4.1f}):")
        print(f"    number of steric clashes (<{steric_clash_cutoff:3.1f}) : {n}")
        print(f"    closest distance   : {c:4.1f}")
        print(f"    translation vector : {trans}")
        
        pdb = PDBfile(self.references[movable]["pdb"])
        pdb.transform(trans=trans)
        if overwrite or \
            (isinstance(outfile, str) and Path(outfile).exists()) or \
            (isinstance(outfile, Path) and outfile.exists()):
            print(f"    writing translated coordinates to {outfile}")
            pdb.write(outfile)
        else:
            print(f"    skipping to write to {outfile} because it already exists.")



class Ternary:
    """Modeling a ternary complex."""

    def __init__(self, 
                 st:parmed.structure.Structure, 
                 config: dict | str | Path, 
                 workdir: str | Path = '.'):
        self.st = st
        self.imap = {(a.residue.number,a.name) : a.idx for a in self.st.atoms}
        if isinstance(config, str):
            self.config = yaml.safe_load(Path(config).read_text())
        elif isinstance(config, Path):
            self.config = yaml.safe_load(config.read_text())
        elif isinstance(config, dict):
            self.config = config
        if isinstance(workdir, Path):
            self.workdir = workdir
        else:
            self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)        

        self.references = self.config['references']
        self.molecule = {}
        for k, amber_mask in self.config['simulation']['topology'].items():
            idx = [self.imap[(b.residue.number, b.name)] for b in self.st[amber_mask].atoms]
            self.molecule[k] = Molecule(self.st, idx)

        self.map = {}
        self.map_data = {}
        self.ligand = {}
        self.pocket = {}
        self.ref = {}
        self.rmsd_ref_coor = np.zeros((len(self.st.atoms), 3))
        self.rmsd_idx_group = {}
        self.centroid = {}
        # reference coordinates for RMSDForce
        # should have identical number of coordinates as the self.st
        # even though all the coordinates are not used

        self.ligand['0'] = self.create_ligand(self.st, 
                                              self.config['ligand']['resname'], "structure-ligand")
        
        for k in self.references:
            self.map_to_reference(k)
            self.rmsd_ref_coor[np.array(self.pocket[k].idx),:] = self.ref[k].pocket.coor()
            self.rmsd_ref_coor[np.array(self.ligand[k].idx),:] = self.ref[k].ligand.coor()
            self.rmsd_idx_group[k] = np.concatenate((self.pocket[k].idx, self.ligand[k].idx,))
            self.map[k] = {
                'ligand': pd.DataFrame(self.map_data[k]['ligand']),
                'pocket': pd.DataFrame(self.map_data[k]['pocket']),
            }

    
    def create_ligand(self, st:parmed.structure.Structure, resname:str, prefix:str) -> Ligand:
        """Create a Ligand class object.

        Args:
            st (parmed.structure.Structure): input structure.
            resname (str): ligand residue name.
            prefix (str): prefix for output pdb file.

        Raises:
            ValueError: when residue name does not exist.

        Returns:
            Ligand: object containing ligand information.
        """
        outfile = (self.workdir / f"{prefix}-{resname}.pdb").as_posix()
        try:
            idx0 = [atom.idx for atom in st.atoms if (atom.residue.name == resname)]
            assert len(idx0) > 0
        except:
            raise ValueError(f"{resname} does not exist")
        try:
            # write the ligand to a pdb file which then is read by Chem.MolFromPDBFile to create a rdkit.Chem.Mol
            # non-hydrogen atoms
            # Order of atoms in the pdb is same as that of the rdmol
            indices = []
            imap = {}
            pdb_idx = 0
            for atom in st.atoms:
                if atom.residue.name == resname and atom.atomic_number > 1:
                    imap[pdb_idx] = atom.idx
                    indices.append(atom.idx)
                    pdb_idx += 1
            st[f':{resname}&!@H='].save(outfile, overwrite=True)
            mol = Chem.MolFromPDBFile(outfile)
            assert mol
            assert imap
        except:
            raise ValueError(f"cannot create molecule for {resname}")
            
        return Ligand(st=st, idx0=idx0, idx=indices, imap=imap, mol=mol)


    def create_pocket(self, st:parmed.structure.Structure, resname:str, dist:float, include:str, exclude:str) -> Pocket:
        """Create a Pocket class object.

        Args:
            st (parmed.structure.Structure): input structure.
            resname (str): residue name for ligand.
            dist (float): distance cutoff from ligand which is regarded as pocket.
            include (str): residue or atom selection to include as pocket.
            exclude (str): residue or atom selection to exlude from pocket.

        Raises:
            ValueError: when pocket cannot be defined.

        Returns:
            Pocket: object containing pocket information.
        """
        try:
            # select all `pocket_atom` atoms within `pocket_dist` from the `ref_resname` residues and save to a pdb.
            # the whole residue is selected if any atom satisfies the distance criteria
            pocket = st[f"(:{resname}<:{dist})&(!:{resname})&({include})&!({exclude})"]
            imap = {(a.residue.number, a.name) : a.idx for a in st.atoms}
            indices = [imap[(b.residue.number, b.name)] for b in pocket.atoms]
            assert len(pocket.residues) > 0
            assert len(pocket.atoms) > 0
        except:
            raise ValueError("cannot define reference ligand binding pocket")

        return Pocket(st=st, idx=indices, residues=pocket.residues, imap=imap)
        

    def create_map_data(self, ref_id:str) -> None:
        """Create a map/dictionary for atom indexes.

        Args:
            ref_id (str): id for reference.
        """
        self.map_data[ref_id] = {
            'ligand' : {
                'ref_idx': self.ref[ref_id].ligand.idx,
                'ref_name': self.ref[ref_id].ligand.name(),
                'ref_atomic_number' : self.ref[ref_id].ligand.atomic_number(), 
                'ref_residue_name' : self.ref[ref_id].ligand.residue_name(),
                'ref_residue_number' :  self.ref[ref_id].ligand.residue_number(),
                'idx': self.ligand[ref_id].idx,
                'name' : self.ligand[ref_id].name(),
                'atomic_number': self.ligand[ref_id].atomic_number(),
                'residue_name' : self.ligand[ref_id].residue_name(),
                'residue_number' : self.ligand[ref_id].residue_number(),
                },
            'pocket' : { 
                'ref_idx' : self.ref[ref_id].pocket.idx,
                'ref_name' : self.ref[ref_id].pocket.name(),
                'ref_residue_name' : self.ref[ref_id].pocket.residue_name(),
                'ref_residue_number' : self.ref[ref_id].pocket.residue_number(),
                'idx' : self.pocket[ref_id].idx,
                'name' : self.pocket[ref_id].name(),
                'residue_name' : self.pocket[ref_id].residue_name(),
                'residue_number' : self.pocket[ref_id].residue_number(),
            }
        }
            
        
    def map_to_reference(self, ref_id: str, quiet:bool=False) -> None:
        """Map structure to the given reference.

        Args:
            ref_id (str): id of reference structure.
            quiet (bool, optional): whether to print details. Defaults to False.

        Raises:
            ValueError: when mapping fails.
        """
        pdb     = self.references[ref_id]["pdb"]
        resname = self.references[ref_id]["ligand"]["resname"]
        dist    = self.references[ref_id]["pocket"]["distance"]
        include = self.references[ref_id]["pocket"]["include"]
        exclude = self.references[ref_id]["pocket"]["exclude"]
        
        if not quiet:
            print(f"Mapping to {ref_id}...")
            
        try:
            if isinstance(pdb, Path):
                ref_st = parmed.load_file(pdb.as_posix())
            else:
                ref_st = parmed.load_file(pdb)
            assert ref_st
        except:
            raise ValueError(f"cannot load reference pdb file: {pdb}") 

        self.map_data[ref_id] = {'ligand': {}, 'pocket': {}}
        
        # Ligand object
        ref_ligand = self.create_ligand(ref_st, resname, "reference-ligand")

        # MCS search
        mcs_res = rdFMCS.FindMCS([self.ligand['0'].mol, ref_ligand.mol])
        query = Chem.MolFromSmarts(mcs_res.smartsString)

        try:
            assert self.ligand['0'].mol.HasSubstructMatch(query)
            assert ref_ligand.mol.HasSubstructMatch(query)
        except:
            print(f"No common substructure between the ligand and the reference {ref_id}")
            sys.exit(1)

        self.references[ref_id]["smarts"] = mcs_res.smartsString

        ligand_match = self.ligand['0'].mol.GetSubstructMatch(query)

        ref_ligand_match = ref_ligand.mol.GetSubstructMatch(query)

        # GetSubstructMatch() returns the indices of the molecule’s atoms 
        # that match a substructure query (only a single match is returned)
        # the ordering of the indices corresponds to the atom ordering in the query. 
        # For example, the first index is for the atom in this molecule 
        # that matches the first atom in the query.
        # create a ligand substructure that matches with the reference in SMARTS

        indices = [self.ligand['0'].imap[i] for i in ligand_match]

        # ligand substructure that matches with the reference
        self.ligand[ref_id] = Ligand(self.st, 
                                     idx=indices, 
                                     mol=self.ligand['0'].mol, 
                                     highlight=ligand_match)
        
        # complete the Ligand object attributes
        ref_ligand.idx = [ref_ligand.imap[i] for i in ref_ligand_match]
        ref_ligand.highlight = ref_ligand_match

        # reference Pocket object
        ref_pocket = self.create_pocket(ref_st, resname, dist, include, exclude)

        # create a Reference object
        self.ref[ref_id] = Reference(st=ref_st, ligand=ref_ligand, pocket=ref_pocket)

        # mapping structure to the reference
        # find matching residue numbers and names with reference
        ref_pocket_residue_numbers = [r.number for r in ref_pocket.residues]
        ref_pocket_residue_names = [r.name for r in ref_pocket.residues]
        offset = [ x-ref_pocket_residue_numbers[0] for x in ref_pocket_residue_numbers ]
        rmap = {} # ref_pocket_residue_number -> st residue.number
        residues = []
        for i in range(len(self.st.residues)):
            if self.st.residues[i].name == ref_pocket_residue_names[0]:
                found = True
                for o, n in zip(offset, ref_pocket_residue_names):
                    st_resname = self.st.residues[i+o].name
                    if st_resname in ['HIE', 'HID']:
                        st_resname = 'HIS'
                    if st_resname != n:
                        found = False
                        break
                if found: # sequence matches
                    for o, r, n in zip(offset, ref_pocket_residue_numbers, ref_pocket_residue_names):
                        st_residx = i + o
                        st_residue_name = self.st.residues[st_residx].name
                        st_residue_number = self.st.residues[st_residx].number
                        rmap[r] = st_residue_number
                        residues.append(self.st.residues[st_residx])
        try:
            assert len(rmap) > 0
        except:
            print("ref_pocket_residue_numbers", ref_pocket_residue_numbers)
            print("ref_pocket_residue_names", ref_pocket_residue_names)
            raise ValueError('No match found for the reference pocket.')
        
        indices = []
        for i in ref_pocket.idx:
            ai = ref_pocket.st.atoms[i]
            matching_residue_number = rmap[ai.residue.number]
            for aj in self.st.atoms:
                if aj.residue.number == matching_residue_number and aj.name == ai.name:
                    indices.append(aj.idx)

        # create a pocket that matches the reference pocket
        self.pocket[ref_id] = Pocket(self.st, idx=indices, residues=residues)

        try:
            assert len(self.ligand[ref_id].idx) == len(self.ref[ref_id].ligand.idx)
            assert len(self.pocket[ref_id].idx) == len(self.ref[ref_id].pocket.idx)
        except:
            raise ValueError('Numbers of indices between reference and structure do not match.')
            
        if not quiet:
            print(f"number of ligand atoms: {len(self.ligand[ref_id].idx)}")
            print(f"number of pocket atoms: {len(self.pocket[ref_id].idx)}")
        
        self.create_map_data(ref_id)

    
    def molview(self, ref_id:str | None = None, width:int=600, height:int=400) -> None:
        """Depict a ligand with substructure highlights.

        Args:
            ref_id (str | None, optional): id of reference structure. Defaults to None.
            width (int, optional): width. Defaults to 600.
            height (int, optional): height. Defaults to 400.
        """
        if ref_id is None:
            highlight = []
            for k in self.references:
                highlight.extend(self.ligand[k].highlight)
            view_svg(self.ligand['0'].mol, highlight=highlight, width=width, height=height)
        else:
            view_svg(self.ref[ref_id].ligand.mol, highlight=self.ref[ref_id].ligand.highlight)

    
    def distance(self) -> Dict[str, float]:
        """Summary of centroid distances.

        Returns:
            dict[str, float]: {ref_id: pocket-to-ligand distance}
        """
        pocket_to_ligand_dist = {}

        [k1, k2] = [ k for k in self.references ]

        for k in self.references:
            self.centroid[k] = {
                "ref_pocket" : self.ref[k].pocket.centroid(),
                "ref_ligand" : self.ref[k].ligand.centroid(),
                "pocket": self.pocket[k].centroid(),
                "ligand": self.ligand[k].centroid(),
            }
            
            print(f"Centroid({k}):")
            
            for kk, v in self.centroid[k].items():
                print(f"    {kk}: {v}")   
            
            dpls = np.linalg.norm(self.centroid[k]["pocket"] - self.centroid[k]["ligand"])
            dplr = np.linalg.norm(self.centroid[k]["ref_pocket"] - self.centroid[k]["ref_ligand"])
            
            print(f"Distance Centroid({k} pocket)-Centroid({k} ligand) {dpls:8.3f}", end=" ")
            print(f"(reference: {dplr:5.3f})")
            
            pocket_to_ligand_dist[k] = dplr

        dpp = np.linalg.norm(self.centroid[k1]["pocket"] - self.centroid[k2]["pocket"])
        dll = np.linalg.norm(self.centroid[k1]["ligand"] - self.centroid[k2]["ligand"])
        
        print(f"Distance Centroid({k1} pocket)-Centroid({k2} pocket) {dpp:8.3f}")
        print(f"Distance Centroid({k1} ligand)-Centroid({k2} ligand) {dll:8.3f}")

        return pocket_to_ligand_dist


    def rmsd(self) -> None:
        """Summary of RMSD."""
        for k in self.references:
            indices = self.rmsd_idx_group[k]
            (_rot, _trans, _centroid, _rmsd) = kabsch_algorithm(self.st.coordinates[indices,:], self.rmsd_ref_coor[indices,:]) 
            print(f"RMSD ({k} pocket and ligand) {_rmsd:8.3f}")


    def move_ligand_to(self, ref_id:str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get transformations to move ligand coordinates to the supposed binding pocket.

        Args:
            ref_id (str): id of reference structure (pocket).

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: (transformed coordinates, centroid, rotation matrix, translation vector)
        """
        # 1. align the reference to the structure using the reference pocket
        # ref_pocket_pos = self.ref_st[ref_id][self.ref_pocket_idx_group[ref_id],:]
        # pocket_pos = self.st.coordinates[self.pocket_idx_group[ref_id],:]
        (rot, trans, centroid, rmsd) = kabsch_algorithm(self.ref[ref_id].pocket.coor(), 
                                                        self.pocket[ref_id].coor())
        # aligned_ref_st_pos = np.dot(self.ref_st[ref_id].coordinates-centroid, rot.T) + centroid + trans
        aligned_ref = np.dot(self.ref[ref_id].coor() - centroid, rot.T) + centroid + trans
                             
        # 2. align the structure ligand to the reference ligand
        aligned_ref_ligand_pos = aligned_ref[self.ref[ref_id].ligand.idx,:]
        # ligand_pos = self.st.coordinates[self.ligand_idx_group[ref_id],:]
        (rot, trans, centroid, rmsd) = kabsch_algorithm(self.ligand[ref_id].coor(), 
                                                        aligned_ref_ligand_pos)
        # ligand = self.st.coordinates[self.ligand_indexes,:]
        # transform whole ligand
        transformed = np.dot(self.ligand['0'].coor() - centroid, rot.T) + centroid + trans
        return (transformed, centroid, rot, trans)


    def move_molecule_to(self, molecule_id:str, ref_id:str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get transformations to move protein coordinates to the supposed ligand.

        Args:
            molecule_id (str): id of molecule.
            ref_id (str): id of reference structure(ligand)

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: _description_
        """
        # 1. align the reference to the structure using the reference ligand
        # ref_ligand_pos = self.ref_st[ref_id].coordinates[self.ref_ligand_idx_group[ref_id],:]
        # ligand_pos = self.st.coordinates[self.ligand_idx_group[ref_id],:]
        (rot, trans, centroid, rmsd) = kabsch_algorithm(self.ref[ref_id].ligand.coor(), 
                                                        self.ligand[ref_id].coor())
        # aligned_ref_st_pos = np.dot(self.ref_st[ref_id].coordinates-centroid, rot.T) + centroid + trans
        aligned_ref = np.dot(self.ref[ref_id].coor() - centroid, rot.T) + centroid + trans
        
        # 2. align the structure pocket to the reference pocket
        aligned_ref_pocket_pos = aligned_ref[self.ref[ref_id].pocket.idx,:]
        # pocket_pos = self.st.coordinates[self.pocket_idx_group[ref_id],:]
        (rot, trans, centroid, rmsd) = kabsch_algorithm(self.pocket[ref_id].coor(), 
                                                        aligned_ref_pocket_pos) 
        transformed = np.dot(self.molecule[molecule_id].coor() -centroid, rot.T) + centroid + trans
        return (transformed, centroid, rot, trans)


    def save_molecule(self, molecule_id:str, path:str | Path, overwrite:bool=True) -> None:
        """Save molecular coordinates to an output file.

        Args:
            molecule_id (str): id of molecule.
            path (str | Path): output file path.
            overwrite (bool, optional): whether to overwrite. Defaults to True.
        """
        if isinstance(path, Path):
            path = path.as_posix()
        st = self.molecule[molecule_id].st
        st[self.molecule[molecule_id].idx].save(path, overwrite=overwrite)