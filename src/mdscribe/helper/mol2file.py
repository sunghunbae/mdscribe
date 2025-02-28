import sys
import re
import pathlib
import numpy as np


class MOL2file:
    """Class for handling MOL2 file."""
    
    # class bound variable
    name_syntax = re.compile(r'(?P<element>[a-zA-Z]+)(?P<serial>\d+)')

    def __init__(self, filename:str | pathlib.Path):
        self.n_atoms = 0
        self.blocks = {}
        self.coords = None # To be modified
        self.pdblines = [] # reserved for PDB conversion

        cursor = None
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                if "@<TRIPOS>" in line:
                    cursor = line.split("@<TRIPOS>")[1].strip().lower()
                    self.blocks[cursor] = []
                    continue
                elif line.startswith('#') or line == '\n':
                    continue
                self.blocks[cursor].append(line)

        try:
            self.n_atoms = len(self.blocks['atom'])
            assert self.n_atoms > 0
        except:
            print("MOL2 block has no atom")
            sys.exit(0)

        try:
            assert len(self.blocks['bond']) > 0
        except:
            print("MOL2 block has no bond")
            sys.exit(0)

        self.coords = np.zeros((self.n_atoms, 3), dtype=np.float32)

        altLoc = ' '
        chainId = ' '
        iCode = ' '
        occupancy = 1.0
        tempFactor = 0.0
        segId = ' '
        for i, a in enumerate(self.blocks['atom']):
            aid, name, x, y, z, atom_type, resid, resname, charge = a.split()[:9]
            x, y, z = float(x), float(y), float(z)
            self.coords[i, :] = x, y, z
            m = MOL2file.name_syntax.match(name)
            element = m.group("element")
            self.pdblines.append(
                f"HETATM {(i+1):>5} {name:<4}{altLoc}{resname:<3} {chainId}{resid:<4}{iCode}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempFactor:6.2f}      {segId:<4}{element:<2}{charge:<2}"
            )
            
        for i, b in enumerate(self.blocks['bond']):
            bid, atom_1, atom_2, bond_type = b.split()[:4]
            _pdb  = f"CONECT {atom_1:>5} {atom_2:>5}"
            self.pdblines.append(_pdb)
        

    def write(self, filename:str | pathlib.Path) -> None:
        """Write to a file.

        Args:
            filename (str | pathlib.Path): output filename.
        """
        with open(filename, 'w') as f:
            i = 0
            for cursor in ['molecule', 'atom', 'bond', 'substructure']:
                if not cursor in self.blocks:
                    continue
                f.write(f"@<TRIPOS>{cursor.upper()}\n")
                i += 1
                if cursor == 'atom':
                    k = 0
                for line in self.blocks[cursor]:
                    if cursor != 'atom':
                        f.write(line)
                    else:    
                        aid, name, x, y, z, atom_type, resid, resname, charge = line.split()[:9]
                        # override x,y,z with self.coords
                        f.write("{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}\n".format(
                            aid, 
                            name, 
                            self.coords[k, 0], 
                            self.coords[k, 1], 
                            self.coords[k, 2],
                            atom_type,
                            resid,
                            resname,
                            float(charge),
                        ))
                        k += 1


    def transform(self, R:np.ndarray, T:np.ndarray) -> None:
        """Apply rotation/translation to the coordinates.

        Args:
            R (np.ndarray): 3x3 rotation matrix
            T (np.ndarray): translation vector
        """
        com = np.mean(self.coords, axis=0)
        self.coords = np.dot(self.coords - com, R) + T