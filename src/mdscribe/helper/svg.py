from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from IPython.display import SVG


def view_svg(rdmol:Chem.Mol, highlight=[], width=300, height=300):
    """SVG depiction for the Jupyter Notebook.

    Reference:
        https://www.linkedin.com/pulse/using-rdkit-jupyter-notebooks-lee-davies/

    Args:
        rdmol (Chem.Mol): rdkit molecule.
        highlight (list, optional): highlighted atom indexes. Defaults to [].
        width (int, optional): width. Defaults to 300.
        height (int, optional): height. Defaults to 300.
    """
    rdmol2d = Chem.Mol(rdmol)
    AllChem.Compute2DCoords(rdmol2d)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(rdmol2d, highlightAtoms=highlight)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    display(SVG(svg.replace("svg:","")))