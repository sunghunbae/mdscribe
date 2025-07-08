from logging import raiseExceptions
import os
import sys
import argparse

from datetime import timedelta, datetime
from mdscribe.helper import PDBfile


def interpret_to_list(expr:str, enumerate_range:bool=False) -> list:
    """Interpret a expression of list of items.

    Args:
        expr (str): expression of single or multiple model(s)/chainId(s)/residues
        enumerate_range (bool, optional): enumerate a range.
            If set True, `3-5,9` will be interpreted and expanded to a list, `[3, 4, 5, 9]`.
            Defaults to False.
    Returns:
        list: 
    """
    if not enumerate_range:
        return [s.upper() for s in expr.split(",")]
    else:
        x = []
        for range_expr in expr.split(","):
            begin_end = list(map(int, range_expr.split("-")))
            if len(begin_end) == 1:
                x.extend(begin_end)
            elif len(begin_end) == 2:
                (begin, end) = begin_end
                x.extend(list(range(begin, end+1, 1)))
            else:
                raise ValueError(f"Invalid expression: {begin_end}")
        return x
    

def pdbrestore():    
    parser = argparse.ArgumentParser(description="PDB restore chainId/resSeq changed by pdb4amber",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r','--restore', dest="restore", type=str, 
                        help="pdb4amber output file xxx_renum.txt")
    parser.add_argument('-o','--outfile', dest="outfile", type=str, 
                        help="output filename")
    parser.add_argument('--overwrite', dest="overwrite", default=False, 
                        action="store_true", help="overwrite")
    parser.add_argument('pdb', help="PDB input file")

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if not os.path.exists(args.pdb):
        print(f"cannot open {args.pdb}")
        sys.exit(1)

    if args.outfile:
        if os.path.exists(args.outfile) and not args.overwrite:
            print(f"outfile exists. use --overwrite to overwrite")
            sys.exit(1)
    else:
        args.outfile = sys.stdout


    pdb = PDBfile(args.pdb)
    
    if args.restore and os.path.exists(args.restore):
        pdb.restore(args.restore)
        pdb.write(args.outfile)



def pdbextract():
    parser = argparse.ArgumentParser(description="PDB extract coordinates by chainId/resSeq",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m','--model', dest="model", type=str, help="model number(s) to extract")
    parser.add_argument('-c','--chain', dest="chain", type=str, help="chainId(s) to extract")
    parser.add_argument('-r','--resseq', dest="resseq", type=str, help="residue number(s) to extract")
    parser.add_argument('-n','--resname', dest="resname", type=str, help="residue name(s) to extract")
    parser.add_argument('-o','--outfile', dest="outfile", type=str, default="pdbextract.pdb", help="output filename")
    parser.add_argument('--overwrite', dest="overwrite", default=False, 
                        action="store_true", help="overwrite")
    parser.add_argument('pdb', help="PDB input file")

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)


    if not os.path.exists(args.pdb):
        print(f"cannot open {args.pdb}")
        sys.exit(1)

    if args.outfile:
        if os.path.exists(args.outfile) and not args.overwrite:
            print(f"outfile exists. use --overwrite to overwrite")
            sys.exit(1)
    else:
        args.outfile = sys.stdout


    if args.model :
        model = interpret_to_list(args.model, enumerate_range=True)
    else:
        model = None

    if args.chain :
        chain = interpret_to_list(args.chain)
    else:
        chain = None

    if args.resseq :
        resseq = interpret_to_list(args.resseq, enumerate_range=True)
    else:
        resseq = None

    if args.resname :
        resname = interpret_to_list(args.resname)
    else:
        resname = None

    pdb = PDBfile(args.pdb)
    pdb.write(args.outfile, model, chain, resname, resseq)