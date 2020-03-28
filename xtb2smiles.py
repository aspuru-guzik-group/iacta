#!/usr/bin/env python
import numpy as np
from rdkit import Chem
from rdkit_utils import fractional2bondtype
from xtb_utils import read_xtb_output

def xtb2smiles(xyzfile, sanitize=False):
    """Transform the result of an xtb run into a smiles.

    Parameters:
    -----------
    xyzfile (str) : path to the xyz file output by an xtb calculation. The
    directory containing xyzfile is assumed to contain wbo and charges xtb
    output files.

    Returns:
    --------
    smiles (str) : The molecular SMILES, as a string.
    """
    atoms,charges,positions,wbo = read_xtb_output(xyzfile)
    em = Chem.EditableMol(Chem.Mol())

    for atom, charge in zip(atoms, charges):
        new_atom = Chem.Atom(atom)
        new_atom.SetFormalCharge(int(round(charge)))
        em.AddAtom(new_atom)

    for i,j,fract_order in wbo:
        em.AddBond(i,j,order=fractional2bondtype(fract_order))

    m = em.GetMol()
    if sanitize:
        m = Chem.RemoveHs(m, sanitize=1)
        Chem.SanitizeMol(m)
    else:
        m = Chem.RemoveHs(m, sanitize=0)
    return Chem.MolToSmiles(m)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Transform the output of an xtb calculation into a SMILES string.")
    parser.add_argument("xyzfile", help="Path to the xyz file output by an xtb"+
    "calculation. The directory containing xyzfile is assumed to contain wbo"+
    "and charges xtb output files.")
    parser.add_argument("-s",
                        help="Flag to sanitize the RdKit molecule. Defaults to"+
                        "false.",action="store_true")
    
    args = parser.parse_args()
    print(xtb2smiles(args.xyzfile, sanitize=args.s))
