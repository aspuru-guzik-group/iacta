import numpy as np
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms, rdmolops
import os
from rdkit_utils import fractional2bondtype
from xtb_utils import read_xtb_output



atoms,charges,positions,wbo = read_xtb_output("test_xtb2smiles/xtbopt.xyz")
em = Chem.EditableMol(Chem.Mol())

for atom, charge in zip(atoms, charges):
    new_atom = Chem.Atom(atom)
    new_atom.SetFormalCharge(int(charge))
    em.AddAtom(new_atom)

for i,j,fract_order in wbo:
    print(fract_order, fractional2bondtype(fract_order), atoms[i], atoms[j])
    em.AddBond(i,j,order=fractional2bondtype(fract_order))
    
m = em.GetMol()
m = Chem.RemoveHs(m)
Chem.SanitizeMol(m)
print(Chem.MolToSmiles(m))
