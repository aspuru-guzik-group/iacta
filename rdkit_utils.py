from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms, rdmolops
import numpy as np
import copy

def rdkit_generate_conformer(SMILES,
                              selenium_bugfix=True):
    """Generate conformer from SMILES using rdkit.
    
    TODO:doc
    """
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    charge = Chem.GetFormalCharge(mol)
    atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    if selenium_bugfix:
        # Convert Se -> S to fix a stupid bug in rdkit.
        for i,symbol in enumerate(atom_symbols):
            if symbol == 'Se':
                mol.GetAtomWithIdx(i).SetAtomicNum(16)

    cid = AllChem.EmbedMolecule(mol)
    res = AllChem.MMFFOptimizeMolecule(mol)
    return mol

def get_bonds(mol, smarts):
    # Match groups
    subs = mol.GetSubstructMatches(
        # CC - - - Cl
        Chem.MolFromSmarts(smarts))
    
    # Get bond
    bonds = []
    for i,j in subs:
        bonds += [mol.GetBondBetweenAtoms(i,j)]
    return bonds

def react(mol, smarts, stretch_factor = 1.2):
    bonds = get_bonds(mol, smarts)
    
    # Elongate bond
    out = []
    for bond in bonds:
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        length = rdMolTransforms.GetBondLength(
            mol.GetConformer(), begin, end)

        nmol = copy.deepcopy(mol)
        rdMolTransforms.SetBondLength(nmol.GetConformer(),
                                      begin, end,
                                      length * stretch_factor)
        out += [dict(mol=nmol,
                     constraints=[(begin, end, length * stretch_factor)])]
    return out


def mix2(mol1, mol2, off=10.0):
    dir=np.array([0,0,1])
    pos1 = mol1.GetConformer().GetPositions()
    pos2 = mol2.GetConformer().GetPositions() + 100*dir/np.sum(dir)
    
    # Distance matrix between 1 and 2
    dmat12 = np.sqrt(np.array([[np.sum(ri - rj)**2
                                for ri in pos2] for rj in pos1]))
    
    # find the two atoms that are the closest
    a1, a2 = np.unravel_index(np.argmin(dmat12), dmat12.shape)
    # vector for offset
    vec = pos2[a2] - pos1[a1]
    vec = vec/np.sum(vec) * off

    # offset is 10 * vec
    offset = Geometry.Point3D(vec[0], vec[1], vec[2])

    # Combine the molecules
    combined = rdmolops.CombineMols(mol1,mol2, offset)
    return combined


