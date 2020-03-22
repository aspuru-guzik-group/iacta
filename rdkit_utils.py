from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms, rdmolops
import numpy as np
import copy

def geom2xyz(atoms, x, y, z, comment=""):
    """Transform a molecular geometry to a .xyz format string."""
    geomstring = "%i\n%s\n"%(len(atoms), comment)
    for i in range(len(atoms)):
        geomstring += " {atom}\t{x}\t{y}\t{z}\n".format(atom=atoms[i],
                                                        x=x[i],
                                                        y=y[i],
                                                        z=z[i])
    geomstring = geomstring[:-1]
    return geomstring

def MolToXYZFile(mol, file, comment=""):
    atoms = [at.GetSymbol() for at in mol.GetAtoms()]
    x,y,z = mol.GetConformer().GetPositions().T
    f = open(file, "w")
    f.write(geom2xyz(atoms, x, y ,z, comment))
    f.close()
    
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
        bonds += [(i+1,j+1,
                   rdMolTransforms.GetBondLength(mol.GetConformer(), i, j))]
        
    return bonds

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

