import numpy as np
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms, rdmolops
import copy
import subprocess
import os

XTB_PATH = "/h/182/clavigne/xtb/bin/xtb"
def xtb_optimize(mol_file):
    dir = os.path.dirname(mol_file)
    fn = os.path.basename(mol_file)
    xtbout = open(dir+"/xtbout", "w")
    xtberr = open(dir+"/xtb.err", "w")
    subprocess.run([XTB_PATH, fn, "--opt"],
                           stderr=xtberr,
                           stdout=xtbout,
                           cwd=dir)
    xtbout.close()
    xtberr.close()
    return dir + "/xtbopt.mol"
    

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

def combine(*molecules):
    molecules2 = [mol if type(mol)==dict
                  else dict(mol=mol,constraints=[])
                  for mol in molecules]
    
    index_offset = 1
    inds = []
    for m in molecules2:
        inds += [index_offset]
        index_offset += m['mol'].GetNumAtoms()
        
    mols = mix2(molecules2[0]['mol'], molecules2[1]['mol'])
    for m in molecules2[2:]:
        mols = mix2(mols, m['mol'])

    constraint_string = "$constrain\n  force constant = 0.5\n"
    for m, off in zip(molecules2, inds):
        for at_i, at_j, l in m['constraints']:
            constraint_string += "  distance: %i, %i, %f\n" %(
                at_i + off,
                at_j + off, l)
    constraint_string += "$end"

    return mols, constraint_string
    
    
##########################################################
# REACTION: Electrophilic aromatic substitution on aniline

# Molecules
mol_smiles = dict(aniline="c1ccccc1(N)",
                 nitronium="[Al](Cl)(Cl)(Cl).ClCl")

# Optimize using xtb
mol_xtbopt = {}
os.makedirs("mols", exist_ok = True)
for name,smiles in mol_smiles.items():
    os.makedirs("mols/" + name, exist_ok = True)
    mol = rdkit_generate_conformer(smiles)
    fn = "mols/" + name + "/init.mol"
    Chem.MolToMolFile(mol,fn)
    
    # Optimize molecule with xtb
    mol_xtbopt[name] = xtb_optimize(fn)

# Set up reactions on aniline
aniline = Chem.MolFromMolFile(mol_xtbopt['aniline'],
                              removeHs=False)
aniline_react = react(aniline, "C[H]")

# Load nitronium
nitronium = Chem.MolFromMolFile(mol_xtbopt['nitronium'],
                                removeHs=False)

# Set up ts searches
os.makedirs("ts/", exist_ok=True)
for reaction_id,aniline_reagent in enumerate(aniline_react):
    folder = "ts/react%i/"% reaction_id
    os.makedirs(folder, exist_ok=True)
    total,constr = combine(aniline_reagent, nitronium)
    Chem.MolToMolFile(total, folder+"init.mol")
    cfile = open(folder + ".constrains", "w")
    cfile.write(constr)
    cfile.close()

    


"""
# Load molecules
ethyl_cl = Chem.MolFromMolFile("ethyl_cl/xtbopt.mol",
                               removeHs=False)
fecl = Chem.MolFromMolFile("fecl/xtbopt.mol",
                           removeHs=False)
benzene = Chem.MolFromMolFile("benzene/xtbopt.mol",
                              removeHs=False)

react(get_bond(ethyl_cl, "CCl"), constraints)


total, const = combine((ethyl_cl, fecl, benzene), constraints)


# make the reaction directory
os.makedirs("react1", exist_ok = True)
Chem.MolToMolFile(total, "react1/init.mol")
cfile = open("react1/.constrains", "w")
cfile.write(const)
cfile.close()
"""        



