import numpy as np
import cclib
import glob

def get_energies(files):
    out = []
    for p in files:
        v = []
        f = open(p, "r")
        for line in f:
            if line[1:5] == "ener":
                v += [float(line[8:25])]

        out += [v[:]]
    return np.array(out)

def get_smiles(file):
    f = open(file, "r")
    out = []
    for line in f:
        out += [line.rstrip()]
    return out

# Get products and reactant smiles
react_smi = get_smiles("mols/reactants.smi")
prod_smi = get_smiles("mols/products.smi")
energies = get_energies(glob.glob("mols/*.xyz"))

