import numpy as np
import cclib
import glob

def get_energies(files):
    v = []
    for p in files:
        f = open(p, "r")
        for line in f:
            if line[1:5] == "ener":
                v += [float(line[8:25])]
    return v

def get_smiles(files):
    out = []
    for p in files:
        f = open(p, "r")
        for line in f:
            out += [line.split("\t")[0]]
    return out

# get energies
mtd_indices = [5]
energies = []
# energies = [get_energies(["reactions/reactants_%2.2i.xyz"% m
#                           for m in mtd_indices])]
k = 0
while True:
    try:
        new = get_energies(["reactions/prop%2.2i_%3.3i.xyz" % (m,k)
                            for m in mtd_indices])
    except FileNotFoundError:
        break

    energies += [new[:]]
    k+=1

energies += [get_energies(["reactions/products_%2.2i.xyz"% m
                          for m in mtd_indices])]    
energies = np.array(energies)

# Get products and reactant smiles
react_smi = get_smiles(["reactions/reactants_%2.2i.smi"% m
                        for m in mtd_indices])
prod_smi = get_smiles(["reactions/products_%2.2i.smi"% m
                          for m in mtd_indices])

reactions = {}
dummies = []
for i in range(len(react_smi)):
    reactants = react_smi[i]
    products = prod_smi[i]
    path = energies[:,i]

    if reactants == products:
        dummies += [path]
    else:
        if (reactants, products) in reactions:
            reactions[reactants,products] += [path]
        else:
            reactions[reactants,products] = [path]
dummies = np.array(dummies)
reaction = []
paths = []
for key, val in reactions.items():
    reaction += [key]
    paths += [np.array(val)]


