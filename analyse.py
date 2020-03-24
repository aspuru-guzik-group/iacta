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

# get energies
mtd_indices = [0,3,8,16]
energies = [get_energies(["reactions/reactants_%2.2i.xyz"% m
                          for m in mtd_indices])]
k = 1
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

# get products and reactants
# for i in mtd_indices:
i= 0
products = cclib.ccopen("reactions/reactants_%2.2i.xyz" % i)
data = products.parse()



