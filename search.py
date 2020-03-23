import numpy as np
from rdkit_utils import *
import vibrations
import xtb_utils
import importlib
importlib.reload(xtb_utils)
importlib.reload(vibrations)
    

xtb = xtb_utils.xtb_driver()

# Molecules
# smiles_react1="CNNC(C)(C)"
# smiles_react2="C1CCCC1([Br])"
smiles_react1="CI"
smiles_react2="CS"

# Build conformers
react1 = rdkit_generate_conformer(smiles_react1)
react2 = rdkit_generate_conformer(smiles_react2)
combined = mix2(react1, react2, off=5.0)
Chem.MolToMolFile(combined, "init_guess.mol")

# Optimize
opt = xtb.optimize("init_guess.mol", "optimized_guess.mol",
                   compute_hessian=True)
opt()
combined = Chem.MolFromMolFile("optimized_guess.mol",
                               removeHs=False)

# Load vibrational quantities
N = combined.GetNumAtoms()
hessian = vibrations.load_xtb_hessian("hessian_optimized_guess.mol")
M = np.array([at.GetMass() for at in combined.GetAtoms()])
atoms = np.array([at.GetSymbol() for at in combined.GetAtoms()])
W2,L = vibrations.normal_modes(M, hessian)
w = np.sqrt(abs(W2))
x0 = combined.GetConformer().GetPositions()
q0 = L.T.dot(vibrations.massweight(M, x0))


# Get the indices of the bond to stretch
bond = get_bonds(combined, "CI")[0]

# Start the run
stretch_factors = 1.1**np.arange(0,10)
tol_sparsify = 0.0

import os
import shutil
os.makedirs("conformers", exist_ok=True)
os.makedirs("stretch", exist_ok=True)
os.makedirs("meta", exist_ok=True)
os.makedirs("optim", exist_ok=True)

# parameters for conformer search, adapted from CREST
nstructures = 50
total_time = 0.5 * N
dumpstep = 1000 * total_time/nstructures        # 10 structures
xmetadyn = ("save=%i"%nstructures,
            "kpush=0.2", "alp=0.8")
xmd = ("shake=2", "step=5",
       "dump=%f"%dumpstep, "time=%f" % total_time)
xwall = ("potential=logfermi",
         "sphere: auto, all")

os.makedirs("conformers/s0", exist_ok=True)
MolToXYZFile(combined, "conformers/s0/00000.xyz")
conformers = [["00000.xyz"]]

# Step 1: stretch conformers
# for stretch_index,stretch in enumerate(stretch_factors):
stretch_index = 0; stretch = stretch_factors[0]

xconstrain = ("force constant = 0.5",
              "distance: %i, %i, %f"% (bond[0],bond[1],
                                       stretch * bond[2]))

print("getting new stretches")
current = conformers[-1]
new = []
ojobs = []
for c in current:
    ojobs += [xtb.optimize("conformers/s%i/"%stretch_index + c,
                          "stretch/" + c,
                           level="loose",
                           xcontrol=dict(
                               wall=xwall,
                               constrain=xconstrain))]
for k,j in enumerate(ojobs):
    print(k, len(ojobs))
    j.start()
    j.close()

print("doing metadyn")
mjobs = []
for c in current:
    mjobs += [xtb.metadyn("stretch/"+c,
                          "meta/"+c,
                          xcontrol=dict(
                              wall=xwall,
                              constrain=xconstrain,
                              metadyn=xmetadyn,
                              md=xmd))]
for k,j in enumerate(mjobs):
    print(k, len(mjobs))    
    j.start()
    j.close()


# Load all the geometries
import cclib
xyzs = []
for c in current:
    # from meta
    parser = cclib.ccopen("meta/"+c)
    data = parser.parse()
    xyzs += [x for x in data.atomcoords]

print("first-pass optimization of metad structures")
os.makedirs("conformers/s%i" % (stretch_index+1), exist_ok=True)
ojobs = []
for k,xyz in enumerate(xyzs):
    f = open("optim/%5.5i.xyz" %k, "w")
    f.write(geom2xyz(atoms, *xyz.T))
    f.close()

    # ojobs += [xtb.optimize(
    #     "optim/%5.5i.xyz" %k,
    #     "optim/%5.5i.xyz" %k,
    #     level="lax",
    #     xcontrol=dict(
    #         wall=xwall,
    #         constrain=xconstrain))]

# for k,job in enumerate(ojobs):
#     job()
#     print(k, len(ojobs))


# Sparsifying
print("doing sparsification of structures")
qs = []
for k in range(len(xyzs)):
    parser = cclib.ccopen("optim/%5.5i.xyz" %k)
    data = parser.parse()
    qs += [L.T.dot(vibrations.massweight(M,x)) - q0 for x in data.atomcoords]

# Convert to coherent state parameter alpha
nwf = len(qs)

# overlap between coherent states
# <a|b> = exp(-0.5 * (b^2 + a^2 - 2 b a))
S = np.zeros((nwf,nwf))
for i in range(nwf):
    for j in range(nwf):
        b = qs[i] * np.sqrt(w)
        a = qs[j] * np.sqrt(w)
        S[i,j] = np.exp(-0.5 * (b**2+a**2 - 2 *b * a).sum())

e, p = np.linalg.eigh(S)
mask = np.cumsum(e)/sum(e)>0.1
new_qs = p[:,mask].T.dot(qs)
new_xyzs = []
print("     sparsified from %i to %i" % (len(e), sum(mask)))
for q in new_qs:
    nx = L.dot(q + q0)              # invert normal mode transform
    nx = vibrations.inv_massweight(M, nx) # invert mass weighting
    new_xyzs += [np.reshape(nx, x0.shape)]

for k,xyz in enumerate(new_xyzs):
    f = open("conformers/s1/%5.5i.xyz" %k, "w")
    f.write(geom2xyz(atoms, *xyz.T))
    f.write("\n")
    f.close()
    
# # print("doing tight optimization of sparsified strauctures")
# # os.makedirs("conformers/s%i" % (stretch_index+1), exist_ok=True)
# # ojobs = []
# # for k,xyz in enumerate(new_xyzs):
# #     f = open("optim/%5.5i.xyz" %k, "w")
# #     f.write(geom2xyz(atoms, *xyz.T))
# #     f.close()

# #     ojobs += [xtb.optimize(
# #         "optim/%5.5i.xyz" %k,
# #         "conformers/s%i/%5.5i.xyz" % (stretch_index+1, k),
# #         level="tight",
# #         xcontrol=dict(
# #             wall=xwall,
# #             constrain=xconstrain))]

# # new = []
# # for k,job in enumerate(ojobs):
# #     print(k, len(ojobs))
# #     new += [os.path.basename(job())]

# conformers += [new]
