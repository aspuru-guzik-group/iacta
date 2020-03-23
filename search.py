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
not_rotvib = abs(W2) > 1e-10
w = np.sqrt(W2[not_rotvib])
L = L[:,not_rotvib]
x0 = combined.GetConformer().GetPositions()
q0 = L.T.dot(vibrations.massweight(M, x0))


# Get the indices of the bond to stretch
bond = get_bonds(combined, "CI")[0]

# Start the run
stretch_factors = 1.2**np.arange(10)

import os
import shutil
os.makedirs("conformers", exist_ok=True)
os.makedirs("stretch", exist_ok=True)
os.makedirs("meta", exist_ok=True)
os.makedirs("optim", exist_ok=True)



# parameters from nanoreactor part of CREST paper
# xcontrol["metadyn"] = ("save=100", "kpush=%f" % (0.02*N),
#                        "alp=0.6")
# xcontrol["md"] = ("shake=0", "step=1",
#                   "time=100")    # you probably want to change that to 10

# parameters for conformer search, adapted from CREST
total_time = 0.5 * N
dumpstep = 1000 * total_time/10        # 10 structures
xmetadyn = ("save=10","kpush=0.2", "alp=0.8")
xmd = ("shake=2", "step=5",
                  "dump=%f"%dumpstep, "time=%f" % total_time)
xwall = ("potential=logfermi", "sphere: auto, all")

os.makedirs("conformers/s0", exist_ok=True)
MolToXYZFile(combined, "conformers/s0/00000.xyz")
conformers = [["00000.xyz"]]

# Step 1: stretch conformers
for stretch_index,stretch in enumerate(stretch_factors):
    xconstrain = ("force constant = 0.5",
                  "distance: %i, %i, %f"% (bond[0],bond[1],
                                           stretch * bond[2]))
    
    print("getting new stretch")
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

        # from stretch
        parser = cclib.ccopen("stretch/"+c)
        data = parser.parse()
        xyzs += [x for x in data.atomcoords]

    print("doing optimization of metad structures")
    os.makedirs("conformers/s%i" % (stretch_index+1), exist_ok=True)
    ojobs = []
    for k,xyz in enumerate(xyzs):
        f = open("optim/%5.5i.xyz" %k, "w")
        f.write(geom2xyz(atoms, *xyz.T))
        f.close()

        ojobs += [xtb.optimize(
            "optim/%5.5i.xyz" %k,
            # "optim/%5.5i.xyz" %k,
            "conformers/s%i/%5.5i.xyz" % (stretch_index+1, k),
            level="normal",
            xcontrol=dict(
                wall=xwall,
                constrain=xconstrain))]


    new = []
    for k,job in enumerate(ojobs):
        print(k, len(ojobs))
        new += [os.path.basename(job())]

    conformers += [new]


# qs = []
# for k in range(len(xyzs)):
#     parser = cclib.ccopen("optim/%5.5i.xyz" %k)
#     data = parser.parse()
#     qs += [L.T.dot(vibrations.massweight(M,x)) - q0 for x in data.atomcoords]

# # Convert to coherent state parameter alpha
# alphas = np.array([q*np.sqrt(w) for q in qs]) 
# nwf = len(alphas)

# # overlap between coherent states
# # <a|b> = exp(-0.5 * (b^2 + a^2 - 2 b a))
# S = np.zeros((nwf,nwf))
# for i in range(nwf):
#     for j in range(nwf):
#         b = alphas[i]
#         a = alphas[j]
#         S[i,j] = np.exp(-0.5 * (b.dot(b) + a.dot(a) - 2 * b.dot(a)))
