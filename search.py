import numpy as np
from rdkit_utils import *
import vibrations
import xtb_utils
import importlib
importlib.reload(xtb_utils)
importlib.reload(vibrations)
    

xtb = xtb_utils.xtb_driver()

# Molecules
smiles_react1="CNNC(C)(C)"
smiles_react2="C1CCCC1([Br])"


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
N = combined.GetNumAtoms()
hessian = vibrations.load_xtb_hessian("hessian_optimized_guess.mol")

"""
# Get the indices of the bond to stretch
bond = get_bonds(combined, "CBr")[0]

# Start the run
stretch_factors = 1.1**np.arange(10)

import os
import shutil
os.makedirs("conformers", exist_ok=True)
os.makedirs("stretch", exist_ok=True)
os.makedirs("meta", exist_ok=True)
MolToXYZFile(combined, "conformers/0.xyz")


xcontrol = {}
# parameters from nanoreactor part of CREST paper
# xcontrol["metadyn"] = ("save=100", "kpush=%f" % (0.02*N),
#                        "alp=0.6")
# xcontrol["md"] = ("shake=0", "step=1",
#                   "time=100")    # you probably want to change that to 10

# parameters for conformer search, adapted from CREST
xcontrol["metadyn"] = ("save=10","kpush=0.1", "alp=0.5")
xcontrol["md"] = ("shake=2", "step=5", "time=%i"%(N//2))
xcontrol["wall"] = ("potential=logfermi", "sphere: auto, all")

conformers = [["0.xyz"]]

# Step 1: stretch conformers
stretch = stretch_factors[0] *1.3
xcontrol["constrain"] = ("force constant = 0.5",
                         "distance: %i, %i, %f"% (bond[0],bond[1],
                                                  stretch * bond[2]))

current = conformers[-1]
new = []
ojobs = []
for c in current:
    ojobs += [xtb.optimize("conformers/" + c,
                          "stretch/" + c,
                           level="loose",
                           xcontrol=xcontrol)]
for j in ojobs:
    j.start()
    j.close()

mjobs = []
for c in current:
    mjobs += [xtb.metadyn("stretch/"+c,
                          "meta/"+c,
                          xcontrol=xcontrol)]
for j in mjobs:
    j.start()
    j.close()
    
# Load all the geometries
import cclib
coords = []
for c in current:
    parser = cclib.ccopen("meta/"+c)
    data = parser.parse()
    coords += [data.atomcoords]

# coords is an array of Nconf x Natom x 3 of all the coordinates generated in
# the metadynamics. Now we want to transform coords to internal coordinates,
# filter it, optimize all the remaining geometries and filter again.
coords = np.concatenate(coords)


"""
