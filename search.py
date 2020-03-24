import numpy as np
from rdkit_utils import *
import xtb_utils
import importlib
import os
import shutil
import subprocess
import cclib
importlib.reload(xtb_utils)

# Initialize the xtb driver
xtb = xtb_utils.xtb_driver()
xtb.extra_args = ["-P 1", "-gfn1"]

# STEP 0: INITIAL CONFORMER BUILDING / SETTING PARAMETERS
# ---------------------------------------------------------------------------- 

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

# Optimize the combined molecule
opt = xtb.optimize("init_guess.mol", "optimized_guess.mol")
combined = Chem.MolFromMolFile(opt(), removeHs=False)

# Get the indices of the bond to stretch
bond = get_bonds(combined, "CI")[0]

# Get additional molecular parameters
N = combined.GetNumAtoms()
atoms = np.array([at.GetSymbol() for at in combined.GetAtoms()])

# Parameters for the search
# -------------------------
# bond stretch factors all
stretch_factors = np.linspace(1.0, 3.0, 20)
# indices to start mtds at
mtd_indices = [0,3,8,16]

nconstraints = len(stretch_factors)
def constraint(i):
    xconstrain = ("force constant = 0.5",
                  "distance: %i, %i, %f"% (bond[0],bond[1],
                                           stretch_factors[i] * bond[2]))
    return xconstrain

# number of mtd structures to generate at each stretch factor
mtd_nstructures = 20

# xTB additional parameters
xwall = ("potential=logfermi",
         "sphere: auto, all")

# Metadynamics parameters (somewhat adapted from CREST)
total_time = 0.2 * N            # TODO: set to 1.0
dumpstep = 1000 * total_time/mtd_nstructures
xmetadyn = ("save=10", 
            "kpush=0.2",
            "alp=0.8")
xmd = ("shake=0",
       "step=1",
       "dump=%f"%dumpstep,
       "time=%f" % total_time)




# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
os.makedirs("metadyn", exist_ok=True)
MolToXYZFile(combined, "metadyn/current.xyz")

print("Stretching initial structure to generate guesses for metadyn...")
print("---------------------------------------------------------------")

for i in range(nconstraints):
    print("    %i out of %i" %(i+1, nconstraints))
    opt = xtb.optimize("metadyn/current.xyz",
                       "metadyn/current.xyz",
                       level="tight",
                       xcontrol=dict(
                           wall=xwall,
                           constrain=constraint(i)))
    opt()
    # add to metadyn starting structures
    shutil.copy("metadyn/current.xyz", "metadyn/in%5.5i.xyz" % i)


# STEP 2: Metadynamics
# ----------------------------------------------------------------------------
print("")
print("Performing metadyn...")
print("---------------------")

# TODO This part can be easily parallelized
metajobs = []
for i in mtd_indices:
    print("    MTD - %i" %(i+1))
    
    metajobs += [xtb.metadyn("metadyn/in%5.5i.xyz" % i,
                             "metadyn/out%5.5i.xyz" % i,
                             xcontrol=dict(
                                 wall=xwall,
                                 constrain=constraint(i),
                                 metadyn=xmetadyn,
                                 md=xmd))]
    metajobs[-1]()

# STEP 3: Take the metadyn results and optimize them through the stretches
# ----------------------------------------------------------------------------
print("")
print("Performing reactions...")
print("---------------------")
os.makedirs("reactions", exist_ok=True)

# for i in mtd_indices:
i = 0
print("MTD states: %i" %(i+1))

shutil.copy("metadyn/out%5.5i.xyz" % i,
            "reactions/current.xyz")

for j in range(i, nconstraints):
    print("     ----> forward  %i out of %i" %(j+1, nconstraints))

    opt = xtb.multi_optimize("reactions/current.xyz",
                             "reactions/current.xyz",
                             level="normal",
                             xcontrol=dict(
                                 wall=xwall,
                                 constrain=constraint(j)))
    opt()
    shutil.copyfile("reactions/current.xyz",
                    "reactions/prop%2.2i_%3.3i.xyz" % (i,j))

print("     ----> forward to products")
opt = xtb.multi_optimize("reactions/current.xyz",
                         "reactions/products_%2.2i.xyz" %i,
                         level="normal",
                         xcontrol=dict(wall=xwall))
opt()

# copy starting point for backward reaction dynamics
shutil.copyfile("reactions/prop%2.2i_%3.3i.xyz" % (i,i),
                "reactions/current.xyz")


for j in range(i-1, 0, -1):
    print("     ----> backward %i out of %i" %(j+1, nconstraints))

    opt = xtb.multi_optimize("reactions/current.xyz",
                             "reactions/current.xyz",
                             level="normal",
                             xcontrol=dict(
                                 wall=xwall,
                                 constrain=constraint(j)))
    opt()
    shutil.copyfile("reactions/current.xyz",
                    "reactions/prop%2.2i_%3.3i.xyz" % (i,j))        


print("     ----> backward to reactants")
opt = xtb.multi_optimize("reactions/current.xyz",
                         "reactions/reactants_%2.2i.xyz" % i,
                         level="normal",
                         xcontrol=dict(wall=xwall))
opt()


