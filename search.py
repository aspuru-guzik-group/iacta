import numpy as np
from rdkit_utils import *
import xtb_utils
import importlib
importlib.reload(xtb_utils)
    

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
opt = xtb.optimize("init_guess.mol", "optimized_guess.mol")
combined = Chem.MolFromMolFile(opt(), removeHs=False)
N = combined.GetNumAtoms()

# Get the indices of the bond to stretch
bond = get_bonds(combined, "CBr")[0]

# Start the run
stretch_factors = 1.1**np.arange(10)

import os
import shutil
os.makedirs("conformers", exist_ok=True)
os.makedirs("stretch", exist_ok=True)
MolToXYZFile(combined, "conformers/0.xyz")


conformers = ["0.xyz"]
stretch = stretch_factors[0]
xcontrol = {"constrain":
            ("force constant = 0.5",
             "distance: %i, %i, %f"% (bond[0],bond[1],
                                        stretch * bond[2]))}

# jobs = []
# for c in conformers:
#     jobs += [xtb.optimize("conformers/" + c, "stretch/" + c,
#                           xcontrol=xcontrol,
#                           log="stretch/" + c + ".log")]


# parameters from nanoreactor part of CREST paper
xcontrol["metadyn"] = ("save=100", "kpush=%f" % (0.02*N),
                       "alp=0.6")
xcontrol["md"] = ("shake=0","step=0.5")
xcontrol["wall"] = ("potential=logfermi",
                    "sphere: auto, all")


metarun = xtb_utils.xtb_run("xtb",
                            "conformers/0.xyz",
                            "--metadyn",
                            "--etemp 6000",
                            xcontrol=xcontrol)
metarun.start(False)

