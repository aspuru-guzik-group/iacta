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

# combined = Chem.MolFromMolFile("optimized_guess.mol",
#                                removeHs=False)

# # Get the indices of the bond to stretch
# bond = get_bonds(combined, "CBr")[0]

# # Start the run
# stretch_factors = 1.1**np.arange(10)

# import os
# import shutil
# os.makedirs("conformers", exist_ok=True)
# os.makedirs("stretch", exist_ok=True)
# MolToXYZFile(combined, "conformers/0.xyz")
# conformers = ["0.xyz"]

# stretch = stretch_factors[0]
# xcontrol = open("xcontrol", "w")
# xcontrol.write("$constrain\n"+
#                "  force constant = 0.5\n" +
#                "  distance: %i, %i, %f\n"% (bond[0],bond[1],
#                                             stretch * bond[2]) +
#                "$end")
# xcontrol.close()

# for c in conformers:
#     xtb.optimize("conformers/" + c, "stretch/" + c,
#                  xcontrol="xcontrol",
#                  log="stretch/" + c + ".log")

# stretch = stretch_factors[0]
# xcontrol = open("xcontrol", "w")
# xcontrol.write("$constrain\n"+
#                "  force constant = 0.5\n" +
#                "  distance: %i, %i, %f\n"% (bond[0],bond[1],
#                                             stretch * bond[2]) +
#                "$end")
# xcontrol.write("$metadyn\n"+
#                "  save=10\n"+
#                "  kpush=1.0\n"+
#                "  alp=0.2\n"+
#                "$end")
# xcontrol.close()

# metarun = xtb_utils.xtb_run("xtb",
#                             "conformers/0.xyz",
#                             "--metadyn",
#                             xcontrol="xcontrol")

