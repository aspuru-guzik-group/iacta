import numpy as np
from rdkit_utils import *
import xtb_utils
import importlib
importlib.reload(xtb_utils)
    

XTB_PATH = "/h/182/clavigne/xtb/bin/xtb"
xtb = xtb_utils.xtb_driver(XTB_PATH)

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
combined = Chem.MolFromMolFile("optimized_guess.mol")

# Get the indices of the bond to stretch
bond = get_bonds(combined, "CBr")[0]

# Start the run
stretch_factors = 1.1**np.arange(10)

conformers = [[]]
stretch = stretch_factors[0]


