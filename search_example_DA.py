import numpy as np
from rdkit_utils import *
import react
import xtb_utils
import os

# Initialize the xtb driver
xtb = xtb_utils.xtb_driver()
xtb.extra_args = ["-gfn2"]

# # STEP 0: INITIAL CONFORMER BUILDING / SETTING PARAMETERS
# # ---------------------------------------------------------------------------- 
# # Molecules
# smiles_react1="C=C-C=O"
# smiles_react2="C1=CCC=C1"

# # Build conformers
# react1 = rdkit_generate_conformer(smiles_react1)
# react2 = rdkit_generate_conformer(smiles_react2)
# combined = mix2(react1, react2, off=5.0)
# Chem.MolToMolFile(combined, "da.mol")

# # Optimize the combined molecule
# opt = xtb.optimize("da.mol", "da.mol")
combined = Chem.MolFromMolFile("da.mol", removeHs=False)
MolToXYZFile(combined, "initial_guess.xyz")

# Get additional molecular parameters
N = 19
params = react.default_parameters(N)

# Constraints for the search
# -------------------------
bond_lengths = np.linspace(1.3, 1.65, 21)
constraints = [("force constant = 0.5",
                "distance: %i, %i, %f"% (1,2, length))
               for length in bond_lengths]
mtd_indices = [0, 5, 10, 15, 20]


# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
os.makedirs("output", exist_ok=True)
react.generate_initial_structures(
    xtb, "output",
    "initial_guess.xyz",
    constraints,
    params)

# STEP 2: Metadynamics
# ----------------------------------------------------------------------------
react.metadynamics_search(
    xtb, "output",
    mtd_indices,
    constraints,
    params,
    nthreads=2)

# STEP 2: Reactions
# ----------------------------------------------------------------------------
react.react(
    xtb, "output",
    mtd_indices,
    constraints,
    params,
    nthreads=4)

