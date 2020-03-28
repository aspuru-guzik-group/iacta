import numpy as np
import react
import xtb_utils
import os

# Initialize the xtb driver
xtb = xtb_utils.xtb_driver()
xtb.extra_args = ["-gfn2"]

# TODO REMOVE RDKIT
# Optimize the combined molecule
# opt = xtb.optimize("da.mol", "da.mol")
# combined = Chem.MolFromMolFile("da.mol", removeHs=False)
# MolToXYZFile(combined, "initial_guess.xyz")

# Get additional molecular parameters
N = 19
params = react.default_parameters(N)

# Constraints for the search
# -------------------------
bond_lengths = np.linspace(1.3, 1.75, 21)
constraints = [("force constant = 1.25",
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

