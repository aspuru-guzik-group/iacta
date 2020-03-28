import numpy as np
import react
import xtb_utils
import os
import pybel

# Initialize the xtb driver
xtb = xtb_utils.xtb_driver()
xtb.extra_args = ["-gfn2"]

# Get additional molecular parameters
mol = next(pybel.readfile("xyz", "./init_DielsAlder.xyz"))
N = len(mol.atoms)
bond = (1,2)
eq_length = eq_length = mol.OBMol.GetBond(*bond).GetEquibLength()
params = react.default_parameters(N)

# Constraints for the search
# -------------------------
stretch_factors = np.linspace(0.95, 1.4, 21)
constraints = [("force constant = 1.25",
                "distance: %i, %i, %f"% (*bond, eq_length * stretch))
               for stretch in stretch_factors]
mtd_indices = [0, 5, 10, 15, 20]


# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
os.makedirs("output", exist_ok=True)
react.generate_initial_structures(
    xtb, "output",
    "init_DielsAlder.xyz",
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

