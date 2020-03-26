import numpy as np
from rdkit_utils import *
from react import *
import xtb_utils
import os

# Initialize the xtb driver
xtb = xtb_utils.xtb_driver()
xtb.extra_args = ["-gfn2"]

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
params = default_parameters(N)

# Constraints for the search
# -------------------------
stretch_factors = np.linspace(1.0, 3.0, 30)
constraints = [("force constant = 0.5",
                "distance: %i, %i, %f"% (bond[0],bond[1],
                                         stretch * bond[2]))
               for stretch in stretch_factors]



# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
os.makedirs("output", exist_ok=True)
MolToXYZFile(combined, "init_guess.xyz")
generate_starting_structures(xtb,
                             "init_guess.xyz",
                             "output",
                             constraints,
                             params)


nmols = 0
for mtd_index in [0,5,10]:
    # STEP 2: Metadynamics search
    # ---------------------------
    metadynamics_search(xtb,
                        mtd_index,
                        "output",
                        constraints,
                        params)


    # STEP 3: Reactions
    # -----------------
    trajs = react(xtb,
                  mtd_index,
                  "output",
                  constraints,
                  params)

    dump_trajectories(trajs, "output", offset=nmols)
    nmols += len(trajs)



