import numpy as np
from react import *
import xtb_utils
import os

# Initialize the xtb driver
xtb = xtb_utils.xtb_driver(scratch="/local-scratch/clavigne")
xtb.extra_args = ["-gfn2", "--etemp 1000"]


# Constraints for the search
# -------------------------
Natoms = 22
bond = (19, 22, 1.06990)
params = default_parameters(Natoms)
stretch_factors = np.linspace(1.0, 3.0, 50)
constraints = [("force constant = 0.5",
                "distance: %i, %i, %f"% (bond[0],bond[1],
                                         stretch * bond[2]))
               for stretch in stretch_factors]

# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
os.makedirs("output_gabes", exist_ok=True)
generate_starting_structures(xtb,
                             "catalyst_p_alkyne.xyz",
                             "output_gabes",
                             constraints,
                             params)

nmols = 0
for mtd_index in [0, 5, 10, 20]:
    # STEP 2: Metadynamics search
    # ---------------------------
    metadynamics_search(xtb,
                        mtd_index,
                        "output_gabes",
                        constraints,
                        params)


    # STEP 3: Reactions
    # -----------------
    trajs = react(xtb,
          mtd_index,
          "output_gabes",
          constraints,
          params)
    
    dump_trajectories(trajs, "output_gabes", offset=nmols)
    nmols += len(trajs)




