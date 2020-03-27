import numpy as np
from rdkit_utils import *
import react
import xtb_utils
import os

import importlib
importlib.reload(xtb_utils)
importlib.reload(react)

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
params = react.default_parameters(N, nmtd=10)

# Constraints for the search
# -------------------------
stretch_factors = np.linspace(1.0, 3.0, 40)
constraints = [("force constant = 0.5",
                "distance: %i, %i, %f"% (bond[0],bond[1],
                                         stretch * bond[2]))
               for stretch in stretch_factors]
mtd_indices = [0, 5, 10, 15]



# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
out = "output"
os.makedirs(out, exist_ok=True)
MolToXYZFile(combined, out + "/initial_geom.xyz")

nthreads = 2
from concurrent.futures import ThreadPoolExecutor
verbose = True

if verbose:
    print("Performing initial stretching...")
    
structures, energies, opt_indices = react.successive_optimization(
    xtb, "init_guess.xyz",
    constraints, params)
react.dump_succ_opt(out + "/init", structures,energies,opt_indices)

with ThreadPoolExecutor(max_workers=nthreads) as pool:
    for mtd_index in mtd_indices:
        pool.submit(
            react.metadynamics_job(
                xtb, mtd_index,
                out+"/init", out+"/metadyn",
                constraints, params))


# load all the structures
meta = out+"/metadyn"
freact = out+"/reacts"
os.makedirs(freact, exist_ok=True)


nreact = 0
with ThreadPoolExecutor(max_workers=nthreads) as pool:
    for mtd_index in mtd_indices:
        structures, energies = react.read_trajectory(
            meta + "/mtd%4.4i.xyz" % mtd_index)
        
        if verbose:
            print("starting MTD job %i..." % mtd_index)


        for s in structures:
            rjob = react.reaction_job(xtb,
                                      s,
                                      mtd_index,
                                      freact + "/react%5.5i/" % nreact,
                                      constraints,
                                      params,
                                      verbose=verbose)
            pool.submit(rjob)
            nreact = nreact + 1
                           


                      



