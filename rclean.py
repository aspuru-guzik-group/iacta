import numpy as np
import react
import xtb_utils
import os
import shutil
import argparse
from ReactionPathway import *

parser = argparse.ArgumentParser(
    description="Driver for the refinement routine.",
    )
parser.add_argument("folder",
                    help="Folder containing the results of the search.",
                    type=str)
parser.add_argument("-b",
                    help="Minimum barrier size between species (in kcal/mol)."
                    +" Defaults to 2 kcal/mol.", default=2.0,
                    type=float)
parser.add_argument("-T",
                    help="Number of threads to use.",
                    type=int, default=1)
"""
also add no opt ts
parser.add_argument("--no-opt",
                    help="Start with an xtb optimization (defaults to true).",
                    action="store_true")
"""
parser.add_argument("--gfn",
                    help="gfn version. Defaults to GFN 2", default="2",
                    type=str)
parser.add_argument("--etemp",
                    help="Electronic temperature. Defaults to 300 K",
                    default="300.0",
                    type=str)
parser.add_argument("--log-level",
                    help="Level of debug printout (see react.py for details).",
                    default=0, type=int)


if "LOCALSCRATCH" in os.environ:
    scratch = os.environ["LOCALSCRATCH"]
else:
    print("warning: $LOCALSCRATCH not set")
    scratch = "."

hartree2ev = 27.2113860217
ev2kcal = 23.061
kcal2hartree = 1/(hartree2ev*ev2kcal)
    
args = parser.parse_args()

# Prepare output files
# --------------------
workdir = args.folder
# TODO: make an overwite? argument
os.makedirs(workdir + "/reactions", exist_ok=True)

# Read all pathways and species
pathways, species = read_all_pathways(workdir, args.b * kcal2hartree)

# Dump:

# 


"""
# Initialize the xtb driver
# -------------------------
xtb = xtb_utils.xtb_driver(scratch=scratch)
xtb.extra_args = ["--gfn " + args.gfn, "--etemp " + args.etemp]

if not args.no_opt:
    xtb.optimize(init, init, level="vtight")

# Get additional molecular parameters
# -----------------------------------
atoms, positions = xtb_utils.read_xyz(init)
N = len(atoms)
bond_length0 = np.sqrt(np.sum((positions[args.atoms[0]-1] -
                               positions[args.atoms[1]-1])**2))
bond = (args.atoms[0], args.atoms[1], bond_length0)

# Initialize parameters
# ---------------------
params = react.default_parameters(N,
                                  shake=args.shake_level,
                                  nmtd=args.mtdn,
                                  log_level=args.log_level)

# Constraints for the search
# -------------------------
stretch_factors = np.linspace(args.s[0], args.s[1], args.sn)
print("Stretching bond between atoms %s%i and %s%i"
      %(atoms[bond[0]-1],bond[0], atoms[bond[1]-1],bond[1]))
print("    with force constant %f" % args.force)
print("    between %7.2f and %7.2f A (%4.2f to %4.2f x bond length)"
      % (min(stretch_factors)*bond[2], max(stretch_factors)*bond[2],
         min(stretch_factors), max(stretch_factors)))
print("    discretized with %i points" % len(stretch_factors))
constraints = [("force constant = %f" % args.force,
                "distance: %i, %i, %f"% (bond[0],bond[1],
                                         stretch * bond[2]))
               for stretch in stretch_factors]
mtd_indices = args.mtdi


# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
react.generate_initial_structures(
    xtb, out_dir,
    init,
    constraints,
    params)

# STEP 2: Metadynamics
# ----------------------------------------------------------------------------
react.metadynamics_search(
    xtb, out_dir,
    mtd_indices,
    constraints,
    params,
    nthreads=args.T)


# STEP 2: Reactions
# ----------------------------------------------------------------------------
react.react(
    xtb, out_dir,
    mtd_indices,
    constraints,
    params,
    nthreads=args.T)

"""
