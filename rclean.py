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
workdir = args.folder

# TODO: refine a specific reactant
reactant = xyz2smiles(workdir + "/init/opt0000.xyz")
threshold = 60.0

# Prepare output files
# --------------------
# TODO: make an overwite? argument
os.makedirs(workdir + "/reactions", exist_ok=True)

# Initialize the xtb driver
# -------------------------
xtb = xtb_utils.xtb_driver(scratch=scratch)
xtb.extra_args = ["--gfn " + args.gfn, "--etemp " + args.etemp]

# Read all pathways and species
# -----------------------------
reactions = read_all_pathways(workdir, args.b * kcal2hartree)


import tempfile
import react_utils
from concurrent.futures import ThreadPoolExecutor

optimizing = set()
nthreads = 4
jobs = []
def opt_a_geom(xtb,folder,index,workdir):
    # Load the XYZ file and save it in the right spot
    fdc, current = tempfile.mkstemp(suffix=".xyz", dir=xtb.scratchdir)
    xyzprod = react_utils.read_trajectory(folder+"/opt.xyz", index)
    with open(current,"w") as f:
        f.write(xyzprod)
    # optimize
    return xtb.optimize(current,
                        workdir+"/reactions/opt_prod_%5.5i.xyz"%k,
                               level="vtight")

with ThreadPoolExecutor(max_workers=nthreads) as pool:
    for k, row in reactions.iterrows():
        print(k,row)
        folder = row.folder
        if (folder,row.prodI) in optimizing:
            pass                    # we already did this one
        else:
            # Add to set
            optimizing.add((folder, row.prodI))
            jobs += [opt_a_geom(xtb, folder, row.prodI, workdir)]
            pool.submit(jobs[-1])

        if (folder,row.reactI) in optimizing:
            pass                    # we already did this one
        else:
            # Add to set
            optimizing.add((folder, row.reactI))
            jobs += [opt_a_geom(xtb, folder, row.reactI, workdir)]
            pool.submit(jobs[-1])

# # reactions forward
# forward = pd.DataFrame(dict(
#     reactSMILES=reactions.reactSMILES,
#     prodSMILES=reactions.prodSMILES,
#     folder=reactions.folder,
#     barrier=(reactions.tsE - reactions.reactE)*hartree2ev*ev2kcal))

# backward = pd.DataFrame(dict(
#     prodSMILES=reactions.reactSMILES,
#     reactSMILES=reactions.prodSMILES,
#     folder=reactions.folder,
#     barrier=(reactions.tsE - reactions.prodE)*hartree2ev*ev2kcal))
# allreacts = forward.append(backward)

# print("Summary of approximate reaction barriers (kcal / mol)")
# pivot = allreacts.pivot_table(values=["barrier"],
#                               index=["reactSMILES", "prodSMILES"],
#                               aggfunc=["min", "mean", "max"])
# print(pivot)


# print("Chosen reactant : ", reactant)
# print("Products with barrier < %9.2f:" % threshold)
# product_pivot = pivot.loc[reactant].reset_index()
# products = product_pivot.prodSMILES[product_pivot[("min","barrier")]<threshold]
# print(products)






"""


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
