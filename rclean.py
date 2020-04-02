import numpy as np
import react
import xtb_utils
import io_utils
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
                    +" Defaults to 0.5 kcal/mol.", default=0.5,
                    type=float)
parser.add_argument("-T",
                    help="Number of threads to use.",
                    type=int, default=1)
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

rmsd_threshold = 0.2
reactions = read_all_pathways(workdir, 0.2)

    











"""
# TODO: refine a specific reactant
reactant = xyz2smiles(workdir + "/init/opt0000.xyz")
threshold = 60.0
nthreads = 4
# Prepare output files
# --------------------

# Initialize the xtb driver
# -------------------------
xtb = xtb_utils.xtb_driver(scratch=scratch)
xtb.extra_args = ["--gfn " + args.gfn, "--etemp " + args.etemp]

# Read all pathways and species
# -----------------------------


# reactions forward
forward = pd.DataFrame(dict(
    reactSMILES=reactions.reactSMILES,
    prodSMILES=reactions.prodSMILES,
    folder=reactions.folder,
    barrier=(reactions.tsE - reactions.reactE)*hartree2ev*ev2kcal))

backward = pd.DataFrame(dict(
    prodSMILES=reactions.reactSMILES,
    reactSMILES=reactions.prodSMILES,
    folder=reactions.folder,
    barrier=(reactions.tsE - reactions.prodE)*hartree2ev*ev2kcal))
allreacts = forward.append(backward)

print("Summary of *approximate* reaction barriers (kcal / mol)")
pivot = allreacts.pivot_table(values=["barrier"],
                              index=["reactSMILES", "prodSMILES"],
                              aggfunc=["min", "mean", "max"])
print(pivot)


# print("Chosen reactant : ", reactant)
# print("Products with barrier < %9.2f:" % threshold)
# product_pivot = pivot.loc[reactant].reset_index()
# products = product_pivot.prodSMILES[product_pivot[("min","barrier")]<threshold]
# print(products)

from concurrent.futures import ThreadPoolExecutor
import react_utils
# Optimize and sort
# -----------------
species = list(set([s for s in reactions.prodSMILES.values]
              +[s for s in reactions.reactSMILES.values]))
print("%i unique SMILES molecules"% len(species))


# Make the SMILES file
with open(workdir+"/smiles", "w") as f:
    for i in range(len(species)):
        f.write(species[i]+"\n")


optimizing = {}
os.makedirs(workdir + "/mols")
for k, row in reactions.iterrows():
    folder = row.folder
    index = row.prodI
    smiles = row.prodSMILES
    for folder,index,smiles in [(row.folder, row.prodI, row.prodSMILES),
                                (row.folder, row.reactI, row.reactSMILES)]:
        if (folder,index) in optimizing:
            pass                    # we already did this one
        else:
            xyzprod = react_utils.read_trajectory(folder+"/opt.xyz", index)

            # Add to set of stuff we found already
            optimizing[(folder, index)] = k

            # also add it to the right file
            with open(workdir+"/mols/%5.5i.xyz"%species.index(smiles), "a") as f:
                f.write(xyzprod)
f.close()

from concurrent.futures import ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=nthreads) as pool:
    for i in range(len(species)):
        pool.submit(
            xtb.screen(workdir+"/mols/%5.5i.xyz" % i,
                       workdir+"/mols/%5.5i.xyz" % i))

# Now that this is done







if not args.no_opt:
    xtb.optimize(init, init, level="vtight")

# Get additional molecular parameters
# -----------------------------------
atoms, positions = io_utils.read_xyz(init)
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
