import numpy as np
import react
import xtb_utils
import io_utils
import os
import shutil
import argparse
from constants import hartree_ev, ev_kcalmol

parser = argparse.ArgumentParser(
    description="Driver for reaction search.",
    )
parser.add_argument("init_xyz",
                    help="Path to file containing the starting geometry.",
                    type=str)
parser.add_argument("atoms",
                    help="Atoms that define the bond to be stretched, numbered according"
                    +" to init_xyz. (NOTE THIS IS 1-INDEXED)",
                    type=int, nargs=2)
parser.add_argument("-o",
                    help="Output folder. defaults to \"output\"",
                    type=str, default="output")
parser.add_argument("-w",
                    help="Overwrite output directory. Defaults to false.",
                    action="store_true")
parser.add_argument("-T",
                    help="Number of threads to use.",
                    type=int, default=1)
parser.add_argument("-s",
                    help="Bond stretch limits. Defaults to (1.0, 3.0)",
                    nargs=2,type=float, default=[1.0,3.0])
parser.add_argument("-sn",
                    help="Number of bond stretches. Defaults to 300.",
                    type=int, default=100)
parser.add_argument("-mtdi",
                    help="Indices of the stretches where MTD is done."
                    +" Defaults $are complicated.",
                    type=int, nargs="+", default=None)
parser.add_argument("-mtdn",
                    help="Number of guesses to generate at each MTD index."
                    +" Defaults to 120.",
                    type=int, default=120)
parser.add_argument("-force",
                    help="Force constant of the stretch."
                    +" Defaults to 1.00 Eₕ / Bohr.",
                    default=1.0,
                    type=float)
parser.add_argument("-gfn",
                    help="gfn version. Defaults to GFN 2", default="2",
                    type=str)
parser.add_argument("-solvent",
                    help="Set GBSA solvent.", 
                    type=str)
parser.add_argument("-chrg",
                    help="Set charge for xtb.", default=None,
                    type=str)
parser.add_argument("-uhf",
                    help="Set spin state for xtb.",
                    default=None,
                    type=str)
parser.add_argument("-etemp",
                    help="Electronic temperature. Defaults to xtb default",
                    default=None,
                    type=str)
parser.add_argument("-opt-level",
                    help="Optimization level. Defaults to vtight.",
                    default="tight",
                    type=str)
parser.add_argument("-threshold",
                    help="Energy threshold for path optimization in kcal/mol."
                    +" Basically, if a barrier is encountered that is higher than"
                    +" this threshold from the optimized reactant energy, the"
                    +" entire path is discarded. Defaults to 50 kcal/mol. Note"
                    +" that some barriers higher than this value will still be"
                    " found, due to implementation details.",
                    default=50.0,
                    type=float)
parser.add_argument("-shake-level",
                    help="If this is 0, the metadynamics run will be performed"
                    +" with parameters that permit bond-breaking, specifically"
                    +" shake=0, but at a slower pace than if shake is set to 1"
                    +" or 2. Defaults to 0.",
                    default=2, type=int)
parser.add_argument("-log-level",
                    help="Level of debug printout (see react.py for details).",
                    default=0, type=int)


if "LOCALSCRATCH" in os.environ:
    scratch = os.environ["LOCALSCRATCH"]
else:
    print("warning: $LOCALSCRATCH not set")
    scratch = "."

args = parser.parse_args()

# Prepare output files
# --------------------
out_dir = args.o
try:
    os.makedirs(out_dir)
except FileExistsError:
    print("😭 Output directory exists:")
    if args.w:
        # Delete the directory, make it and restart
        print("   -w flag is on -> %s is overwritten ✏."% args.o)
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
    else:
        print("   -w flag is off -> exiting! 🚪")
        raise SystemExit(1)
        
init = shutil.copy(args.init_xyz,
                   out_dir)
atoms, positions, comment = io_utils.read_xyz(init)
Natoms = len(atoms)

# Initialize the xtb driver
# -------------------------
if args.log_level>1:
    delete=False
else:
    delete=True
xtb = xtb_utils.xtb_driver(scratch=scratch,
                           delete=delete)
xtb.extra_args = ["--gfn",args.gfn]
if args.etemp:
    xtb.extra_args += ["--etemp", args.etemp]
if args.chrg:
    xtb.extra_args += ["--chrg", args.chrg]
if args.uhf:
    xtb.extra_args += ["--uhf", args.uhf]    
if args.solvent:
    xtb.extra_args += ["--gbsa", args.solvent]

# Initialize parameters
# ---------------------
ethreshold = args.threshold / (hartree_ev * ev_kcalmol)
params = react.default_parameters(Natoms,
                                  shake=args.shake_level,
                                  nmtd=args.mtdn,
                                  ethreshold=ethreshold,
                                  optlevel=args.opt_level,
                                  log_level=args.log_level)

# Temporarily set -P to number of threads for the next, non-parallelizable two
# steps.
xtb.extra_args += ["-P", str(args.T)]

# Optimize starting geometry including wall
# -----------------------------------------
print("Optimizing initial geometry 📐...")
opt = xtb.optimize(init, init,
                   level=args.opt_level,
                   xcontrol={"wall":params["wall"]})
opt()

# Read result of optimization
atoms, positions, comment = io_utils.read_xyz(init)
E = io_utils.comment_line_energy(comment)
print("👍 Done!    E₀ = %15.7f Eₕ" % E)

# Get bond parameters
# -------------------
bond_length0 = np.sqrt(np.sum((positions[args.atoms[0]-1] -
                               positions[args.atoms[1]-1])**2))
bond = (args.atoms[0], args.atoms[1], bond_length0)

# Constraints for the search
# -------------------------
stretch_factors = np.linspace(args.s[0], args.s[1], args.sn)
print("Stretching bond between atoms %s%i and %s%i"
      %(atoms[bond[0]-1],bond[0], atoms[bond[1]-1],bond[1]))
print("    with force constant 💪💪 %f" % args.force)
print("    between 📏 %7.2f and %7.2f A (%4.2f to %4.2f x bond length)"
      % (min(stretch_factors)*bond[2], max(stretch_factors)*bond[2],
         min(stretch_factors), max(stretch_factors)))
print("    discretized with %i points" % len(stretch_factors))
print("    with a maximum barrier height %7.5f Eₕ (%7.3f kcal/mol) " %
      (ethreshold, args.threshold))
constraints = [("force constant = %f" % args.force,
                "distance: %i, %i, %f"% (bond[0],bond[1],
                                         stretch * bond[2]))
               for stretch in stretch_factors]

# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
react.generate_initial_structures(
    xtb, out_dir,
    init,
    constraints,
    params)

# reset threading
xtb.extra_args = xtb.extra_args[:-2]

mtd_indices = args.mtdi
if mtd_indices is None:
    # Read the successive optimization, then set mtd points to ground and TS
    # geometries.
    from read_reactions import read_reaction, xyz2smiles
    out = read_reaction(out_dir + "/init")
    reactant=xyz2smiles(args.init_xyz)[0]
    init = xyz2smiles(out_dir + "/init/opt.xyz")
    mtd_indices = [i for i,smi in enumerate(init) if smi==reactant][::3]
    print("Reactant 👉", reactant)

# Sort the indices, do not do the same point twice.
mtd_indices = sorted(list(set(mtd_indices)))
if len(mtd_indices) == 0:
    raise SystemExit(-1)


# STEP 2: Metadynamics
# ----------------------------------------------------------------------------
react.metadynamics_search(
    xtb, out_dir,
    mtd_indices,
    constraints,
    params,
    nthreads=args.T)

react.metadynamics_refine(
    xtb, out_dir,
    init,
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

