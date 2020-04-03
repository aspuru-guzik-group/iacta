import numpy as np
import react
import xtb_utils
import os
import shutil
import argparse

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
                    +" Defaults to one point every 10.",
                    type=int, nargs="+", default=None)
parser.add_argument("-mtdn",
                    help="Number of guesses to generate at each MTD index.",
                    type=int, default=80)
parser.add_argument("-force",
                    help="Force constant of the stretch."
                    +" Defaults to 1.00 Eh/Bohr.",
                    default=1.0,
                    type=float)
parser.add_argument("-no-opt",
                    help="Start with an xtb optimization (defaults to true).",
                    action="store_true")
parser.add_argument("-gfn",
                    help="gfn version. Defaults to GFN 2", default="2",
                    type=str)
parser.add_argument("-etemp",
                    help="Electronic temperature. Defaults to 300 K",
                    default="300.0",
                    type=str)
parser.add_argument("-opt-level",
                    help="Optimization level. Defaults to vtight.",
                    default="vtight",
                    type=str)
parser.add_argument("-shake-level",
                    help="If this is 0, the metadynamics run will be performed"
                    +" with parameters that permit bond-breaking, specifically shake=0,"
                    +" but at a slower pace than if shake is set to 1 or 2. Defaults to 0.",
                    default=0, type=int)
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
os.makedirs(out_dir)
init = shutil.copy(args.init_xyz,
                   out_dir)

# Initialize the xtb driver
# -------------------------
if args.log_level>1:
    delete=False
else:
    delete=True
xtb = xtb_utils.xtb_driver(scratch=scratch,
                           delete=delete)
xtb.extra_args = ["--gfn " , args.gfn, "--etemp " , args.etemp]

# Get additional molecular parameters
# -----------------------------------
atoms, positions = xtb_utils.read_xyz(init)
N = len(atoms)

# Initialize parameters
# ---------------------
params = react.default_parameters(N,
                                  shake=args.shake_level,
                                  nmtd=args.mtdn,
                                  optlevel=args.opt_level,
                                  log_level=args.log_level)

if not args.no_opt:
    print("optimizing initial geometry...")
    opt = xtb.optimize(init, init, level="vtight",
                       xcontrol={"wall":params["wall"]})
    opt()


bond_length0 = np.sqrt(np.sum((positions[args.atoms[0]-1] -
                               positions[args.atoms[1]-1])**2))
bond = (args.atoms[0], args.atoms[1], bond_length0)


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
if mtd_indices is None:
    mtd_indices = np.arange(0, args.sn, 10)


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
