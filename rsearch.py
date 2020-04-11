import numpy as np
import react
import xtb_utils
import io_utils
import os
import shutil
import argparse
from constants import hartree_ev, ev_kcalmol
import yaml

def rsearch(out_dir, log_level=0, nthreads=1):
    if "LOCALSCRATCH" in os.environ:
        scratch = os.environ["LOCALSCRATCH"]
    else:
        print("warning: $LOCALSCRATCH not set")
        scratch = "."

    # load parameters
    with open(out_dir + "/user.yaml", "r") as f:
        user_params = yaml.load(f, Loader=yaml.Loader)
    with open(out_dir + "/default.yaml", "r") as f:
        params = yaml.load(f, Loader=yaml.Loader)        

    # Merge, replacing defaults with user parameters
    for key,val in user_params.items():
        params[key] = val

    # Interpret log level
    if log_level>1:
        delete=False
    else:
        delete=True

    # Command log file
    if log_level >0:
        logfile = open(out_dir + "/commandlog", "a")
        logfile.write("--------------------------"
                      +"--------------------------------------\n")
    else:
        logfile = None



    # Initialize the xtb driver
    # -------------------------
    xtb = xtb_utils.xtb_driver(scratch=scratch,
                               delete=delete,
                               logfile=logfile)

    xtb.extra_args = ["--gfn",params["gfn"]]
    if params["etemp"]:
        xtb.extra_args += ["--etemp", params["etemp"]]
    if params["chrg"]:
        xtb.extra_args += ["--chrg", params["chrg"]]
    if params["uhf"]:
        xtb.extra_args += ["--uhf", params["uhf"]]    
    if params["solvent"]:
        xtb.extra_args += ["--gbsa", params["solvent"]]
        
    # Temporarily set -P to number of threads for the next, non-parallelizable
    # two steps.
    xtb.extra_args += ["-P", str(nthreads)]

    # Optimize starting geometry including wall
    # -----------------------------------------
    init1 = out_dir + "/init.xyz"
    with open(init1, "w") as f:
        f.write(params["xyz"])
    

    print("Optimizing initial geometry 📐...")
    opt = xtb.optimize(init1, init1,
                       level=params["optim"],
                       xcontrol={"wall":params["wall"]})
    opt()

    # Read result of optimization
    atoms, positions, E = io_utils.traj2npy(init1, index=0)
    print("    E₀    = %15.7f Eₕ" % E)
    Emax = E + params["ewin"] / (hartree_ev * ev_kcalmol)
    print("    max E = %15.7f Eₕ  (E₀ + %5.1f kcal/mol)" %
          (Emax,params["threshold"]))
     params["emax"] = Emax           # update parameters
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple driver for reaction search. Builds a parameter file in an output folder.",
        )

    # These parameters do not have defaults
    parser.add_argument("init_xyz",
                        help="Path to file containing the starting geometry.",
                        type=str)
    parser.add_argument("atoms",
                        help="Atoms that define the bond to be stretched, numbered according"
                        +" to init_xyz. (NOTE THIS IS 1-INDEXED)",
                        type=int, nargs=2)

    # These are run specific parameters
    parser.add_argument("-o",
                        help="Output folder. Defaults to \"output\"",
                        type=str, default="output")
    parser.add_argument("-w",
                        help="Overwrite output directory. Defaults to false.",
                        action="store_true")
    parser.add_argument("-t", "--threads",
                        help="Number of threads to use.",
                        type=int, default=1)
    parser.add_argument("--log-level",
                        help="Level of debug printout (see react.py for details).",
                        default=0, type=int)
    parser.add_argument("-p", "--params", help="File containing numerical parameters.",
                        type=str, default="parameters/default.yaml")


    # These parameters (and some more!) have defaults in parameters/default.yaml.
    parser.add_argument("--optim", help="Optimization level.", type=str)

    parser.add_argument("-s","--stretch-limits", help="Bond stretch limits.", nargs=2,type=float)
    parser.add_argument("-r","--stretch-resolution", help="Number of bond stretches.", type=int)
    parser.add_argument("-k","--force-constant", help="Force constant of the stretch.", type=float)

    parser.add_argument("--gfn", help="gfn version.", type=str)
    parser.add_argument("--etemp", help="Electronic temperature.", type=str)
    parser.add_argument("--solvent", help="GBSA solvent.", type=str)
    parser.add_argument("-c", "--chrg", help="Charge.", type=str)
    parser.add_argument("-u", "--uhf", help="Spin state", type=str)

    args = parser.parse_args()

    # Prepare output files
    # --------------------
    out_dir = args.o
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        print("Output directory exists:")
        if args.w:
            # Delete the directory, make it and restart
            print("   👍 but that's fine! -w flag is on.")
            print("   📁 %s is overwritten."% args.o)
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)
        else:
            print("   👎 -w flag is off -> exiting! 🚪")
            raise SystemExit(-1)


    # Get default parameters
    params_file = args.params
    with open(params_file, "r") as f:
        default_params = yaml.load(f, Loader=yaml.Loader)

    # Save user-set parameters for reproducibility.
    user_params = {}
    args_as_dict = vars(args)
    for p in default_params.keys():
        argvalue = args_as_dict.get(p, None)
        if argvalue:
            user_params[p] = argvalue

    # Load the xyz file
    xyz,E = io_utils.traj2str(args.init_xyz, index=0)
    with open(out_dir + "/user.yaml", "w") as f:
        # begin with some metadata
        yaml.dump(io_utils.metadata(),f)
        yaml.dump(user_params,f)
        # write xyz at the beginning by hand so that it's formatted nicely.
        f.write("xyz: |\n")
        for line in xyz.split("\n"):
            # indent properly
            f.write("  " + line.lstrip() + "\n")
        f.write("\n")

    # Additionally save defaults for reproducibility
    shutil.copy(params_file, out_dir)

    rsearch(out_dir, log_level=args.log_level,
            nthreads=args.threads)
    
"""    



# Get bond parameters
# -------------------
bond_length0 = np.sqrt(np.sum((positions[params["atoms"][0]-1] -
                               positions[params["atoms"][1]-1])**2))
bond = (params["atoms"][0], params["atoms"][1], bond_length0)

# Constraints for the search
# -------------------------
stretch_factors = np.linspace(params["s"][0], params["s"][1], params["sn"])
print("Stretching bond between atoms %s%i and %s%i"
      %(atoms[bond[0]-1], bond[0], atoms[bond[1]-1], bond[1]))
print("    with force constant 💪💪 %f" % params["force"])
print("    between 📏 %7.2f and %7.2f A (%4.2f to %4.2f x bond length)"
      % (min(stretch_factors)*bond[2], max(stretch_factors)*bond[2],
         min(stretch_factors), max(stretch_factors)))
print("    discretized with %i points" % len(stretch_factors))
constraints = [("force constant = %f" % params["force"],
                "distance: %i, %i, %f"% (bond[0],bond[1],
                                         stretch * bond[2]))
               for stretch in stretch_factors]

# STEP 1: Initial generation of guesses
# ----------------------------------------------------------------------------
if do_stretch:
    react.generate_initial_structures(
        xtb, out_dir,
        init1,
        constraints,
        params)

# reset threading
xtb.extra_args = xtb.extra_args[:-2]

mtd_indices = params["mtdi"]
if mtd_indices is None:
    # Read the successive optimization, then set mtd points to ground and TS
    # geometries.
    reactant, E = io_utils.traj2smiles(init0, index=0)
    init, E = io_utils.traj2smiles(out_dir + "/init/opt.xyz")
    mtd_indices = [i for i,smi in enumerate(init) if smi==reactant][::3]
    print("Reactant 👉", reactant)

# Sort the indices, do not do the same point twice.
mtd_indices = sorted(list(set(mtd_indices)))
if len(mtd_indices) == 0:
    print("Reactant not found in initial stretch! 😢")
    print("Optimization probably reacted. Alter geometry and try again.")
    raise SystemExit(-1)


# STEP 2: Metadynamics
# ----------------------------------------------------------------------------
if do_mtd:
    react.metadynamics_search(
        xtb, out_dir,
        mtd_indices,
        constraints,
        params,
        nthreads=args.T)

if do_mtd_refine:
    react.metadynamics_refine(
        xtb, out_dir,
        init1,
        mtd_indices,
        constraints,
        params,
        nthreads=args.T)

# STEP 2: Reactions
# ----------------------------------------------------------------------------
if do_reactions:
    react.react(
        xtb, out_dir,
        mtd_indices,
        constraints,
        params,
        nthreads=args.T)


if logfile:
    logfile.close()
"""
