import numpy as np
import react
from react_utils import stretch
from analysis import postprocess_reaction
import xtb_utils
import io_utils
from io_utils import pybel
import os
import shutil
import argparse
from constants import hartree_ev, ev_kcalmol, bohr_ang
import yaml
from datetime import datetime

def cval(mol, atoms_i):
    atoms = [mol.GetAtom(i) for i in atoms_i]
    if len(atoms)==2:
        return mol.GetBond(*atoms).GetLength()
    if len(atoms)==3:
        return mol.GetAngle(*atoms)
    if len(atoms)==4:
        return mol.GetTorsion(*atoms)

def init_xtb_driver(params, log_level=0):
    # todo : move this stuff to xtb_driver
    if "LOCALSCRATCH" in os.environ:
        scratch = os.environ["LOCALSCRATCH"]
    else:
        print("warning: $LOCALSCRATCH not set")
        scratch = "."

    # Interpret log level
    if log_level>1:
        delete=False
    else:
        delete=True

    # Command log file
    if log_level >0:
        logfile = None          # TODO fix?
        # logfile = open(out_dir + "/commandlog", "a")
        # logfile.write("--------------------------"
        #               +"--------------------------------------\n")
    else:
        logfile = None

    # Initialize the xtb driver
    # -------------------------
    xtb = xtb_utils.xtb_driver(scratch=scratch,
                               delete=delete,
                               logfile=logfile)

    xtb.extra_args = ["--gfn",str(params["gfn"])]
    if params["etemp"]:
        xtb.extra_args += ["--etemp", str(params["etemp"])]
    if params["chrg"]:
        xtb.extra_args += ["--chrg", str(params["chrg"])]
    if params["uhf"]:
        xtb.extra_args += ["--uhf", str(params["uhf"])]
    if params["solvent"]:
        xtb.extra_args += ["--gbsa", params["solvent"]]
    return xtb

def rsearch(out_dir, defaults,
            log_level=0, nthreads=1):

    time_start = datetime.today().ctime()


    # load parameters
    with open(out_dir + "/user.yaml", "r") as f:
        user_params = yaml.load(f, Loader=yaml.Loader)
    with open(defaults, "r") as f:
        params = yaml.load(f, Loader=yaml.Loader)

    # Merge, replacing defaults with user parameters
    for key,val in user_params.items():
        params[key] = val

    xtb = init_xtb_driver(params, log_level=log_level)

    # Temporarily set -P to number of threads for the next, non-parallelizable
    # two steps.
    xtb.extra_args += ["-P", str(nthreads)]

    # Optimize starting geometry including wall
    # -----------------------------------------
    init0 = out_dir + "/init_raw.xyz"
    with open(init0, "w") as f:
        f.write(params["xyz"])
    init1 = out_dir + "/init_opt.xyz"


    if params["cavity_radius"]:
        radius_bohr = params["cavity_radius"] / bohr_ang
    else:
        # Load the molecule and compute its radius for the wall size
        at, pos = io_utils.xyz2numpy(params["xyz"])
        # Compute all interatomic distances
        distances = []
        for i in range(len(at)):
            for j in range(i):
                distances += [np.sqrt(np.sum((pos[i] - pos[j])**2))]

        # Cavity is 1.5 x maximum distance
        radius_bohr = max(distances) * 1.5 / bohr_ang

    print("Size of constraining cavity: %f A" % (radius_bohr * bohr_ang))

    params["wall"] = ["potential=logfermi",
                        "sphere=%f, all" % radius_bohr]



    print("Optimizing initial geometry...")
    opt = xtb.optimize(init0, init1,
                       level=params["optim"],
                       xcontrol={"wall":params["wall"],
                                 # move to center of mass
                                 "cma":""})
    opt()

    # Read result of optimization and set the maximum energy.
    mol, E = io_utils.traj2mols(init1, index=0)
    print("    E₀    = %15.7f Eₕ" % E)
    params["E0"] = E
    Emax = E + params["emax_global"] / (hartree_ev * ev_kcalmol)
    print("    max E = %15.7f Eₕ  (E₀ + %5.1f kcal/mol)" %
          (Emax,params["emax_global"]))


    # Get constraints parameters
    # --------------------------
    atoms = params["atoms"]
    ob_at = [mol.GetAtom(at) for at in atoms]
    current = cval(mol, atoms)

    try:
        low, high = params["driving_limits"]
    except TypeError:
        high = params["driving_limits"]
        low = current
    npts = params["driving_num"]
    print("\n")
    print("+--------------------------+")
    print("|   *Coordinate Driving*   |")

    if len(atoms)==2:
        print("|   Interatomic distance   |")
        print("+--------------------------+")
        print("Atoms: %s#%i --- %s#%i" % (ob_at[0].GetType(), atoms[0],
                                        ob_at[1].GetType(), atoms[1]) )
        print("\n            from: %6.2f Å" % low)
        print("              to: %6.2f Å" % high)
        print("             opt: %6.2f Å" % current)
        print("          nsteps: %i" % npts)
    if len(atoms)==3:
        print("|      Bending angle       |")
        print("+--------------------------+")
        print("Atoms: %s#%i      %s#%i" % (ob_at[0].GetType(), atoms[0],
                                         ob_at[2].GetType(), atoms[2]) )
        print("          \     /")
        print("            %s#%i" % (ob_at[1].GetType(), atoms[1]) )
        print("\n            from: %6.2f°" % low)
        print("              to: %6.2f°" % high)
        print("             opt: %6.2f°" % current)
        print("          nsteps: %i" % npts)
    if len(atoms)==4:
        print("|      Torsion angle       |")
        print("+--------------------------+")
        print("Atoms: %s#%i     " % (ob_at[0].GetType(), atoms[0]))
        print("          \     ")
        print("           %s#%i -- %s#%i" % (ob_at[1].GetType(), atoms[1],
                                           ob_at[2].GetType(), atoms[2]) )
        print("                     \    ")
        print("                       %s#%i" % (ob_at[3].GetType(), atoms[3]))
        print("\n            from: %6.2f°" % low)
        print("              to: %6.2f°" % high)
        print("             opt: %6.2f°" % current)
        print("          nsteps: %i" % npts)


    # Constraints for the search
    # -------------------------
    if not params['force']:
        # we do so quite simply from a 4 points polynomial fit
        params['force'] = 5.0
        if len(atoms) == 2:
            x0 = np.linspace(current - 0.05, current + 0.05, 5)
            structs, y = stretch(
                xtb, init1,
                atoms,
                x0[0], x0[-1], len(x0),
                params,
                verbose=True)

            mols = [pybel.readstring("xyz", s.lower()).OBMol for s in structs]
            x = [abs(cval(mol, atoms)) for mol in mols]
            x = np.array(x)

            y = np.array(y)
            p = np.polyfit(x, y, 2)
            k = 2*p[0]
            params["force"] = float(k * bohr_ang)
            print("    computed force constant 💪💪 %f" % params["force"])
        else:
            params["force"] = 1.0
            print("     default force constant 💪💪 %f" % params["force"])
    else:
        print("    with force constant 💪💪 %f" % params["force"])


    # STEP 1: Initial generation of guess conformers
    # ----------------------------------------------------------------------------
    react.generate_initial_structures(
        xtb, out_dir, init1,
        atoms, low, high, npts,
        params)

    # reset threading
    xtb.extra_args = xtb.extra_args[:-2]

    # Refinement and selection
    mtd_indices = react.select_initial_structures(
        xtb, out_dir, init1,
        atoms, low, high, npts,
        params, nthreads=nthreads)

    # STEP 2: Metadynamics
    # ----------------------------------------------------------------------------
    react.metadynamics_search(
        xtb, out_dir,
        mtd_indices,
        atoms, low, high, npts,
        params,
        nthreads=nthreads)

    react.metadynamics_refine(
        xtb, out_dir,
        init1,
        mtd_indices,
        atoms, low, high, npts,
        params,
        nthreads=nthreads)

    # STEP 3: Reactions
    # ----------------------------------------------------------------------------
    react.react(
        xtb, out_dir,
        mtd_indices,
        atoms, low, high, npts,
        params,
        nthreads=nthreads)

    # todo: re-integrate
    # if logfile:
    #     logfile.close()


    time_end = datetime.today().ctime()
    with open(out_dir + "/run.yaml", "w") as f:
        # begin with some metadata
        meta = io_utils.metadata()
        meta["start"] = time_start
        meta["end"] =time_end
        yaml.dump(io_utils.metadata(),f)
        # Every parameter and then some
        yaml.dump(params,f)
        # dump extra stuff
        yaml.dump({"nthreads":nthreads,
                   "done_metadynamics_pts":list(mtd_indices)})





# ========================== CLI INTERFACE ================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple driver for reaction search. Builds a parameter file in an output folder.",
        )

    # These parameters do not have defaults
    parser.add_argument("init_xyz",
                        help="Path to file containing the starting geometry.",
                        type=str)
    parser.add_argument("atom1_atom2",
                        help="Atoms that define the coordinate to be driven. Two"
                        +" atoms define a stretch, three atoms define a bend and 4"
                        +" atoms define a torsion. NOTE: Atoms are numbered"
                        +" starting from 1, as is standard in chemistry.",
                        type=int, nargs=2)
    parser.add_argument("atom3",
                        type=int, nargs="?")
    parser.add_argument("atom4",
                        type=int, nargs="?")
    parser.add_argument("driving_to",
                        type=float,
                        help="driving-to should be the bond length in angstrom or the angle in degrees of the driving coordinate at the end of driving. The start of driving is the corresponding value at equilibrium for init_xyz (if --driving-from is not given).")

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
    parser.add_argument("-p", "--params", help="File containing default numerical parameters. Defaults to parameters/default.yaml in the minigabe directory.",
                        type=str, default=None)
    parser.add_argument("-d", "--dump",
                        help="Make output directory and save user parameters, but do not"
                        +" search for reaction. Such files can be run using rsearch-restart.py, or scripted in python using the rsearch() function.",
                        action="store_true")

    # These parameters (and some more!) have defaults in parameters/default.yaml.
    parser.add_argument("--driving-from",
                        help="Minimum value of the coordinate to be driven.", type=float)
    parser.add_argument("--optim",
                        help="Optimization level used during coordinate driving.", type=str)
    parser.add_argument("--no-initial-mtd",
                        help="Do not initialize from a metadynamics-derived set of structures.",
                        action="store_true")
    parser.add_argument("-n","--driving-num", help="Number of points for coordinate driving.", type=int)
    parser.add_argument("-k","--force",
                        help="Force constant of the driving (Hartree/bohr or Hartree/rad)."
                        +" Defaults to 1.0 for bends and torsion, and to the bond strength as"
                        +" calculated from a five point relaxed scan for stretches.",
                        type=float)
    parser.add_argument("--gfn", help="gfn version.", type=str)
    parser.add_argument("--solvent", help="GBSA solvent.", type=str)
    parser.add_argument("-c", "--chrg", help="Charge.", type=str)
    parser.add_argument("-u", "--uhf", help="Spin state", type=str)
    parser.add_argument("--etemp", help="Electronic temperature.", type=str)

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
    if params_file is None:
        # use parameters/default.yaml
        folder = os.path.dirname(__file__)
        if not folder:
            folder = "."
        params_file = folder \
            + "/parameters/default.yaml"

    with open(params_file, "r") as f:
        default_params = yaml.load(f, Loader=yaml.Loader)

    # Save user-set command line parameters for reproducibility.
    user_params = {}
    args_as_dict = vars(args)
    for p in default_params.keys():
        # arguments named the same as in the file
        argvalue = args_as_dict.get(p, None)
        if argvalue:
            user_params[p] = argvalue

    # other arguments that do not have the same structure has in the
    # default.yaml file
    if args.no_initial_mtd:
        user_params['imtd'] = False

    # atoms defining the driving coordinate
    user_params["atoms"] = args.atom1_atom2
    if not args.atom3 is None:
        user_params["atoms"] += [args.atom3]
    if not args.atom4 is None:
        user_params["atoms"] += [args.atom4]

    # driving limits
    if args.driving_from is None:
        user_params["driving_limits"] = args.driving_to
    else:
        user_params["driving_limits"] = [args.driving_from, args.driving_to]

    # Load the xyz file
    xyz,E = io_utils.traj2str(args.init_xyz, index=0)
    with open(out_dir + "/user.yaml", "w") as f:
        yaml.dump(user_params,f)
        # write xyz at the beginning by hand so that it's formatted nicely.
        f.write("xyz: |\n")
        for line in xyz.split("\n"):
            # indent properly
            f.write("  " + line.lstrip() + "\n")
        f.write("\n")

    if not args.dump:
        rsearch(out_dir, params_file,
                log_level=args.log_level,
                nthreads=args.threads)
