import os
import shutil
import argparse
import yaml

# ========================== CLI INTERFACE ================================
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
                        type=str, default=None)
    parser.add_argument("-d", "--dump",
                        help="Make output directory and save user parameters, but do not"
                        +" search for reaction. Such files can be run using rsearch-restart.py, or scripted in python using the rsearch() function.",
                        action="store_true")

    # These parameters (and some more!) have defaults in parameters/default.yaml.
    parser.add_argument("--optim", help="Optimization level.", type=str)

    parser.add_argument("-s","--stretch-limits", help="Bond stretch limits.", nargs=2,type=float)
    parser.add_argument("-n","--stretch-num", help="Number of bond stretches.", type=int)    
    parser.add_argument("-k","--force", help="Force constant of the stretch.", type=float)

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
    if params_file is None:
        # use parameters/default.yaml
        params_file = os.path.dirname(__file__)\
            + "/parameters/default.yaml"
    
    with open(params_file, "r") as f:
        default_params = yaml.load(f, Loader=yaml.Loader)

    user_params = {}
    args_as_dict = vars(args)
    for p in default_params.keys():
        argvalue = args_as_dict.get(p, None)
        if argvalue:
            user_params[p] = argvalue

    # Load the xyz file
    xyz,E = io_utils.traj2str(args.init_xyz, index=0)

    # Save everything in the yaml file
    with open(out_dir + "/user.yaml", "w") as f:
        yaml.dump(user_params,f)
        # write xyz at the beginning by hand so that it's formatted nicely.
        f.write("xyz: |\n")
        for line in xyz.split("\n"):
            # indent properly
            f.write("  " + line.lstrip() + "\n")
        f.write("\n")
        






