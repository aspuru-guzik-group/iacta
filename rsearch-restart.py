import os
import argparse
from rsearch import rsearch
import yaml
import shutil

# ========================== CLI INTERFACE ================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple driver for reaction search. Runs reactions for a  user.yaml parameter file.")

    # These parameters do not have defaults
    parser.add_argument("user_params",
                        help="Path to user.yaml file or path to run directory "+
                        "containing the parameter file.",
                        type=str)

    # These are run specific parameters
    parser.add_argument("-o",
                        help="Output folder. Defaults to the same directory as user_params",
                        type=str, default=None)
    parser.add_argument("-w",
                        help="Overwrite output directory. Defaults to false.",
                        action="store_true")
    parser.add_argument("-t", "--threads",
                        help="Number of threads to use.",
                        type=int, default=1)
    parser.add_argument("--log-level",
                        help="Level of debug printout (see react.py for details).",
                        default=0, type=int)
    parser.add_argument("-p", "--params", help="File containing default parameters, i.e. those not set in user.yaml. Defaults to parameters/default.yaml in the minigabe directory.",
                        type=str, default=None)

    args = parser.parse_args()

    # Load user parameters (or try at least)
    try:
        pfile = args.user_params
        with open(pfile, "r") as f:
            user_params = yaml.load(f, Loader=yaml.Loader)
    except IsADirectoryError:
        pfile = args.user_params + "/user.yaml"
        with open(pfile, "r") as f:
            user_params = yaml.load(f, Loader=yaml.Loader)

    # save user parameters
    with open(pfile, "r") as f:
        upfile = f.read()

    # Prepare output files
    # --------------------
    if args.o:
        out_dir = args.o
    else:
        # TODO: This and the -w flag is bad, we should fix it
        out_dir = os.path.dirname(pfile)

    try:
        os.makedirs(out_dir)
    except FileExistsError:
        print("Output directory exists:")
        if args.w:
            # Delete the directory, make it and restart
            print("   👍 but that's fine! -w flag is on.")
            print("   📁 %s is overwritten."% out_dir)
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)
        else:
            print("   👎 -w flag is off -> exiting! 🚪")
            raise SystemExit(-1)

    # copy user params
    with open(out_dir +"/user.yaml", "w") as f:
        f.write(upfile)

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

    out = rsearch(out_dir, params_file,
                  log_level=args.log_level,
                  nthreads=args.threads)
