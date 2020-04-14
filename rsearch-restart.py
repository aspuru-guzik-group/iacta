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
                        help="Path to user.yaml file.",
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
    parser.add_argument("-p", "--params", help="File containing default parameters, i.e. those not set in user.yaml. Defaults to parameters/default.yaml.",
                        type=str, default=None)

    args = parser.parse_args()

    # Load user parameters (or try at least)
    with open(args.user_params, "r") as f:
        user_params = yaml.load(f, Loader=yaml.Loader)
        
    # Prepare output files
    # --------------------
    if args.o:
        out_dir = args.o
    else:
        # TODO: This and the -w flag is bad, we should fix it
        out_dir = os.path.dirname(args.user_params)
        
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

    # copy user params
    shutil.copy(args.user_params, out_dir + "/user.yaml")

    # Get default parameters
    params_file = args.params
    if params_file is None:
        # use parameters/default.yaml
        params_file = os.path.dirname(__file__)\
            + "parameters/default.yaml"
    
    with open(params_file, "r") as f:
        default_params = yaml.load(f, Loader=yaml.Loader)
        

    rsearch(out_dir, params_file,
            log_level=args.log_level,
            nthreads=args.threads)
    





