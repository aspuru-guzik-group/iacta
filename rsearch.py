import numpy as np
import os
import react
import read_reactions
import xtb_utils
import io_utils
import argparse
from constants import hartree_ev, ev_kcalmol
import yaml
from datetime import datetime

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

    # Initialize the xtb driver
    # -------------------------
    xtb = xtb_utils.xtb_driver(scratch=scratch,
                               delete=delete)


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

    # Put conformers on disk
    conformers_file = out_dir + "/conformers.xyz"
    with open(conformers_file, "w") as f:
        f.write(params["xyz"])

    # Load individual conformers
    conformers, energies = io_utils.traj2str(conformers_file)

    # Main loop over conformers
    for iconformer, conformer in enumerate(conformers):
        print("="*30 + " React. conf. %i of %i " %(iconformer+1,len(conformers))
              + "="*30)
        # Make the conformer directory
        conformer_dir = out_dir + "/conformer-%i" % (iconformer+1)
        os.makedirs(conformer_dir)
        
        # Temporarily set -P to number of threads for the next, non-parallelizable
        # two steps.
        xtb.extra_args += ["-P", str(nthreads)]

        # Optimize starting geometry including wall
        # -----------------------------------------
        init0 = conformer_dir + "/init_raw.xyz"
        with open(init0, "w") as f:
            f.write(conformer)
        init1 = conformer_dir + "/init_opt.xyz"

        print("Optimizing initial geometry 📐...")
        opt = xtb.optimize(init0, init1,
                           level=params["optim"],
                           xcontrol={"wall":params["wall"],
                                     # move to center of mass
                                     "cma":""})
        opt()

        # Read result of optimization
        atoms, positions, E = io_utils.traj2npy(init1, index=0)
        print("    E₀    = %15.7f Eₕ" % E)
        Emax = E + params["ewin"] / (hartree_ev * ev_kcalmol)
        print("    max E = %15.7f Eₕ  (E₀ + %5.1f kcal/mol)" %
              (Emax,params["ewin"]))
        params["emax"] = Emax           # update parameters


        # Get bond parameters
        # -------------------
        bond_length0 = np.sqrt(np.sum((positions[params["atoms"][0]-1] -
                                       positions[params["atoms"][1]-1])**2))

        # Constraints for the search
        # -------------------------
        slow, shigh = params["stretch_limits"]
        npts = params["stretch_num"]
        low = slow * bond_length0
        high = shigh * bond_length0
        atom1, atom2 = params["atoms"]

        print("Stretching bond between atoms %s%i and %s%i"
              %(atoms[atom1-1], atom1, atoms[atom2-1], atom2))
        print("    with force constant 💪💪 %f" % params["force"])
        print("    between 📏 %7.2f and %7.2f A (%4.2f to %4.2f x bond length)"
              % (low, high, slow, shigh))
        print("    discretized with %i points" % npts)


        # STEP 1: Initial generation of guesses
        # ----------------------------------------------------------------------------
        react.generate_initial_structures(
            xtb, conformer_dir, init1,
            atom1, atom2, low, high, npts,
            params, verbose=False)

        # reset threading
        xtb.extra_args = xtb.extra_args[:-2]

        # Read the successive optimization, then set mtd points to ground and TS
        # geometries.
        reactant, E0 = io_utils.traj2smiles(init1, index=0)
        init, E = io_utils.traj2smiles(conformer_dir + "/init/opt.xyz")
        reaction = read_reactions.read_reaction(conformer_dir + "/init")
        E = np.array(E)
        print("Reactant 👉", reactant)
        print("Molecules 👇")
        for i in range(len(reaction["E"])):
            if reaction["is_stable"][i]:
                print("%3i  %+7.3f -> %s" % (reaction["stretch_points"][i],
                                              reaction["E"][i],
                                              reaction["SMILES"][i]))
            else:
                print("%3i  %+7.3f     ...  ⛰  ..." %
                      (reaction["stretch_points"][i], reaction["E"][i]))

        if params["mtdi"]:
            mtd_indices = params["mtdi"]
        else:
            mtd_indices = [k for k in reaction["stretch_points"]]

            # additional indices at repeated intervals
            step = params["mtd_step"]
            if step:
                mtd_indices += list(np.arange(0,len(E), step))

            if params["mtd_only_reactant"]:
                mtd_indices = [i for i in mtd_indices if init[i] == reactant]
                print("     ... metadynamics performed only for reactants")

                if len(mtd_indices) == 0:
                    print("Reactant not found in initial stretch! 😢")
                    print("Optimization probably reacted. Alter geometry and try again.")
                    raise SystemExit(-1)

                # Also do the steps just before and just after
                mtd_indices += [max(mtd_indices) + 1,
                                min(mtd_indices) - 1]

            # Sort the indices, do not do the same point twice, make sure the points
            # are in bound
            mtd_indices = sorted(list(set([i for i in mtd_indices
                                           if (i >= 0 and i < len(E))])))



        # STEP 2: Metadynamics
        # ----------------------------------------------------------------------------
        react.metadynamics_search(
            xtb, conformer_dir,
            mtd_indices,
            atom1, atom2, low, high, npts,
            params,
            nthreads=nthreads)

        react.metadynamics_refine(
            xtb, conformer_dir,
            init1,
            mtd_indices,
            atom1, atom2, low, high, npts,
            params,
            nthreads=nthreads)

        # STEP 3: Reactions
        # ----------------------------------------------------------------------------
        react.react(
            xtb, conformer_dir,
            mtd_indices,
            atom1, atom2, low, high, npts,        
            params,
            nthreads=nthreads)
    print("="*60)
    print("... dumping run parameters ...")
    time_end = datetime.today().ctime()
    with open(out_dir + "/run.yaml", "w") as f:
        # begin with some metadata
        meta = io_utils.metadata()
        meta["start"] = time_start
        meta["end"] =time_end
        meta["nthreads"] = nthreads
        yaml.dump(io_utils.metadata(),f)
        
        # dump every other parameter
        yaml.dump(params,f)


    print("minigabe OUT. No more work to do! 🧔🍻\n\n")
