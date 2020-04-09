from concurrent.futures import ThreadPoolExecutor
import react_utils
from io_utils import read_trajectory, dump_succ_opt
from math import inf
import os
from glob import glob
import tempfile
from constants import hartree_ev, ev_kcalmol

"""
This file contains a bunch of user-friendly, multithreaded drivers for
react_utils.
"""

def default_parameters(Natoms,
                       nmtd=80,
                       optlevel="tight",
                       ethreshold=inf,
                       rthreshold=0.5,
                       shake=0,
                       log_level=0):
    """Generate a dictionary of default parameters for the reaction search.

    The parameters are loosely based on those described in the CREST paper.

    Parameters:
    -----------

    Natoms (int) : number of atoms.

    Optional Parameters:
    --------------------
    nmtd (int) : number of structures to extract from the metadynamics
    portion.

    optlevel (str) : Optimization level to used throughout. Defaults to
    "tight".

    ethreshold (float) : Energy threshold for the successive optimizations. If
    a barrier ever reaches this threshold, the optimization goes no further.

    shake (int) : SHAKE level (0, 1 or 2). If SHAKE is > 0, the time
    parameters are adjusted to go faster.

    log_level (int) : Logging level. 0 -> normal, 1-> intermediate
    optimization steps. 2 -> keep all temp files (TODO: not implemented yet).

    Returns:
    --------

    dict : parameters dictionary to pass to other functions in react.py

    """
    parameters = {}
    parameters["optlevel"] = optlevel
    parameters["ethreshold"] = ethreshold
    parameters["rthreshold"] = rthreshold

    # Logging
    if log_level > 0:
        parameters["log_opt_steps"] = True
    else:
        parameters["log_opt_steps"] = False

    
    # xTB additional parameters
    parameters["wall"] = ("potential=logfermi",
                          "sphere: auto, all")

    # Metadynamics parameters (from CREST)
    total_time = int(round(0.4 * Natoms))
    SAVE="save=%i" % nmtd

    # Parameters from mquick
    parameters["CREST"] = [
        (SAVE,"kpush=0.0200","alpha=1.000"),
        (SAVE,"kpush=0.0100","alpha=1.000"),
        (SAVE,"kpush=0.0200","alpha=0.500"),
        (SAVE,"kpush=0.0100","alpha=0.500"),
        (SAVE,"kpush=0.0200","alpha=0.250"),
        (SAVE,"kpush=0.0100","alpha=0.250"),]
    
    # parameters from normal
    # parameters["CREST"] = [
    #     (SAVE, "kpush=0.033", "alpha=1.300"),
    #     (SAVE, "kpush=0.017", "alpha=1.300"),
    #     (SAVE, "kpush=0.008", "alpha=1.300"),
    #     (SAVE, "kpush=0.033", "alpha=0.780"),
    #     (SAVE, "kpush=0.017", "alpha=0.780"),
    #     (SAVE, "kpush=0.008", "alpha=0.780"),
    #     (SAVE, "kpush=0.033", "alpha=0.468"),
    #     (SAVE, "kpush=0.017", "alpha=0.468"),
    #     (SAVE, "kpush=0.008", "alpha=0.468"),
    #     (SAVE, "kpush=0.033", "alpha=0.281"),
    #     (SAVE, "kpush=0.017", "alpha=0.281"),
    #     (SAVE, "kpush=0.008", "alpha=0.281"),
    #     (SAVE, "kpush=0.011", "alpha=0.100"),
    #     (SAVE, "kpush=0.055", "alpha=0.800")]

    if shake == 0:
        parameters["md"] = ("shake=0",
                            "step=1",
                            "dump=100.0",
                            "time=%f" % total_time)
    else:
        parameters["md"] = ("shake=%i" % shake,
                            "step=5",
                            "dump=100.0",
                            "time=%f" % total_time)

    return parameters

def generate_initial_structures(xtb_driver,
                                workdir,
                                guess_xyz_file,
                                constraints,
                                parameters,
                                verbose=True):
    
                                
    if verbose:
        print("-----------------------------------------------------------------")
        print("Performing initial stretching 💪😎💪...")

    outputdir = workdir + "/init"
    os.makedirs(outputdir)
    
    structures, energies, opt_indices = react_utils.successive_optimization(
        xtb_driver, guess_xyz_file,
        constraints, parameters, # no barrier for initial run
        failout=outputdir + "/FAILED",
        verbose=verbose)

    computed = len(opt_indices)

    dump_succ_opt(outputdir,
                  structures,
                  energies,
                  opt_indices,
                  extra=parameters["log_opt_steps"],
                  split=True)
    if verbose:
        print("Done!  %i steps were evaluated. \n"%computed)
    return computed

def metadynamics_search(xtb_driver,
                        workdir,
                        mtd_indices,
                        constraints,
                        parameters,
                        verbose=True,
                        nthreads=1):

    if verbose:
        print("-----------------------------------------------------------------")
        print("Performing metadynamics jobs 👷...")
        print(mtd_indices)
        print("with %i threads. Working..." % nthreads)

        
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        futures = []

        for mtd_index in mtd_indices:
            crest_jobs = react_utils.crest_jobs(
                xtb_driver, mtd_index,
                workdir +"/init", workdir + "/metadyn",
                constraints, parameters)

            for ji,j in enumerate(crest_jobs):
                futures += [pool.submit(j)]
                
                if verbose:
                    futures[-1].add_done_callback(
                        lambda _: print("🔨",end="",flush=True))

        for f in futures:
            f.result()
                    
    if verbose:
        print("\nDone!\n")

def quick_opt_job(xtb, xyz, level, xcontrol):
    # TODO comment
    with tempfile.NamedTemporaryFile(suffix=".xyz",
                                     dir=xtb.scratchdir) as T:
        T.write(bytes(xyz, 'ascii'))
        T.flush()
        opt = xtb.optimize(T.name,
                           T.name,
                           level=level,
                           xcontrol=xcontrol)
        opt()
        xyz, E = read_trajectory(T.name)
    return xyz, E

def metadynamics_refine(xtb_driver,
                        workdir,
                        reference,
                        mtd_indices,
                        constraints,
                        parameters,
                        verbose=True,
                        nthreads=1):
    refined_dir = workdir + "/CRE"
    mtd_dir = workdir + "/metadyn"
    os.makedirs(refined_dir, exist_ok=True)
    
    for mtd_index in mtd_indices:
        structures = []
        Es = []
        for file in glob(mtd_dir + "/mtd%4.4i_*.xyz"%mtd_index):
            structs, E = read_trajectory(file)
            structures += structs
            Es += E
            
        # Loosely optimize the structures in parallel
        if verbose:
            print("MTD%i> loaded %i structures, optimizing 📐..."
                  %(mtd_index, len(structures)))
        with ThreadPoolExecutor(max_workers=nthreads) as pool:
            futures = []
            for s in structures:
                future = pool.submit(
                    quick_opt_job, xtb_driver, s, parameters["optlevel"],
                    dict(wall=parameters["wall"],
                         constrain=constraints[mtd_index]))
                futures += [future]
                
        converged = []
        errors = []
        for f in futures:
            exc = f.exception()
            if exc:
                errors += [f]
            else:
                converged += [f.result()]

        if verbose:
            print("        converged 👍: %i"% len(converged))
            print("        errors 👎: %i"%len(errors))

        if verbose:
            print("        carefully selecting conformers 🔎...")

        fn = refined_dir + "/mtd%4.4i.xyz" % mtd_index
        f = open(fn, "w")
        for s, E in converged:
            f.write(s[0])
        f.close()

        # convert to kcal and use use 10 000 if its higher than 10 000
        Et = min(10000.0, parameters["ethreshold"] * hartree_ev * ev_kcalmol)
        # same
        Rt = min(10000.0, parameters["rthreshold"])
        
        cre = xtb_driver.cregen(reference,
                                fn, fn,
                                ethreshold=Et,
                                rthreshold=Rt)
        error = cre()
        s, E = read_trajectory(fn)
        if verbose:
            print("  → %i structures selected for reactions 🔥" % len(s))

def react(xtb_driver,
          workdir,
          mtd_indices,
          constraints,
          parameters,
          verbose=True,
          nthreads=1):
    if verbose:
        print("-----------------------------------------------------------------")
        print("Performing reactions ⏣ → 🔥...")
    
    # load all the structures
    meta = workdir+"/CRE"

    nreact = 0
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        futures = []

        for mtd_index in mtd_indices:
            structures, energies = read_trajectory(
                meta + "/mtd%4.4i.xyz" % mtd_index)

            if verbose:
                print("   MTD%i>\tn(⏣ 🡒 🔥) = %i"
                      %(mtd_index, len(structures)))
                
            for s in structures:
                futures += [pool.submit(
                    react_utils.reaction_job(
                        xtb_driver,
                        s,
                        mtd_index,
                        workdir + "/react%5.5i/" % nreact,
                        constraints,
                        parameters))]
                nreact = nreact + 1

    if verbose:
        print("No more work to do. 🛌📺\n\n")
