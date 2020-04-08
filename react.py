from concurrent.futures import ThreadPoolExecutor
import react_utils
from io_utils import read_trajectory, dump_succ_opt
from math import inf
import os

"""
This file contains a bunch of user-friendly, multithreaded drivers for
react_utils.
"""

def default_parameters(Natoms,
                       nmtd=80,
                       optlevel="tight",
                       ethreshold=inf,
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
    parameters["nmtd"] = nmtd
    parameters["optlevel"] = optlevel
    parameters["ethreshold"] = ethreshold

    # Logging
    if log_level > 0:
        parameters["log_opt_steps"] = True
    else:
        parameters["log_opt_steps"] = False

    
    # xTB additional parameters
    parameters["wall"] = ("potential=logfermi",
                          "sphere: auto, all")

    # Metadynamics parameters (from CREST)
    total_time = 0.5 * Natoms
    parameters["CREST"] = [
        ("save=10", "kpush=0.033", "alpha=1.300"),
        ("save=10", "kpush=0.017", "alpha=1.300"),
        ("save=10", "kpush=0.008", "alpha=1.300"),
        ("save=10", "kpush=0.033", "alpha=0.780"),
        ("save=10", "kpush=0.017", "alpha=0.780"),
        ("save=10", "kpush=0.008", "alpha=0.780"),
        ("save=10", "kpush=0.033", "alpha=0.468"),
        ("save=10", "kpush=0.017", "alpha=0.468"),
        ("save=10", "kpush=0.008", "alpha=0.468"),
        ("save=10", "kpush=0.033", "alpha=0.281"),
        ("save=10", "kpush=0.017", "alpha=0.281"),
        ("save=10", "kpush=0.008", "alpha=0.281"),
        ("save=10", "kpush=0.011", "alpha=0.100"),
        ("save=10", "kpush=0.055", "alpha=0.800")]

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
        print("Performing initial stretching...")

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
        print("Performing metadynamics job on indices...")
        print(mtd_indices)
        print("with %i threads." % nthreads)

    futures = []
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        for mtd_index in mtd_indices:
            crest_jobs = react_utils.crest_jobs(
                xtb_driver, mtd_index,
                workdir +"/init", workdir + "/metadyn",
                constraints, parameters)

            if verbose:
                print("   index %3i : submitting %i jobs..."
                      %(mtd_index, len(crest_jobs)))
                
            for ji,j in enumerate(crest_jobs):
                future=pool.submit(j)
                
                if verbose:
                    future.add_done_callback(
                        lambda _: print("X",end="",flush=True))
                    
    if verbose:
        print("\nDone!\n")

def react(xtb_driver,
          workdir,
          mtd_indices,
          constraints,
          parameters,
          verbose=True,
          nthreads=1):
    if verbose:
        print("-----------------------------------------------------------------")
        print("Performing reactions...")
    
    # load all the structures
    meta = workdir+"/metadyn"

    nreact = 0
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        futures = []

        for mtd_index in mtd_indices:
            structures, energies = read_trajectory(
                meta + "/mtd%4.4i.xyz" % mtd_index)

            if verbose:
                print("   mtdi = %4i    n(react) = %i"
                      %(mtd_index+1, len(structures)))
                
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
        print("Done!\n\n")
