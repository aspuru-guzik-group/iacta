from concurrent.futures import ThreadPoolExecutor
import react_utils
"""
This file contains a bunch of user-friendly, multithreaded drivers for
react_utils.
"""

def default_parameters(Natoms,
                       nmtd=80,
                       optlevel="tight",
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

    # Logging
    if log_level > 0:
        parameters["log_opt_steps"] = True
    else:
        parameters["log_opt_steps"] = False

    
    # xTB additional parameters
    parameters["wall"] = ("potential=logfermi",
                          "sphere: auto, all")

    # Metadynamics parameters (somewhat adapted from CREST)
    total_time = 0.5 * Natoms
    dumpstep = 1000 * total_time/nmtd
    parameters["metadyn"] = ("save=%i" % nmtd, 
                             "kpush=0.2",
                             "alp=0.8")

    if shake == 0:
        parameters["md"] = ("shake=0",
                            "step=2",
                            "dump=%f" % dumpstep,
                            "time=%f" % total_time)
    else:
        parameters["md"] = ("shake=%i" % shake,
                            "step=10",
                            "dump=%f" % dumpstep,
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

    structures, energies, opt_indices, opt_grads = react_utils.successive_optimization(
        xtb_driver, guess_xyz_file,
        constraints, parameters, verbose=verbose)

    react_utils.dump_succ_opt(workdir + "/init",
                              structures,energies,opt_indices,
                              extra=parameters["log_opt_steps"])
    if verbose:
        print("Done!\n")

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
    
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        for mtd_index in mtd_indices:
            pool.submit(
                react_utils.metadynamics_job(
                    xtb_driver, mtd_index,
                    workdir +"/init", workdir + "/metadyn",
                    constraints, parameters))

    if verbose:
        print("Done!\n")

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
        for mtd_index in mtd_indices:
            structures, energies = react_utils.read_trajectory(
                meta + "/mtd%4.4i.xyz" % mtd_index)

            for s in structures:
                rjob = react_utils.reaction_job(xtb_driver,
                                                s,
                                                mtd_index,
                                                workdir + "/react%5.5i/" % nreact,
                                                constraints,
                                                parameters)
                pool.submit(rjob)
                nreact = nreact + 1
                
        if verbose:
            print("  running %i jobs on %i threads" % (nreact, nthreads))

    if verbose:
        print("Done!\n\n")
