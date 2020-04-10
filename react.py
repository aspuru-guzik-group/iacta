from concurrent.futures import ThreadPoolExecutor
import react_utils
from io_utils import traj2str
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
                       optlevel="tight",
                       ethreshold=inf,
                       rthreshold=0.5,
                       emax=inf,
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
    parameters["emax"] = emax

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
    nmtd = Natoms               # TODO parameterize
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

    react_utils.dump_succ_opt(outputdir,
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
            structs, E = traj2str(file)
            structures += structs
            Es += E
            
        # Loosely optimize the structures in parallel
        if verbose:
            print("MTD%i>\tloaded %i structures, optimizing 📐..."
                  %(mtd_index, len(structures)))
        with ThreadPoolExecutor(max_workers=nthreads) as pool:
            futures = []
            for s in structures:
                future = pool.submit(
                    react_utils.quick_opt_job,
                    xtb_driver, s, parameters["optlevel"],
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
            f.write(s)
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
        s, E = traj2str(fn)
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
        print("\n")
        print("-----------------------------------------------------------------")
        print("  miniGabe 🧔,")
        print("           a reaction search ⏣ → 🔥 algorithm")
        print("-----------------------------------------------------------------")
    
    # load all the structures
    meta = workdir+"/CRE"

    # Make a large list of all the reactions
    if verbose:
        print("Loading metadynamics structures...")

    all_structures = []
    for mtd_index in mtd_indices:
        structures, energies = traj2str(
            meta + "/mtd%4.4i.xyz" % mtd_index)
        # CREGEN has sorted these too so the first element is the lowest
        # energy one.
        
        if verbose:
            print("   MTD%i> N(⏣ 🡒 🔥) = %i"
                  %(mtd_index, len(structures)))
            
        all_structures += [[mtd_index,structures]]
    if verbose:
        print("building round-robin worklist...")

    worklist = []
    while True:
        for mtd_index, ls_structs in all_structures:
            if len(ls_structs):
                worklist += [(mtd_index,ls_structs.pop(0))]
            else:
                all_structures.remove([mtd_index,ls_structs])
                
        if len(all_structures) == 0:
            break
        
    if verbose:
        print("📜 = %i reactions to compute"% len(worklist))
        print("    with the help of 🧔 × %i threads" % nthreads)

    nreact = 0
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        futures = []

        for mtd_index, structure in worklist:
            futures += [pool.submit(
                react_utils.reaction_job(
                    xtb_driver,
                    structure,
                    mtd_index,
                    workdir + "/react%5.5i/" % nreact,
                    constraints,
                    parameters))]
            nreact = nreact + 1

    if verbose:
        print("No more work to do! 🧔🍻\n\n")
