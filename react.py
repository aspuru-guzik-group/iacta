from concurrent.futures import ThreadPoolExecutor
import react_utils
from io_utils import traj2str, read_xtb_hessian
import numpy as np
from math import inf
import os
from glob import glob
import tempfile
from constants import hartree_ev, ev_kcalmol

"""
This file contains a bunch of user-friendly, multithreaded drivers for
react_utils.
"""

def generate_initial_structures(xtb_driver,
                                workdir,
                                guess_xyz_file,
                                atoms,
                                low,high,npts,
                                parameters,
                                verbose=True):


    if verbose:
        print("-----------------------------------------------------------------")
        print("Performing initial driving 💪😎💪...")

    outputdir = workdir + "/init"
    os.makedirs(outputdir)

    structures, energies = react_utils.stretch(
        xtb_driver, guess_xyz_file,
        atoms,
        low, high, npts,
        parameters,
        failout=outputdir + "/FAILED",
        verbose=True)

    react_utils.dump_succ_opt(outputdir,
                              structures,
                              energies,
                              split=True)
    if verbose:
        print("Done!")

def metadynamics_search(xtb_driver,
                        workdir,
                        mtd_indices,
                        atoms, low, high, npts,
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
            mtd_jobs = react_utils.metadynamics_jobs(
                xtb_driver, mtd_index,
                atoms, low, high, npts,
                workdir +"/init", workdir + "/metadyn", parameters)

            for ji,j in enumerate(mtd_jobs):
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
                        atoms, low, high, npts,
                        parameters,
                        verbose=True,
                        nthreads=1):
    refined_dir = workdir + "/CRE"
    mtd_dir = workdir + "/metadyn"
    os.makedirs(refined_dir, exist_ok=True)

    # points
    points = np.linspace(low, high, npts)

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
                    xtb_driver, s, parameters["optcregen"],
                    dict(wall=parameters["wall"],
                         cma="",
                         constrain = react_utils.make_constraint(
                             atoms,
                             points[mtd_index], parameters["force"])
                         )
                    )
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

        cre = xtb_driver.cregen(reference,
                                fn, fn,
                                ewin=parameters["ewin"],
                                rthr=parameters["rthr"],
                                ethr=parameters["ethr"],
                                bthr=parameters["bthr"])

        error = cre()
        s, E = traj2str(fn)
        if verbose:
            print("  → %i structures selected for reactions 🔥" % len(s))

def react(xtb_driver,
          workdir,
          mtd_indices,
          atoms, low, high, npts,
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
    os.makedirs(workdir + "/reactions/")
    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        futures = []

        for mtd_index, structure in worklist:
            futures += [pool.submit(
                react_utils.reaction_job(
                    xtb_driver,
                    structure,
                    mtd_index,
                    atoms, low, high, npts,
                    workdir + "/reactions/%5.5i/" % nreact,
                    parameters))]
            nreact = nreact + 1

        for f in futures:
            # crash if f raised exception
            f.result()

    if verbose:
        print("No more work to do! 🧔🍻\n\n")
