import os
import shutil
import subprocess
import numpy as np
import tempfile
import re
from analysis import read_reaction
from io_utils import traj2str
import json

# ------------------- utility routines -----------------------------------#
def dump_succ_opt(output_folder, structures, energies,
                  split=False):

    os.makedirs(output_folder, exist_ok=True)
    # Dump the optimized structures in one file
    with open(output_folder + "/opt.xyz", "w") as f:
        for s in structures:
            f.write(s)

    if split:
        # Also dump the optimized structures in many files
        for stepi, s in enumerate(structures):
            with open(output_folder + "/opt%4.4i.xyz" % stepi, "w") as f:
                f.write(s)

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
        xyz, E = traj2str(T.name, 0)
    return xyz, E

def stretch_constraint(atom1, atom2, val, force):
    return ("force constant=%f" % force,
            "distance: %i, %i, %f" % (atom1, atom2, val))

# ------------------------------------------------------------------------------#
def stretch(xtb, initial_xyz,
            atom1, atom2, pts,
            parameters,
            failout=None,
            intermediate_steps=False,
            verbose=True):
    """ TODO
    """
    # Make scratch files
    fdc, current = tempfile.mkstemp(suffix=".xyz", dir=xtb.scratchdir)
    fdl, log = tempfile.mkstemp(suffix="_log.xyz", dir=xtb.scratchdir)
    fdr, restart = tempfile.mkstemp(dir=xtb.scratchdir)

    # prepare the current file
    shutil.copyfile(initial_xyz, current)
    structures = []
    energies = []

    if verbose:
        print("  %i successive optimizations" % len(pts))

    for i in range(len(pts)):
        direction = "->"
        if verbose:
            print("      " + direction + "  %3i " % i, end="")

        opt = xtb.optimize(current,
                           current,
                           log=log,
                           failout=failout,
                           level=parameters["optim"],
                           restart=restart,
                           xcontrol=dict(
                               constrain=stretch_constraint(
                                   atom1, atom2,
                                   pts[i], parameters["force"]),
                               wall=parameters["wall"]))

        error = opt()
        if error != 0:
            # An optimization failed, we get out of this loop.
            break

        news, newe = traj2str(log)
        nsteps = len(newe)
        minsteps = parameters["steps_before_save"]
        if intermediate_steps:
            if nsteps > minsteps:
                structures += news[minsteps:]
                energies += newe[minsteps:]
            else:
                structures += [news[-1]]
                energies += [newe[-1]]

        else:
            structures += [news[-1]]
            energies += [newe[-1]]

        if verbose:
            print("   step👣=%4i    energy💡= %9.5f Eₕ"%(nsteps,
                                                         energies[-1]))

        # todo put back
        # if newe[-1] > parameters["emax"] and barrier:
        #     if verbose:
        #         print("   ----- energy threshold exceeded -----")
        #     break

    # Got to make sure that you close the filehandles!
    os.close(fdc)
    os.close(fdl)
    os.close(fdr)
    os.remove(current)
    os.remove(log)
    os.remove(restart)
    return structures, energies

def stretch_xtbscan(xtb, initial_xyz,
            atom1, atom2, low, high, npts,
            parameters,
            failout=None,
            verbose=True):
    """Optimize a structure through successive constraints.

    TODO FIX

    Parameters:
    -----------

    xtb (xtb_driver) : driver object for xtb.

    initial_xyz (str): path to initial structure xyz file.

    parameters (dict) : additional parameters, as obtained from
    default_parameters() above. TODO: Describe parameters in more details

    Optional Parameters:
    --------------------

    failout (str) : Path where logging information is output if the
    optimization fails.

    verbose (bool) : print information about the run. defaults to True.

    Returns:
    --------

    structures : list of .xyz formatted strings that include every structure
    through the multiple optimization.

    energies : list of floats of xtb energies (in Hartrees) for the structures.

    """
    # Make scratch files
    fdc, current = tempfile.mkstemp(suffix=".xyz", dir=xtb.scratchdir)

    # prepare the current file
    shutil.copyfile(initial_xyz, current)

    opt = xtb.optimize(current,
                       current,
                       failout=failout,
                       level=parameters["optim"],
                       xcontrol=dict(
                           wall=parameters["wall"],
                           constrain=stretch_constraint(atom1, atom2,
                                                        low, parameters["force"]),
                           scan=("1: %f, %f, %i" % (low, high, npts),)))

    error = opt()
    structs, energies = traj2str(current)

    if verbose:
        for k, E in enumerate(energies):
            print("   👣=%4i    energy💡= %9.5f Eₕ"%(k, E))

    # Got to make sure that you close the filehandles!
    os.close(fdc)
    os.remove(current)
    return structs, energies

def metadynamics_jobs(xtb,
                      mtd_index,
                      atom1, atom2, pts,
                      input_folder,
                      output_folder,
                      parameters):
    """Return a metadynamics search job for other "transition" conformers.

    mtd_index is the index of the starting structure to use as a starting
    point in the metadynamics run. Returns an unevaluated xtb_job to be used
    in ThreadPool.

    Parameters:
    -----------

    xtb (xtb_driver) : driver object for xtb.

    mtd_index (int) : index of the starting structure, as generated from
    generate_starting_structures. Correspond to the index of the element in
    the constratins list for which metadynamics is run.

    input_folder (str) : folder containing the MTD input structures.

    output_folder (str) : folder where results are stored.

    constraints (list): list of constraints. Each constraints should be a
    tuple of parameters to be put into a $constrain group in xcontrol.

    parameters (dict) : additional parameters, as obtained from
    default_parameters() above. TODO: Describe parameters in more details

    Optional Parameters:
    --------------------

    verbose (bool) : print information about the run. defaults to True.

    Returns:
    --------

    None

    """
    os.makedirs(output_folder, exist_ok=True)

    mjobs = []
    meta = parameters["metadynamics"]
    inp = input_folder + "/opt%4.4i.xyz" % mtd_index
    # Set the time of the propagation based on the number of atoms.
    with open(inp, "r") as f:
        Natoms = int(f.readline())
    md = meta["md"] + ["time=%f" % (meta["time_per_atom"] * Natoms),
                       "dump=%f" % (meta["time_per_atom"] * Natoms
                                    * 1000.0/meta["nmtd"])]
    S = "save=%i" % meta["save"]

    for metadyn_job, metadyn_params in enumerate(meta["jobs"]):
        outp = output_folder + "/mtd%4.4i_%2.2i.xyz" % (mtd_index,metadyn_job)
        mjobs += [
            xtb.metadyn(inp, outp,
                        failout=output_folder +
                        "/FAIL%4.4i_%2.2i.xyz" % (mtd_index, metadyn_job),
                        xcontrol=dict(
                            wall=parameters["wall"],
                            metadyn=metadyn_params + [S],
                            md=md,
                            cma="",
                            constrain=stretch_constraint(atom1, atom2,
                                                         pts[mtd_index],
                                                         parameters["mtd_force"])))]
    return mjobs


def reaction_job(xtb,
                 initial_xyz,
                 mtd_index,
                 atom1, atom2, points,
                 output_folder,
                 parameters):
    """Take structures generated by metadynamics and build a reaction trajectory.

    Takes a structure generated by metadynamics_search() and successively
    optimize it by applying the sequence of $constrain objects in constraints
    forward from mtd_index to obtain products, and backward from mtd_index to
    obtain reactants.

    This is the final step in the reaction space search, and it generates
    molecular trajectories. It should be noted that this function returns an
    unevaluated job, to be fed to ThreadPool.

    Parameters:
    -----------

    xtb (xtb_driver) : driver object for xtb.

    initial_xyz (str) : initial xyz as a string.

    mtd_index (int) : Index of the element in the constraints list that
    generated the starting structures.

    TODO

    output_folder (str) : folder where results are stored.

    parameters (dict) : additional parameters, as obtained from
    default_parameters() above. TODO: Describe parameters in more details

    Optional Parameters:
    --------------------

    verbose (bool) : print information about the run. defaults to True.

    Returns:
    --------

    react_job() : function which, when evaluated, computes the trajectory

    """

    def react_job():
        os.makedirs(output_folder, exist_ok=True)

        with open(output_folder + "/initial.xyz", "w") as f:
            f.write(initial_xyz)

        forw = points[mtd_index+1:]
        back = points[:mtd_index][::-1]
        start = points[mtd_index]
        # note: we want back to start at the same point as forward, otherwise
        # we get a lot more stretch on the backward trajectory and weird stuff
        # happens. In this way, we stretch to start on the first optimization,
        # then move onward from there forward, and backward from there too.
        # This avoids way overstretching the first step of the backward
        # propagation, or doing optimizations twice.

        # We want to make sure to optimize the initial xyz so that both
        # forward and backward start from optimized structures.
        opt = xtb.optimize(output_folder + "/initial.xyz",
                           output_folder + "/start.xyz",
                           failout=output_folder + "/FAILED_OPT",
                           level=parameters["optim"],
                           xcontrol=dict(
                               cma="",
                               wall=parameters["wall"],
                               constrain=stretch_constraint(
                                   atom1, atom2,
                                   start, parameters["force"])))
        opt()
        sstructs, se = traj2str(output_folder+"/start.xyz", index=0)

        # Forward reaction
        if len(forw)>1:
            fstructs, fe = stretch(
                xtb, output_folder + "/start.xyz",
                atom1, atom2, forw,
                parameters,
                failout=output_folder + "/FAILED_FORWARD",
                intermediate_steps=True,
                verbose=False)          # otherwise its way too verbose
        else:
            fstructs = []
            fe = []

        # Backward reaction
        if len(back)>1:
            bstructs, be = stretch(
                xtb, output_folder + "/start.xyz",
                atom1, atom2, back,
                parameters,
                failout=output_folder + "/FAILED_BACKWARD",
                intermediate_steps=True,
                verbose=False)          # otherwise its way too verbose
        else:
            bstructs = []
            be = []

        # Dump forward reaction and backward reaction quantities
        dump_succ_opt(output_folder,
                      bstructs[::-1] + [sstructs] + fstructs,
                      be[::-1] + [se] + fe,
                      split=False)

        # Now read the reaction we dumped
        isomeric = read_reaction(output_folder, resolve_chiral=True)
        nisomeric = read_reaction(output_folder, resolve_chiral=False)

        isomeric['mtdi'] = int(mtd_index)
        nisomeric['mtdi'] = int(mtd_index)

        # summarized
        with open(output_folder + "/react.json","w") as f:
            json.dump(nisomeric, f)
        with open(output_folder + "/react-iso.json","w") as f:
            json.dump(isomeric, f)


    return react_job
