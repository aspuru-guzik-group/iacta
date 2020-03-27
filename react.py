import os
import shutil
import subprocess
import numpy as np
import tempfile

def default_parameters(Natoms,
                       nmtd=80,
                       optlevel="tight"):
    """Generate a dictionary of default parameters for the reaction search.

    The parameters are loosely based on those described in the CREST paper.

    Parameters:
    -----------

    Natoms (int) : number of atoms.

    nmtd (int) : number of structures to extract from the metadynamics
    portion.

    Optional Parameters:
    --------------------

    optlevel (str) : Optimization level to used throughout. Defaults to
    "tight".

    Returns:
    --------

    dict : parameters dictionary to pass to other functions in react.py

    """
    parameters = {}
    parameters["nmtd"] = nmtd
    parameters["optlevel"] = optlevel
    # xTB additional parameters
    parameters["wall"] = ("potential=logfermi",
                          "sphere: auto, all")

    # Metadynamics parameters (somewhat adapted from CREST)
    total_time = 0.5 * Natoms
    dumpstep = 1000 * total_time/nmtd
    parameters["metadyn"] = ("save=100", 
                             "kpush=0.2",
                             "alp=0.8")
    parameters["md"] = ("shake=0",
                        "step=2",
                        "dump=%f"%dumpstep,
                        "time=%f" % total_time)

    return parameters

def successive_optimization(xtb,
                            initial_xyz,
                            constraints,
                            parameters,
                            verbose=True):
    """Optimize a structure through successive constraints.

    Takes an initial structure and successively optimize it applying the
    sequence of $constrain objects in constraints. These can be generated by
    code such as this one, which stretch a bond between atoms 10 and 11 in
    initial_xyz from 1.06 A to 3 x 1.06 A in 80 steps,

    stretch_factors = np.linspace(1.0, 3.0, 80)
    constraints = [("force constant = 0.5",
                    "distance: %i, %i, %f" % (10, 11, stretch * 1.06)
                   for stretch in stretch_factors]

    Parameters:
    -----------

    xtb (xtb_driver) : driver object for xtb.

    initial_xyz (str): path to initial structure xyz file.

    constraints (list): list of constraints. Each constraints should be a
    tuple of parameters to be put into a $constrain group in xcontrol.

    parameters (dict) : additional parameters, as obtained from
    default_parameters() above. TODO: Describe parameters in more details

    Optional Parameters:
    --------------------

    verbose (bool) : print information about the run. defaults to True.

    Returns:
    --------

    structures : list of .xyz formatted strings that include every structure
    through the multiple optimization.

    energies : list of floats of xtb energies (in Hartrees) for the structures.

    opt_indices : list of integers representing the indices of the optimized
    structures at each step.

    """

    # Make scratch files
    fdc, current = tempfile.mkstemp(suffix=".xyz", dir=xtb.scratchdir)
    fdl, log = tempfile.mkstemp(suffix="_log.xyz", dir=xtb.scratchdir)

    # prepare the current file
    shutil.copyfile(initial_xyz, current)
    structures = []
    energies = []
    opt_indices = []
    

    if verbose:
        print("  %i successive optimizations" % len(constraints))

    for i in range(len(constraints)):
        direction = "->"
        if verbose:
            print("      " + direction + "  %3i " % i, end="")
            
        opt = xtb.optimize(current,
                           current,
                           log=log,
                           level=parameters["optlevel"],
                           xcontrol=dict(
                               wall=parameters["wall"],
                               constrain=constraints[i]))
        opt()
        
        news, newe = read_trajectory(log)
        structures += news
        energies += newe
        opt_indices += [len(structures)-1]
        if verbose:
            print("   nsteps=%4i   Energy=%9.5f Eh"%(len(news), newe[-1]))

    os.remove(current)
    os.remove(log)
    return structures, energies, opt_indices
        

        

def metadynamics_search(xtb,
                        mtd_index,
                        output_folder,
                        constraints,
                        parameters,
                        verbose=True):
    """Perform a metadynamics search for other "transition" conformers.

    This function is to be called after generate_starting_structures().
    mtd_index is the index of the starting structure to use as a starting
    point in the metadynamics run.

    Parameters:
    -----------

    xtb (xtb_driver) : driver object for xtb.

    mtd_index (int) : index of the starting structure, as generated from
    generate_starting_structures. Correspond to the index of the element in
    the constratins list for which metadynamics is run.

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
    
    if verbose:
        print("Performing metadynamics on constraint %i" % (mtd_index+1))

    metapath = output_folder + "/metadyn"
    
    mjob = xtb.metadyn(metapath + "/in%5.5i.xyz" % mtd_index,
                       metapath + "/out%5.5i.xyz" % mtd_index,
                       xcontrol=dict(
                           wall=parameters["wall"],
                           metadyn=parameters["metadyn"],
                           md=parameters["md"],
                           constrain=constraints[mtd_index]))
    mjob()


def react(xtb,
          mtd_index,
          output_folder,
          constraints,
          parameters,
          verbose=True):
    """Take structures generated by metadynamics and build reaction trajectories.

    Takes an ensemble of structures generated by metadynamics_search() and
    successively optimize it by applying the sequence of $constrain objects in
    constraints forward from mtd_index to obtain products, and backward from
    mtd_index to obtain reactants. 

    This is the final step in the reaction space search, and it generates
    molecular trajectories.

    Parameters:
    -----------

    xtb (xtb_driver) : driver object for xtb.

    mtd_index (int) : index of the starting structure, as generated from
    generate_starting_structures. Correspond to the index of the element in
    the constratins list for which metadynamics is run.

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
    if verbose:
        print("Performing reactions on constraint %i" % (mtd_index+1))

    react_path = output_folder + "/reactions"
    trajectories = {}
    os.makedirs(react_path, exist_ok=True)

    # current file
    curr = react_path + "/current.xyz"
    shutil.copy(output_folder + "/metadyn/out%5.5i.xyz" % mtd_index, curr)

    for j in range(mtd_index, len(constraints)):
        if verbose:
            print("     ----> forward  %i out of %i" %(j+1, len(constraints)))

        opt = xtb.multi_optimize(curr,
                                 curr,
                                 level=parameters["optlevel"],
                                 xcontrol=dict(
                                     wall=parameters["wall"],
                                     constrain=constraints[j]))
        opt()
        read_ensemble_xyz(trajectories, curr)
        
        if j == mtd_index:
            # Copy starting point for backward trajectory
            shutil.copyfile(curr,
                            react_path + "/prop%2.2i.xyz" % (mtd_index))

    if verbose:
        print("     ----> forward to products")
        
    opt = xtb.multi_optimize(curr,
                             curr,
                             level=parameters["optlevel"],
                             xcontrol=dict(wall=parameters["wall"]))
    opt()
    read_ensemble_xyz(trajectories, curr)

    # copy starting point for backward
    shutil.copyfile(react_path + "/prop%2.2i.xyz" % (mtd_index),
                    curr)


    for j in range(mtd_index-1, -1, -1):
        if verbose:
            print("     ----> backward %i out of %i" %(j+1, len(constraints)))

        opt = xtb.multi_optimize(curr,
                                 curr,
                                 level=parameters["optlevel"],
                                 xcontrol=dict(
                                     wall=parameters["wall"],
                                     constrain=constraints[j]))
        opt()
        read_ensemble_xyz(trajectories, curr,
                          prepend=True)        
        
    if verbose:
        print("     ----> backward to reactants")
        
    opt = xtb.multi_optimize(curr,
                             curr,
                             level=parameters["optlevel"],
                             xcontrol=dict(wall=parameters["wall"]))
    opt()
    read_ensemble_xyz(trajectories, curr,
                      prepend=True)            
        
    return trajectories


def read_ensemble_xyz(trajs, filepath, prepend=False):
    """Read an ensemble xyz file and append or prepend to molecular trajectories."""
    
    file_in = open(filepath, 'r')
    mol_index = 0
    # Load molecules
    with open(filepath, 'r') as f:
        while True:
            first_line = file_in.readline()
            # EOF -> blank line
            if not first_line:
                break
                
            this_mol = first_line
            natoms = int(first_line.rstrip())
        
            for i in range(natoms+1):
                this_mol += file_in.readline()

            if mol_index in trajs:
                pass
            else:
                trajs[mol_index] = []
            
            if prepend:
                trajs[mol_index].insert(0, this_mol)
            else:
                trajs[mol_index].append(this_mol)
                
            mol_index += 1

    return trajs

def read_trajectory(filepath):
    """Read an xyz file containing a trajectory."""
    structures = []
    energies = []
    with open(filepath, 'r') as f:
        while True:
            first_line = f.readline()
            # EOF -> blank line
            if not first_line:
                break
                
            this_mol = first_line
            natoms = int(first_line.rstrip())

            comment_line = f.readline()
            this_mol += comment_line
            energies += [float(comment_line[8:25])]
        
            for i in range(natoms):
                this_mol += f.readline()

            structures += [this_mol]
    return structures,energies

def dump_succ_opt(output_folder, structures, energies, opt_indices, full=True):
    os.makedirs(output_folder, exist_ok=True)
    # Dump the optimized structures
    for stepi, oi in enumerate(opt_indices):
        with open(output_folder + "/opt%4.4i.xyz" % stepi, "w") as f:
            f.write(structures[oi])

    # Dump indices of optimized structures, energies of optimized structures,
    # all energies and all structures
    np.savetxt(output_folder + "/Eopt", np.array(energies)[opt_indices], fmt="%15.8f")

    if full:
        np.savetxt(output_folder + "/indices", opt_indices, fmt="%i")
        np.savetxt(output_folder + "/E", energies, fmt="%15.8f")
        with open(output_folder + "/log.xyz", "w") as f:
            for s in structures:
                f.write(s)





