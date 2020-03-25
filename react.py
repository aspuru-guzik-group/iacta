import os
import shutil
import subprocess
import numpy as np

def default_parameters(N,
                       nmtd=80,
                       optlevel="tight"):
    parameters = {}
    parameters["nmtd"] = nmtd
    parameters["optlevel"] = optlevel
    # xTB additional parameters
    parameters["wall"] = ("potential=logfermi",
                          "sphere: auto, all")

    # Metadynamics parameters (somewhat adapted from CREST)
    total_time = 0.5 * N
    dumpstep = 1000 * total_time/nmtd
    parameters["metadyn"] = ("save=100", 
                             "kpush=0.2",
                             "alp=0.8")
    parameters["md"] = ("shake=0",
                        "step=2",
                        "dump=%f"%dumpstep,
                        "time=%f" % total_time)

    return parameters

def generate_starting_structures(xtb,
                                 initial_xyz,
                                 output_folder,
                                 constraints,
                                 parameters,
                                 verbose=True):
    metapath = output_folder + "/metadyn"
    os.makedirs(metapath, exist_ok=True)
    shutil.copy(initial_xyz, metapath + "/current.xyz")

    if verbose:
        print("Stretching initial structure to generate guesses for metadyn...")
        print("---------------------------------------------------------------")

    for i in range(len(constraints)):
        if verbose:
            print("    %i out of %i" %(i+1, len(constraints)))
        opt = xtb.optimize(metapath + "/current.xyz",
                           metapath + "/current.xyz",
                           level=parameters["optlevel"],
                           xcontrol=dict(
                               wall=parameters["wall"],
                               constrain=constraints[i]))
        opt()
        # add to metadyn starting structures
        shutil.copy(metapath + "/current.xyz", metapath + "/in%5.5i.xyz" % i)

def metadynamics_search(xtb,
                        mtd_index,
                        output_folder,
                        constraints,
                        parameters,
                        verbose=True):
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
    
    if verbose:
        print("Performing reactions on constraint %i" % (mtd_index+1))

    react_path = output_folder + "/reactions"
    os.makedirs(react_path, exist_ok=True)

    shutil.copy(output_folder + "/metadyn/out%5.5i.xyz" % mtd_index,
                react_path + "/current.xyz")

    for j in range(mtd_index, len(constraints)):
        print("     ----> forward  %i out of %i" %(j+1, len(constraints)))

        opt = xtb.multi_optimize(react_path + "/current.xyz",
                                 react_path + "/current.xyz",
                                 level=parameters["optlevel"],
                                 xcontrol=dict(
                                     wall=parameters["wall"],
                                     constrain=constraints[j]))
        opt()
        shutil.copyfile(react_path + "/current.xyz",
                        react_path + "/prop%2.2i_%3.3i.xyz" % (mtd_index,j))

    print("     ----> forward to products")
    opt = xtb.multi_optimize(react_path + "/current.xyz",
                             react_path + "/products_%2.2i.xyz" %mtd_index,
                             level=parameters["optlevel"],
                             xcontrol=dict(wall=parameters["wall"]))
    opt()

    # copy starting point for backward reaction dynamics
    shutil.copyfile(react_path + "/prop%2.2i_%3.3i.xyz" % (mtd_index,mtd_index),
                    react_path + "/current.xyz")


    for j in range(mtd_index-1, -1, -1):
        print("     ----> backward %i out of %i" %(j+1, len(constraints)))

        opt = xtb.multi_optimize(react_path + "/current.xyz",
                                 react_path + "/current.xyz",
                                 level=parameters["optlevel"],
                                 xcontrol=dict(
                                     wall=parameters["wall"],
                                     constrain=constraints[j]))
        opt()
        shutil.copyfile(react_path + "/current.xyz",
                        react_path + "/prop%2.2i_%3.3i.xyz" % (mtd_index,j))


    print("     ----> backward to reactants")
    opt = xtb.multi_optimize(react_path + "/current.xyz",
                             react_path + "/reactants_%2.2i.xyz" % mtd_index,
                             level=parameters["optlevel"],
                             xcontrol=dict(wall=parameters["wall"]))
    opt()
