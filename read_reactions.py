import pybel
import numpy as np
import glob

class ReactionPathway:
    def __init__(self, species, minima, barriers, id):
        self.species = tuple(species)
        self.minima = minima
        self.barriers = barriers
        self.id = id

def xyz2smiles(xyzfile):
    output = []
    for m in pybel.readfile("xyz",xyzfile):
        output+= [m.write(format="smi", opt={"c":1,"n":1}).rstrip()]
    if len(output) == 1:
        return output[0]
    else:
        return output

def analyze_reaction(react_folder):
    """Extract chemical quantities from a reaction.

    Parameters:
    -----------
    react_folder (str) : folder storing the reaction data, obtained from the
    reaction_job() routine in react.py.

    Returns:
    -------
    A populated ReactionPathway object
    """
    try:
        smiles = xyz2smiles(react_folder + "/opt.xyz")
    except OSError:
        # run did not converge
        return None
    
    E = np.loadtxt(react_folder + "/Eopt")
    pathway = [smiles[0]]
    minima = []
    barriers = []
    last = 0

    # loop through smiles and detect changes
    for si,s in enumerate(smiles):
        if s == pathway[-1]:
            pass
        else:
            # we have a state change at index=si
            # minimum energy point for the current state
            minE_point = np.argmin(E[last:si]) + last
            minima += [E[minE_point]]
            # Energy barrier height
            # maximum of energy between the minimum and the state change
            barriers += [np.max(E[minE_point:si])]

            pathway += [s]           # save the new state
            last = si
    minima += [np.min(E[last:])]
    return ReactionPathway(pathway, minima, barriers, react_folder)


def read_all_reactions(output_folder, verbose=True):
    folders = glob.glob(output_folder + "/react[0-9]*")
    pathways = {}
    if verbose:
        print("Found %i reactions. Loading..." % len(folders), end="")

    failed = []
    for f in folders:
        rpath = analyze_reaction(f)
        if rpath:
            if rpath.species in pathways:
                pathways[rpath.species] += [rpath]
            else:
                pathways[rpath.species] = [rpath]
        else:
            failed += [f]
            
    if verbose:
        print("  done.")
        print("  successful: %i    failed: %i "
              % (len(folders)-len(failed), len(failed)))
        print("  unique reaction pathways: %i"% len(pathways))
        print("  chemical species: %i" % len(species))
    return pathways

pathways = read_all_reactions("output")

