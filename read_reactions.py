import pybel
import numpy as np

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
    mols (list of str) : smiles for each of the reactant, product and
    intermediate species in the reaction.

    energies (ndarray) : energy values for each of the species.

    barriers (ndarray) : barrier energy between each of the species, defined
    as the maximum energy reached between species.
    """


    smiles = xyz2smiles(react_folder + "/opt.xyz")
    E = np.loadtxt(react_folder + "/Eopt")
    
    mols = [smiles[0]]
    energies = []
    barriers = []
    last = 0

    # loop through smiles and detect changes
    for si,s in enumerate(smiles):
        if s == mols[-1]:
            pass
        else:
            # we have a state change at index=si
            # minimum energy point for the current state
            minE_point = np.argmin(E[last:si]) + last
            energies += [E[minE_point]]
            # Energy barrier height
            # maximum of energy between the minimum and the state change
            barriers += [np.max(E[minE_point:si])]

            mols += [s]           # save the new state
            last = si
    energies += [np.min(E[last:])]
    return mols, energies, barriers

react_folder = "./output/react00000"
mols,energies,barriers = analyze_reaction(react_folder)
