import pybel
import numpy as np
import glob
import constants
import io_utils

def xyz2smiles(xyzfile, chiral=False):
    output = []
    for m in pybel.readfile("xyz",xyzfile):
        output+= [m.write(format="smi", opt={"c":1,"n":1,
                                             "i":chiral*1}).rstrip()]
    if len(output) == 1:
        return output[0]
    else:
        return output

def read_reaction(react_folder):
    """Extract chemical quantities from a reaction.

    Parameters:
    -----------
    react_folder (str) : folder storing the reaction data, obtained from the
    reaction_job() routine in react.py.

    Returns:
    -------
    A populated ReactionPathway object
    """
    opt = react_folder + "/opt.xyz"
    smiles = xyz2smiles(opt)
    
    E = np.loadtxt(react_folder + "/Eopt")
    mols = [smiles[0]]
    regions = []
    rstart = 0

    # loop through smiles and detect changes in smiles
    for si,s in enumerate(smiles):
        if s == mols[-1]:
            pass
        else:
            mols += [s]
            regions += [(rstart,si)]
            rstart = si

    regions += [(rstart, len(smiles))]

    imins = []
    for start,end in regions:
        imins+= [np.argmin(E[start:end]) + start]

    # Important potential energy surface points
    ipots = [imins[0]]
    for k in range(1,len(imins)):
        # maximum between minima
        imax = np.argmax(E[imins[k-1]:imins[k]]) + imins[k-1]
        # if imax is different from both, we add it to the pot curve too
        if imax != imins[k] and imax != imins[k-1]:
            ipots += [imax]
        ipots += [imins[k]]

    energies = np.array([E[k] for k in ipots])
    smiles = [smiles[k] for k in ipots]
    structures = [(opt,k) for k in ipots]
    return energies, smiles, structures
        
def read_all_reactions(output_folder, verbose=True):
    """Read and parse all reactions in a given folder."""
    folders = glob.glob(output_folder + "/react[0-9]*")
    if verbose:
        print("Parsing folder <%s>, with" % output_folder)
        print("   %6i trajectories..." % len(folders))

    failed = []
    noreact = []
    pathways = []
    
    for f in folders:
        try:
            energies, smiles, structures = read_reaction(f)
        except OSError:
            # Convergence failed
            failed += [f]

        if len(smiles)>1:
            pathways += [(energies, smiles, structures)]
        else:
            # No reaction in this pathway
            noreact += [f]

            
    if verbose:
        print(" - %6i that did not converge" % len(failed))
        print(" - %6i with no reactions" % len(noreact))
        print("--------------")
        print(" = %6i reactions found" % len(pathways))
        # print("    with %6i distinct stable species" % len(species))
    return pathways


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Read reaction products."
        )

    parser.add_argument("folder", help="Folder containing the react*** files.")
    args = parser.parse_args()
    folder =args.folder
    pathways = read_all_reactions(folder)
