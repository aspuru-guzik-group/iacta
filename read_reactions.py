import pybel
import numpy as np
import glob
from constants import hartree_ev, ev_kcalmol
import io_utils
import os

def xyz2smiles(xyzfile, chiral=False):
    output = []
    if chiral:
        flags = {"c":1,"n":1}
    else:
        flags = {"c":1,"n":1, "i":1}
    for m in pybel.readfile("xyz",xyzfile):
        output+= [m.write(format="smi", opt=flags).rstrip()]
    return output

def read_reaction(react_folder,chirality=False):
    """Extract chemical quantities from a reaction.

    Parameters:
    -----------
    react_folder (str) : folder storing the reaction data, obtained from the
    reaction_job() routine in react.py.

    chirality (bool) : whether to separate enantiomers.

    Returns:
    -------
    TODO
    """
    opt = react_folder + "/opt.xyz"
    smiles = xyz2smiles(opt, chirality)
    E = np.loadtxt(react_folder + "/Eopt", ndmin=1)
    
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
        imin = np.argmin(E[start:end]) + start
        stable = True
        # Check if imin is a local minima
        if imin < len(smiles)-1:
            if E[imin]>E[imin+1]:
                stable = False
        if imin > 0:
            if E[imin]>E[imin-1]:
                stable = False

        if stable:
            imins += [imin]
                
    # Important potential energy surface points
    ipots = [imins[0]]
    # determine if these are stable, ts or errors
    stable = [1]
    for k in range(1,len(imins)):
        # maximum between minima
        imax = np.argmax(E[imins[k-1]:imins[k]]) + imins[k-1]
        # if imax is different from both, we add it to the pot curve too
        if imax != imins[k] and imax != imins[k-1]:
            ipots += [imax]
            stable += [0]
        ipots += [imins[k]]
        stable += [1]

    energies = np.array([E[k] for k in ipots])
    smiles = [smiles[k] for k in ipots]
    structures = [(opt,k) for k in ipots]
    out = {"E":energies,
           "SMILES":smiles,
           "files":structures,
           "is_stable":stable,
           "stretch_points":ipots}
    return out
        
def read_all_reactions(output_folder, verbose=True, chirality=False):
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
            read_out = read_reaction(f,chirality)
        except OSError:
            # Convergence failed
            failed += [f]
        else:
            if len(read_out["SMILES"])>1:
                pathways += [(read_out["E"],
                              read_out["SMILES"],
                              read_out["files"],
                              read_out["is_stable"],
                              "react"+f[-5:])]
            else:
                # No reaction in this pathway
                noreact += [f]

            
    if verbose:
        print(" - %6i that did not converge" % len(failed))
        print(" - %6i with no reactions" % len(noreact))
        print("--------------")
        print(" = %6i reactions found" % len(pathways))
    return pathways

def prune_pathways(pathways):
    reactions = {}

    for energies, species, files, stable, where in pathways:
        stables = []
        barriers = []
        current = 0.0
        for i in range(len(species)):
            if stable[i]:
                stables += [species[i]]
                barriers += [current]
                estable = energies[i]
                current = 0.0
            else:
                barrier = energies[i] - estable
                if barrier > current:
                    current = barrier
        barriers = np.array(barriers[1:])
        stables = tuple(stables)

        if stables in reactions:
            if (np.max(barriers)< np.max(reactions[stables][0])):
                # if all the barriers are smaller, we make this the new best
                # reaction.
                reactions[stables] = (barriers, energies, species,
                                      files, stable, where)
        else:
            reactions[stables] = (barriers, energies, species,
                                  files, stable, where)

    return reactions



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Read reaction products."
        )
    parser.add_argument("folder", help="Folder containing the react*** files.")
    parser.add_argument("-o", help="Folder where final results are stored.",
                        default=None)
    parser.add_argument("-c", help="Resolve chiral compounds. Defaults to no.",
                        action="store_true")
    args = parser.parse_args()
    folder =args.folder
    outfolder = args.o
    if outfolder is None:
        outfolder = folder + "/results"
    os.makedirs(outfolder, exist_ok=True)
    
    pathways = read_all_reactions(folder, chirality=args.c)
    reactions = prune_pathways(pathways)
    print("   %6i reaction pathways are unique." % len(reactions))
    print("\n\n")
    mol_index = {}
    index = 1

    # Sort the reaction by barrier
    reactions = sorted(reactions.items(), key=lambda x: x[1][0].max())
    
    for key,val in reactions:
        for molecule in key:
            if molecule not in mol_index.keys():
                mol_index[molecule] = index
                mol = pybel.readstring("smi", molecule)
                # save to file
                mol.draw(show=False,
                         filename=outfolder + "/mol%3.3i.png"%index)
                index += 1


    start = "|--=====-----==-- REACTIONS AND BARRIERS (kcal/mol) --==------====---|\n"
    start+= "  #     File         Reaction\n"
    index = 1
    f = open(outfolder + "/summary", "w")
    print(start, end="")
    f.write(start)
    for key, val in reactions:
        upper = "%3i >   " % index + val[5]+"   "+  key[0]
        for i in range(1,len(key)):
            new = " == %.2f ==> " % (val[0][i-1]*hartree_ev * ev_kcalmol) + key[i]
            upper += new
        upper += "\n"

        print(upper, end="")
        f.write(upper)

        with open(outfolder + "/reactions%i.xyz" % index, "w") as ftraj:
            for file, ind in val[3]:
                xyz = io_utils.read_trajectory(file, ind)
                ftraj.write(xyz)
        index += 1

    table = ("\n\nTable of molecules\n"
              +  "------------------\n")
    for smile, index in mol_index.items():
        table += "mol%3.3i.png\t"% index + smile + "\n"

    print(table)
    f.write(table)
        
    f.close()
    

        
    
                    
        

        

            

