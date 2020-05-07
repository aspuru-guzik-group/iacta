import pybel
import numpy as np
import glob
from constants import hartree_ev, ev_kcalmol
import io_utils
import analysis
import pandas as pd
import os
import importlib

importlib.reload(analysis)

def read_all_reactions(output_folder,
                       verbose=True,
                       restart=True,
                       save=True,
                       resolve_chiral=False):
    """Read and parse all reactions in a given folder."""
    folders = glob.glob(output_folder + "/react[0-9]*")
    if verbose:
        print("Parsing folder <%s>, with" % output_folder)
        print("   %6i trajectories..." % len(folders))

    failed = []
    noreact = []
    pathways = []
    if restart:
        try:
            old_df = pd.read_pickle(output_folder+"/results_raw.pkl")
        except FileNotFoundError:
            old_df = pd.DataFrame()
        else:
            if verbose:
                print(" - %6i trajectories in restart file" % len(old_df))

    old_indices = old_df.index
    new_indices = []
    for f in folders:
        # nasty parsing...
        index = int(f[-9:-4])
        if index in old_indices:
            # already in restart file
            continue
        
        imtd = int(f[-3:])
        try:
            read_out = analysis.read_reaction(f,
                                              resolve_chiral=resolve_chiral)
        except OSError:
            # Convergence failed
            failed += [f]
        else:
            read_out['folder'] = "react"+f[-9:]
            read_out['iMTD'] = imtd
            new_indices += [index]
            pathways += [read_out]

            
    if verbose:
        print(" - %6i that did not converge" % len(failed))
        print("--------------")

    new = pd.DataFrame(pathways, index=new_indices)
    data = old_df.append(new)

    if verbose:
        if len(old_df):
            print(" = %6i new pathways loaded" % len(pathways))
            print("   %6i pathways including restart" % len(data))
        else:
            print(" = %6i pathways" % len(pathways))
            
    if save:
        if len(new):
            if verbose:
                print("       ... saving new data ...")
            data.to_pickle(output_folder+"/results_raw.pkl")

    if verbose:
        print("                                      done.")
    
    return data

def get_species_table(pathways, verbose=True):
    if verbose:
        print("\nBuilding table of chemical species, from")
        print("   %6i reaction pathways" % len(pathways))

    species = {}
    for irow, row in pathways.iterrows():
        for k in range(len(row.is_stable)):
            if row.is_stable[k]:
                smi = row.SMILES[k]
                if smi in species:
                    if row.E[k] < species[smi]['E']:
                        species[smi] = {
                            'smiles':smi,
                            'E':row.E[k],
                            'folder':row.folder,
                            'position':row.stretch_points[k]}
                else:
                    species[smi] = {
                        'smiles':smi,
                        'E':row.E[k],
                        'folder':row.folder,
                        'position':row.stretch_points[k]}

    k,v = zip(*species.items())
    out = pd.DataFrame(v).sort_values('E').reset_index(drop=True)
    if verbose:
        print("   %6i stable-ish species found" % len(out))
        evspan = (out.E.max() - out.E.min()) * hartree_ev
        print("          with energies spanning %3.1f eV" % (evspan))
        print("          saving structures...")
        print("                                      done.")

    return out

def reaction_network_layer(pathways, reactant):
    to_smiles = []
    ts_i = []
    ts_E = []
    folder = []
    imtd = []
    for k, rowk in pathways.iterrows():
        if reactant in rowk.SMILES:
            i = rowk.SMILES.index(reactant)
            E = rowk.E
            stable = rowk.is_stable

            # forward
            for j in range(i+1, len(stable)):
                if stable[j]:
                    # we have a stable -> stable reaction
                    to_smiles += [rowk.SMILES[j]]
                    tspos = np.argmax(E[i:j]) + i
                    ts_i += [rowk.stretch_points[tspos]]
                    ts_E += [E[tspos]]
                    folder += [rowk.folder]
                    imtd += [rowk.iMTD]

            # backward
            for j in range(i-1, -1, -1):
                if stable[j]:
                    # we have a stable -> stable reaction
                    to_smiles += [rowk.SMILES[j]]
                    tspos = np.argmax(E[j:i]) + j
                    ts_i += [rowk.stretch_points[tspos]]
                    ts_E += [E[tspos]]
                    folder += [rowk.folder]
                    imtd += [rowk.iMTD]                    

    out = pd.DataFrame({
        'from':[reactant] * len(to_smiles),
        'to':to_smiles,
        'E_TS':ts_E,
        'i_TS':ts_i,
        'folder':folder,
        'iMTD':imtd})

    return out


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

    pathways = read_all_reactions(folder, resolve_chiral=args.c)
    species = get_species_table(pathways)


    reactant, E = io_utils.traj2smiles(folder + "/init_opt.xyz", index=0)
    verbose=True
    if verbose:
        print("\nBuilding reaction network iteratively,")
        print("Reactant: %s" % reactant)

    layer = reaction_network_layer(pathways, reactant)
    E0 = species[species.smiles == reactant].E
    if verbose:
        products = set(layer.to)
        print("   %6i reactions found, with" % len(layer))
        print("   %6i products:" %len(products))
        for p in products:
            print("    %s" % p)
            # print("    dH:%8.4f kcal/mol" % 

    
    



    

    """
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


    start = "|--=====-----==--  ⏣ ⮂ 🔥 REACTIONS AND ⛰  BARRIERS (kcal/mol) --==------====---|\n"
    start+= "  #     File         Reaction\n"
    index = 1
    f = open(outfolder + "/summary", "w")
    print(start, end="")
    f.write(start)
    for key, val in reactions:
        upper = "%3i >   " % index + val[5]+"   "+  key[0]
        for i in range(1,len(key)):
            new = " = ⛰ %.2f => " % (val[0][i-1]*hartree_ev * ev_kcalmol) + key[i]
            upper += new
        upper += "\n"

        print(upper, end="")
        f.write(upper)

        with open(outfolder + "/reactions%i.xyz" % index, "w") as ftraj:
            for file, ind in val[3]:
                xyz, E = io_utils.traj2str(file, ind)
                ftraj.write(xyz)
        index += 1

    table = ("\n\nTable of molecules\n"
              +  "------------------\n")
    for smile, index in mol_index.items():
        table += "📁 mol%3.3i.png\t"% index + smile + "\n"

    print(table)
    f.write(table)
        
    f.close()
    """    

        
    
                    
        

        

            

