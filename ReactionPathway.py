import pybel
import numpy as np
import glob
import pandas as pd
from rmsd import rmsd

def xyz2smiles(atoms, xyz, chiral=False):
    string = "%i\n\n" % len(atoms)
    for i in range(len(atoms)):
        string += "%s\t%f\t%f\t%f\n" % (atoms[i], *xyz[i])

    mol = pybel.readstring("xyz",string)
    return mol.write(format="smi", opt={"c":1,"n":1,"i":chiral*1}).rstrip()


def split_output(folder, rmsd_threshold):
    output = np.load(folder + "/opt_raw.npz")    
    structures = output["structures"]
    x1 = structures[0]
    energies = output["energies"]
    atoms = output["atoms"]
    grads = abs(np.sum(output["gradients"], (1,2)))
    # Take the lowest energy structure as the reference (doesn't really matter
    # though I think)
    x0 = structures[0]

    rmsd_split = []
    start = 0
    rmsds = []
    # Split into structures based on rmsd
    for i in range(len(structures)):
        rmsds += [rmsd(x0,structures[i])]

    start = 0
    for i in range(1,len(structures)):
        if abs(rmsds[i] - rmsds[i-1]) > rmsd_threshold:
            rmsd_split += [(start,i)]
            start = i
    rmsd_split += [(start,len(structures))]

    # Minima of the gradients are those states we keep. These are either TS or
    # minima
    indices = [start + np.argmin(grads[start:end]) for start,end in rmsd_split]
    
    smiles = []
    E = []
    structs =[]
    for i in indices:
        smiles += [xyz2smiles(atoms,structures[i])]
        E += [energies[i]]
        structs +=[structures[i]]

    out = dict(cmpds=smiles,
               energies=E,
               structures=structs,
               folder=folder)
               

    return out


def read_all_pathways(output_folder, rmsd_threshold, verbose=True):
    """Read and parse all reactions in a given folder."""
    folders = glob.glob(output_folder + "/react[0-9]*")
    if verbose:
        print("Parsing folder <%s>, with" % output_folder)
        print("   %6i trajectories..." % len(folders))

    failed = []
    noreact = []
    pathways = []
    i = 0
    for f in folders:
        try:
            rpath = split_output(f, rmsd_threshold)
        except FileNotFoundError:
            # Convergence failed
            failed += [f]
        else:
            if len(rpath["cmpds"])>0:
                pathways += [rpath]
                i+= 1
            else:
                # No reaction in this pathway
                noreact += [f]
            
    if verbose:
        print(" - %6i that did not converge" % len(failed))
        print(" - %6i with only one species" % len(noreact))
        print("--------------")
        print(" = %6i multi-molecules pathways" % len(pathways))

    return pathways

# # kbT / hbar = 1/beta hbar = 1/beta (Eh^-1) * 1/hbar(eV fs) * 27 eV / hartree
# # EYRING_PREF = 1/constants.hbar * constants.hartree_ev # fs^-1
# def build_rate_matrix(pathways, species, beta):
#     """Build the matrix of rate from pathways."""
#     N = len(species)
#     R = np.zeros((N,N))
#     Z = np.zeros((N,N))             # number of pathways contributing to i->j
#     pref = EYRING_PREF / beta
#     for p in pathways:
#         for i in range(len(p.species)-1):
#             # species index into R
#             ri = species[p.species[i]]
#             rj = species[p.species[i+1]]
            
#             # minima
#             Ei = p.minima[i]
#             Ej = p.minima[i+1]
            
#             # barrier
#             barrier = p.barriers[i]

#             if ri == rj:
#                 # This is a pathway between two conformers or a reaction that
#                 # went to an intermeidate and failed. so we don't include it.
#                 # TODO: Find a better approach?
#                 pass
#             else:
#                 # Rate from i -> j
#                 k_i2j= np.exp(-(barrier-Ei) * beta)
#                 R[rj, ri] += k_i2j
#                 Z[rj, ri] += 1

#                 # # Reverse rate j<-i by detailed balance.
#                 k_j2i = np.exp(-(barrier-Ej) * beta)

#                 R[ri, rj] += k_j2i
#                 Z[ri, rj] += 1

#     for i in range(N):
#         for j in range(N):
#             if Z[i,j] > 0:
#                 R[i,j] *= pref/Z[i,j]

#     for i in range(N):
#         R[i,i] = 0
#         R[i,i] = -np.sum(R[:,i])

#     return R



# if __name__ == "__main__":
#     from scipy.integrate import odeint
#     import argparse
#     parser = argparse.ArgumentParser(
#         description="Summarize reaction products"
#         )

#     parser.add_argument("folder", help="Folder containing the react*** files.")
#     parser.add_argument("--temp", type=float, default=298.15,
#                         help="Temperature for the kinetic rate model.")
#     args = parser.parse_args()
#     folder =args.folder
#     temperature = args.temp

#     kbT = constants.kb * temperature
#     beta = 1/(kbT / constants.hartree_ev)

#     # TODO: make this an argument
#     reactant = xyz2smiles(folder + "/init/opt0000.xyz")
#     pathways, species = read_all_reactions(folder, 0.5 * kbT/constants.hartree_ev)
    
#     print("Setting up a system of rate equations at T=%f K..." % temperature)
#     R = build_rate_matrix(pathways, species, beta)
#     print(np.max(R), np.min(R))
#     print("Chosen reactant : ", reactant)
#     reactant_id = species[reactant]
#     N = len(species)

#     # set initial vector
#     y0 = np.zeros(N)
#     y0[reactant_id] = 1.0

#     # define DE system
#     def f(v,t):
#         return R.dot(v)

#     # Time axis
#     ts = np.array([0, 100, 10000, 1e6, 1e9, 1e12])
#     print(" .... solving rate equations ....")
#     y = odeint(f,y0,ts)
#     spec2 = dict(zip(species.values(), species.keys()))
#     def print_summary(y):
#         # sort the smiles by population
#         pops_smiles = sorted(zip(y,
#             [spec2[i] for i in range(N)]),
#             reverse=True)

#         print("   % total     % products       Molecule")
#         for pop, smi in pops_smiles:
#             if smi == reactant:
#                 ratio = " -------"
#             else:
#                 ratio = "%8.4f"%(100*pop/(1-y[reactant_id]))
#             print("  %8.4f" % (100*pop) +
#                   "       " + ratio+
#                   "       " + smi)
#     # print out summary
#     print("\n\n ======== SUMMARY OF REACTION PRODUCTS =========")
#     print("... at 100 fs ...")
#     print_summary(y[1])

#     print("... at 10 ps ...")
#     print_summary(y[2])

#     print("... at 1 ns ...")
#     print_summary(y[3])

#     print("... at 1 μs ...")
#     print_summary(y[4])

#     print("... at 1 ms ...")
#     print_summary(y[4])



