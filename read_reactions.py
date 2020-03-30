import pybel
import numpy as np
import glob
import constants

class ReactionPathway:
    def __init__(self, species, minima, barriers, id):
        self.species = tuple(species)
        self.minima = minima
        self.barriers = barriers
        self.id = id
        self.nreactions = len(set(self.species))

def xyz2smiles(xyzfile, chiral=False):
    output = []
    for m in pybel.readfile("xyz",xyzfile):
        output+= [m.write(format="smi", opt={"c":1,"n":1,
                                             "i":chiral*1}).rstrip()]
    if len(output) == 1:
        return output[0]
    else:
        return output

def read_reaction(react_folder, barrier_tol):
    """Extract chemical quantities from a reaction.

    Parameters:
    -----------
    react_folder (str) : folder storing the reaction data, obtained from the
    reaction_job() routine in react.py.

    Optional parameters:
    --------------------
    smallest_barrier (float) : smallest barrier height to consider a molecule
    a true metastable intermediate.

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
    mols = [smiles[0]]
    regions = []
    rstart = 0

    # loop through smiles and detect changes
    for si,s in enumerate(smiles):
        if s == mols[-1]:
            pass
        else:
            mols += [s]
            regions += [(rstart,si)]
            rstart = si

    regions += [(rstart, len(smiles))]

    # Refine regions by merging those that have barriers less than
    # barrier_tol. This is nasty but it works.
    refined = []
    while True:
        if len(regions) == 0:
            break
           
        if len(regions) == 1:
            refined += [regions.pop()]
            break

        current = regions.pop(0)
        current_imin = np.argmin(E[current[0]:current[1]]) + current[0]
        current_Emin = E[current_imin]
        # iterate over remaining regions
        while regions:
            future = regions[0]
            future_imin = np.argmin(E[future[0]:future[1]]) + future[0]
            future_Emin = E[future_imin]

            # barrier is the max between the two minima (inclusive
            barrier = np.max(E[current_imin:future_imin+1])

            # if the barrier is to small, we merge current and future regions.
            if ((barrier < (current_Emin + barrier_tol))
                or (barrier < (future_Emin + barrier_tol))):
                current = (current[0], future[1])

                # current_imin = np.argmin(E[current[0]:current[1]]) + current[0]
                # current_Emin = E[current_imin]
                
                # if the minima of future is below the minima of current, we
                # switch the minima of current to the minima of future.
                if future_Emin < current_Emin:
                    current_Emin = future_Emin
                    current_imin = future_imin
                    
                # pop the one we just melded in
                regions.pop(0)


            # Otherwise, we break out of the inner loop and save current
            else:
                break
        refined +=[current[:]]

    # Save all those quantities
    Emin = []
    imins = []
    mols = []
    for r in refined:
        imin = np.argmin(E[r[0]:r[1]]) + r[0]
        imins += [imin]
        Emin += [E[imin]]
        mols += [smiles[imin]]
    barrier = []
    for k in range(len(refined)-1):
        barrier += [np.max(E[imins[k]:imins[k+1]+1])]
            
    return ReactionPathway(mols, Emin, barrier, react_folder)

def read_all_reactions(output_folder, Etol, verbose=True):
    """Read and parse all reactions in a given folder."""
    folders = glob.glob(output_folder + "/react[0-9]*")
    if verbose:
        print("Parsing folder <%s>, with" % output_folder)
        print("   %6i trajectories..." % len(folders))

    failed = []
    noreact = []
    pathways = []
    species = {}
    index = 0
    for f in folders:
        rpath = read_reaction(f, Etol)
        if rpath:
            if rpath.nreactions>1:
                for spec in rpath.species:
                    if spec in species:
                        pass
                    else:
                        species[spec] = index
                        index+=1
                pathways += [rpath]
            else:
                # No reaction in this pathway
                noreact += [f]
        else:
            # Convergence failed
            failed += [f]
            
    if verbose:
        print(" - %6i that did not converge" % len(failed))
        print(" - %6i with no reactions" % len(noreact))
        print("--------------")
        print(" = %6i reactions found" % len(pathways))
        print("    with %6i distinct stable species" % len(species))
    return pathways,species

# kbT / hbar = 1/beta hbar = 1/beta (Eh^-1) * 1/hbar(eV fs) * 27 eV / hartree
EYRING_PREF = 1/constants.hbar * constants.hartree_ev # fs^-1
def build_rate_matrix(pathways, species, beta):
    """Build the matrix of rate from pathways."""
    N = len(species)
    R = np.zeros((N,N))
    Z = np.zeros((N,N))             # number of pathways contributing to i->j
    pref = EYRING_PREF / beta
    for p in pathways:
        for i in range(len(p.species)-1):
            # species index into R
            ri = species[p.species[i]]
            rj = species[p.species[i+1]]
            
            # minima
            Ei = p.minima[i]
            Ej = p.minima[i+1]
            
            # barrier
            barrier = p.barriers[i]

            if ri == rj:
                # This is a pathway between two conformers or a reaction that
                # went to an intermeidate and failed. so we don't include it.
                # TODO: Find a better approach?
                pass
            else:
                # Rate from i -> j
                k_i2j= np.exp(-(barrier-Ei) * beta)
                R[rj, ri] += k_i2j
                Z[rj, ri] += 1

                # # Reverse rate j<-i by detailed balance.
                k_j2i = np.exp(-(barrier-Ej) * beta)

                R[ri, rj] += k_j2i
                Z[ri, rj] += 1

    for i in range(N):
        for j in range(N):
            if Z[i,j] > 0:
                R[i,j] *= pref/Z[i,j]

    for i in range(N):
        R[i,i] = 0
        R[i,i] = -np.sum(R[:,i])

    return R



if __name__ == "__main__":
    from scipy.integrate import odeint
    import argparse
    parser = argparse.ArgumentParser(
        description="Summarize reaction products"
        )

    parser.add_argument("folder", help="Folder containing the react*** files.")
    parser.add_argument("--temp", type=float, default=298.15,
                        help="Temperature for the kinetic rate model.")
    args = parser.parse_args()
    folder =args.folder
    temperature = args.temp

    kbT = constants.kb * temperature
    beta = 1/(kbT / constants.hartree_ev)

    # TODO: make this an argument
    reactant = xyz2smiles(folder + "/init/opt0000.xyz")
    pathways, species = read_all_reactions(folder, 0.5 * kbT/constants.hartree_ev)
    
    print("Setting up a system of rate equations at T=%f K..." % temperature)
    R = build_rate_matrix(pathways, species, beta)
    print(np.max(R), np.min(R))
    print("Chosen reactant : ", reactant)
    reactant_id = species[reactant]
    N = len(species)

    # set initial vector
    y0 = np.zeros(N)
    y0[reactant_id] = 1.0

    # define DE system
    def f(v,t):
        return R.dot(v)

    # Time axis
    ts = np.array([0, 100, 10000, 1e6, 1e9, 1e12])
    print(" .... solving rate equations ....")
    y = odeint(f,y0,ts)
    spec2 = dict(zip(species.values(), species.keys()))
    def print_summary(y):
        # sort the smiles by population
        pops_smiles = sorted(zip(y,
            [spec2[i] for i in range(N)]),
            reverse=True)

        print("   % total     % products       Molecule")
        for pop, smi in pops_smiles:
            if smi == reactant:
                ratio = " -------"
            else:
                ratio = "%8.4f"%(100*pop/(1-y[reactant_id]))
            print("  %8.4f" % (100*pop) +
                  "       " + ratio+
                  "       " + smi)
    # print out summary
    print("\n\n ======== SUMMARY OF REACTION PRODUCTS =========")
    print("... at 100 fs ...")
    print_summary(y[1])

    print("... at 10 ps ...")
    print_summary(y[2])

    print("... at 1 ns ...")
    print_summary(y[3])

    print("... at 1 μs ...")
    print_summary(y[4])

    print("... at 1 ms ...")
    print_summary(y[4])



