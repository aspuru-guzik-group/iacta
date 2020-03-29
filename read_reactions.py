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

def xyz2smiles(xyzfile):
    output = []
    for m in pybel.readfile("xyz",xyzfile):
        output+= [m.write(format="smi", opt={"c":1,"n":1}).rstrip()]
    if len(output) == 1:
        return output[0]
    else:
        return output

def analyze_reaction(react_folder,
                     smallest_barrier=1e-3):
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
    pathway = []
    minima = []
    barriers = []
    last = 0
    last_smiles = smiles[0]

    # loop through smiles and detect changes
    for si,s in enumerate(smiles):
        if s == last_smiles:
            pass
        else:
            # we have a state change at index=si
            # minimum energy point for the current state
            minE_point = np.argmin(E[last:si]) + last
            minE = E[minE_point]
            barrier_height = np.max(E[minE_point:si])
            # We only record the change if the barrier exists.
            if barrier_height - minE > smallest_barrier:
                minima += [E[minE_point]]
                # Energy barrier height
                # maximum of energy between the minimum and the state change
                barriers += [np.max(E[minE_point:si])]
                pathway += [smiles[minE_point]]

            last_smiles = s
            last = si

    # If the last point was not itself a reaction, then we add the last
    # molecule to the list.
    if last == si:
        pass
    else:
        minE_point = np.argmin(E[last:si]) + last
        minima += [E[minE_point]]
        pathway += [smiles[minE_point]]

    if len(pathway) == 0:
        print("warning: problem with trajectory in ", react_folder)
        return None
        
    # If the same molecule appears twice in a row, that means we have a
    # conformation change with a barrier, but still no reaction
    npathway = [pathway[0]]
    nminima = [minima[0]]

    if len(barriers)>0:
        # we have some sort of reaction or conformer change with barriers
        nbarriers = [barriers[0]]

        for i in range(1,len(pathway)):
            if pathway[i] == pathway[i-1]:
                # No transformation, just another conformer

                # update barriers with the current max barrier and minima with the
                # current min minima.
                nminima[-1] = min(nminima[-1], minima[i])
                if i< len(pathway)-1:                
                    nbarriers[-1] = max(nbarriers[-1], barriers[i])
            else:
                npathway += [pathway[i]]
                nminima += [minima[i]]

                if i< len(pathway)-1:
                    nbarriers += [barriers[i]]
    else:
        nbarriers = []



    return ReactionPathway(npathway, nminima, nbarriers, react_folder)


def read_all_reactions(output_folder, verbose=True):
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
        rpath = analyze_reaction(f)
        if rpath:
            if len(rpath.species)>1:
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
        print("  with %i unique chemical species" % len(species))
    return pathways,species



pathways, species = read_all_reactions("output")

T = 300.0
print("Setting up a system of rate equations at T=%f..." % T)

kbT = constants.kb * T # (eV)
pref = kbT/constants.hbar * 1e6        # ns^-1

N = len(species)
R = np.zeros((N,N))
Z = np.zeros((N,N))             # number of pathways contributing to i->j
for p in pathways:
    for i in range(len(p.species)-1):
        ri = species[p.species[i]]
        Ei = p.minima[i]
        rj = species[p.species[i+1]]
        Ej = p.minima[i+1]
        barrier = p.barriers[i]

        if ri == rj:
            print("wat")

        # Rate from i -> j
        k_ij= np.exp(-(barrier-Ei) * constants.hartree_ev/kbT)
        R[rj, ri] += k_ij
        Z[rj, ri] += 1
        
        # Reverse rate j<-i by detailed balance.
        k_ji= k_ij * np.exp(-(Ei-Ej) * constants.hartree_ev/kbT)
        R[ri, rj] += k_ji
        Z[ri, rj] += 1

for i in range(N):
    for j in range(N):
        if Z[i,j] > 0:
            R[i,j] *= pref/Z[i,j]
            
for i in range(N):
    R[i,i] = 0
    R[i,i] = -np.sum(R[:,i])

reactant = xyz2smiles("./output/init/opt0000.xyz")
print("Chosen reactant : ", reactant)
reactant_id = species[reactant]

y0 =np.zeros(N)
y0[reactant_id] = 1.0

def f(v,t):
    return R.dot(v)

from scipy.integrate import odeint
ts = np.logspace(-1,6, 1000)

print(" .... solving rate equations ....")
y = odeint(f,y0,ts)

# short time
y_short = y[1]
y_long = y[-1]


spec2 = dict(zip(species.values(), species.keys()))
# silly...
pops_smiles_short = sorted(zip(y_short,
                               [spec2[i] for i in range(N)]),
                           reverse=True)

pops_smiles_long = sorted(zip(y_long,
                               [spec2[i] for i in range(N)]),
                           reverse=True)


print("\n\n ======== SUMMARY OF REACTION PRODUCTS =========")
print("... t → 0 products ...")
print("   % total     % products        Molecule")
for pop, smi in pops_smiles_short:
    if smi == reactant:
        ratio = "  ---  "
    else:
        ratio = "%7.4f"%(100*pop/(1-y_short[reactant_id]))
    print("   %7.4f" % (100*pop) +
          "        " + ratio+
          "        " + smi)

print("\n")
print("... → ∞ products ...")
print("   % total     % products        Molecule")
for pop, smi in pops_smiles_long:
    if smi == reactant:
        ratio = "       "
    else:
        ratio = "%7.4f"%(100*pop/(1-y_long[reactant_id]))
    print("   %7.4f" % (100*pop) +
          "        " + ratio+
          "        " + smi)
