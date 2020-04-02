# Various IO routines
import numpy as np
import re

# File reading/writing routines
def read_wbo(filepath):
    bonds = []
    with open(filepath, 'r') as f:
        for line in f:
            i,j, wbo = line.split()
            # These are 1-indexed
            bonds += [(int(i)-1, int(j)-1, float(wbo))]
    return bonds

def read_charges(filepath):
    charges = []
    with open(filepath, 'r') as f:
        for line in f:
            charges += [float(line.rstrip())]
    return charges

def read_xtb_output(xyzfile):
    dir = os.path.dirname(xyzfile)
    if not dir:
        dir = "."
    # Read in the xyz file
    atoms, positions = read_xyz(xyzfile)
    # Read in wbo file
    wbo = read_wbo(dir  + "/wbo")
    # Read in partial charges
    charges = read_charges(dir + "/charges")
    return atoms, charges, positions, wbo

def read_xtb_gradient(filen):
    output = []
    with open(filen, "r") as f:
        for line in f:
            if "$gradient" in line:
                # start of the gradient group
                break
            
        gradient = ""
        n = 0
        for line in f:
            if "$end" in line:
                break
            
            # The non-gradient line have atoms at the end
            char = re.search(" [A-z]", line)
            if char:
                pass
            else:
                gradient += line
                n+=1
                
    return np.fromstring(gradient, sep=" ").reshape((n,3))

def read_trajectory(filepath, index=None):
    """Read an xyz file containing a trajectory."""
    structures = []
    energies = []
    k = 0
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
            # first number on comment_line
            m = re.search('-?[0-9]*\.[0-9]*', comment_line)       
            energies += [float(m.group())]
        
            for i in range(natoms):
                this_mol += f.readline()

            if index is None:
                structures += [this_mol]
            else:
                if k == index:
                    return this_mol
                
            k+=1
    return structures,energies

def read_xyz(filepath, index=0):
    """Read an xyz file."""
    with open(filepath, 'r') as f:
        curr = 0
        while True:
            first_line = f.readline()
            # EOF -> blank line
            if not first_line:
                break
            natoms = int(first_line.rstrip())
            comment_line = f.readline()

            if curr == index:
                atoms = []
                positions = np.zeros((natoms, 3))
                for i in range(natoms):
                    line = f.readline()
                    positions[i,:] = np.fromstring(line[2:],
                                                        count=3, sep=" ")
                    atoms += [line[0:2].rstrip().lstrip()]
                return atoms, positions
            else:
                for i in range(natoms):
                    f.readline()
                curr += 1

