# Various IO routines
import numpy as np
import re
import os
import openbabel.pybel as pybel
import subprocess
import re

def metadata():
    # Return a dictionary with some metadata to improve reproducibility
    xtbv = subprocess.check_output(["xtb", "--version"],
                                   stderr=subprocess.DEVNULL).strip().decode()
    for line in xtbv.split("\n"):
        if re.search("version", line):
            xtbvl = line.strip()
    return {
        "hostname":subprocess.check_output(
            ["hostname"]).strip().decode(),
        "xtb":xtbvl,
        "commit":subprocess.check_output(
            ["git", "describe", "--always"],
            cwd=os.path.dirname(__file__)).strip().decode()}

# =================== xTB output  reading/writing routines =====================
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

# =================== xyz trajectory files reading/writing routines =============================
def traj2str(filepath, index=None, as_list=False):
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
            E = comment_line_energy(comment_line)

            for i in range(natoms):
                this_mol += f.readline()

            if index is None:
                structures += [this_mol]
                energies += [E]
            else:
                if k == index:
                    if as_list:
                        return [this_mol], [E]
                    else:
                        return this_mol, E
                
            k+=1
    return structures,energies

def traj2smiles(filepath, index=None, chiral=False):
    """Read an xyz file and convert to a list of SMILES ."""
    # Read the trajectory
    strs, E = traj2str(filepath, index=index, as_list=True)
    output = []
    
    if chiral:
        flags = {"c":1,"n":1}
    else:
        flags = {"c":1,"n":1, "i":1}
    for s in strs:
        output+= [pybel.readstring("xyz", s).write(format="smi", opt=flags).rstrip()]

    if index is None:
        return output, E
    else:
        return output[0], E[0]

def traj2npy(filepath, index=None):
    """Read an xyz file and convert to numpy arrays."""
    # Read the trajectory
    strs, E = traj2str(filepath, index=index, as_list=True)
    
    atoms = []
    positions = []
    for s in strs:
        at, r = xyz2numpy(s)
        atoms += [at]
        positions += [r]
        
    if index is None:
        return atoms, positions, E
    else:
        return atoms[0], positions[0], E[0]

def comment_line_energy(comment_line):
    m = re.search('-?[0-9]*\.[0-9]*', comment_line)
    if m:
        E = float(m.group())
    else:
        E = np.nan
    return E

def xyz2numpy(string):
    """Convert xyz file as a string to a numpy array."""
    # remove the two lines of the header and any empty lines at the end
    lines = [l for l in string.split("\n")[2:] if l]
    atoms = []
    positions = np.zeros((len(lines), 3))
    for iatom, line in enumerate(lines):
        positions[iatom,:] = np.fromstring(line[2:], count=3, sep=" ")
        atoms += [line[0:2].rstrip().lstrip()]
    return atoms, positions
