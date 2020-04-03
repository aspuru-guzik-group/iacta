# Various IO routines
import numpy as np
import re
import os

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

def comment_line_energy(comment_line):
    m = re.search('-?[0-9]*\.[0-9]*', comment_line)
    if m:
        E = float(m.group())
    else:
        E = np.nan
    return E

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
            E = comment_line_energy(comment_line)

            for i in range(natoms):
                this_mol += f.readline()

            if index is None:
                structures += [this_mol]
                energies += [E]
            else:
                if k == index:
                    return this_mol
                
            k+=1
    return structures,energies

def xyz2numpy(input, natoms):
    """Read an xyz file."""
    if type(input) == str:
        lines = input.split("\n")
    elif type(input) == list:
        lines = input

    atoms = []
    positions = np.zeros((natoms, 3))
    iatom = 0
    for line in input:
        positions[iatom,:] = np.fromstring(line[2:],
                                       count=3, sep=" ")
        atoms += [line[0:2].rstrip().lstrip()]
        iatom+=1
        if iatom == natoms:
            break
    return atoms, positions

def read_xyz(filepath, index=0):
    """Read an xyz file and convert to numpy arrays."""
    atoms = []
    xyzs = []
    comments = []
    with open(filepath, 'r') as f:
        curr = 0
        while True:
            first_line = f.readline()
            # EOF -> blank line
            if not first_line:
                break
            natoms = int(first_line.rstrip())
            comment_line = f.readline()

            if index=="all":
                at, r = xyz2numpy(f,natoms)
                atoms += [at]
                xyzs += [r]
                comments += [comment_line]

            elif curr == index:
                atoms, pos = xyz2numpy(f, natoms)
                return atoms, pos, comment_line
            else:
                for i in range(natoms):
                    f.readline()
                curr += 1
    return atoms,xyzs,comments

def dump_succ_opt(output_folder, structures, energies, opt_indices,
                  concat=False,
                  extra=False):
    os.makedirs(output_folder, exist_ok=True)
    

    if concat:
        # Dump the optimized structures in one file                
        with open(output_folder + "/opt.xyz", "w") as f:
            for oi in opt_indices:
                f.write(structures[oi])
    else:
        # Dump the optimized structures in many files            
        for stepi, oi in enumerate(opt_indices):
            with open(output_folder + "/opt%4.4i.xyz" % stepi, "w") as f:
                f.write(structures[oi])

    # Dump indices of optimized structures, energies of optimized structures,
    # all energies and all structures
    np.savetxt(output_folder + "/Eopt", np.array(energies)[opt_indices], fmt="%15.8f")

    if extra:
        np.savetxt(output_folder + "/indices", opt_indices, fmt="%i")
        np.savetxt(output_folder + "/E", energies, fmt="%15.8f")
        with open(output_folder + "/log.xyz", "w") as f:
            for s in structures:
                f.write(s)





