import pybel
import openbabel as ob
import argparse

parser = argparse.ArgumentParser(
    description="Magically make some reactants from a SMILES string and a SMARTS match.")

parser.add_argument("SMILES",
                    help="SMILES string for the reactant, converted to an xyz"
                    +" file using openbabel.", type=str)
parser.add_argument("SMARTS",
                    help="SMARTS string to identify the atom numbers in"
                    " the resulting file.", type=str)
parser.add_argument("-o",
                    help="Output xyz file. Defaults reactants.xyz.", type=str,
                    default="reactants.xyz")
parser.add_argument("--ff",
                    help="Forcefield to use. Defaults to UFF. If \"no\", no optimization is performed at all.", type=str,
                    default="uff")
parser.add_argument("--steps",
        help="Number of steps of optimization to perform. defaults to 50.",
        type=int, default=50)

args = parser.parse_args()
reactants = args.SMILES
reaction = pybel.Smarts(args.SMARTS)

# Make a conformer
reactants_mol = pybel.readstring("smi",reactants)
if args.ff == "no":
    reactants_mol.make3D(steps=0)
else:
    reactants_mol.make3D(forcefield=args.ff, steps=args.steps)

# write out to stdout
print("--------------------- XYZ FILE -------------------")
xyz = reactants_mol.write("xyz", filename=args.o, overwrite=True)
print("         writing to <%s>..."%args.o)
print("")
print("SMARTS string: %s"%args.SMARTS)
match = reaction.findall(reactants_mol)

if match:
    print("The following atoms were found to match:")

    for atoms in match:
        str =   "Atoms:    "
        index = "Indices:  "
        for i in atoms:
            new = reactants_mol.atoms[i-1].type + "    "
            str += new
            new2 = "%i" % i
            index += new2 + " " * (len(new) - len(new2))
        print(str)
        print(index)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~")



else:
    print("No match!")



