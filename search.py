import numpy as np
from rdkit_utils import *
from xtb_utils import *


    
XTB_PATH = "/h/182/clavigne/xtb/bin/xtb"
def xtb_optimize(mol):
    with tempfile.NamedTemporaryFile(suffix=".mol") as t:
        Chem.MolToMolFile(mol,t.name)
        opt = xtb_run(XTB_PATH, t.name, "--opt", cwd=".")
        opt.proc.wait()
        out = Chem.MolFromMolFile(opt.dir + "/xtbopt.mol")
        opt.close()
        return out

xtb = xtb_driver(XTB_PATH)

# Molecules
smiles_react1="CNNC(C)(C)"
smiles_react2="C1CCCC1([Br])"


# Build conformers
react1 = rdkit_generate_conformer(smiles_react1)
react2 = rdkit_generate_conformer(smiles_react2)
combined = mix2(react1,react2, off=5.0)

# xtb parameters
combined = xtb_optimize(combined)
Chem.MolToMolFile(combined, "out.mol")



              

