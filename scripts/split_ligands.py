import os
from rdkit import Chem
home = os.path.expanduser("~")
sdf = f"{home}/benchmark_project/protein-ligand-benchmark/data/p38/02_ligands/ligands.sdf"
out = f"{home}/benchmark_project/docking/ligands"
os.makedirs(out, exist_ok=True)
for m in Chem.SDMolSupplier(sdf, removeHs=False):
    name = m.GetProp("_Name").strip()
    Chem.MolToMolFile(m, f"{out}/{name}.mol")
    print(name, m.GetNumAtoms())
