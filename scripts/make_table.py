import os, sys, yaml, math
import pandas as pd
from rdkit import Chem

target = sys.argv[1] if len(sys.argv) > 1 else "p38"
home = os.path.expanduser("~")
src = f"{home}/benchmark_project/protein-ligand-benchmark/data/{target}/00_data/ligands.yml"
data = yaml.safe_load(open(src))

to_molar = {"nM": 1e-9, "uM": 1e-6, "µM": 1e-6, "pM": 1e-12, "mM": 1e-3, "M": 1.0}
RT = 0.0019872 * 298.15  # kcal/mol at 298.15 K

rows = []
for name, d in data.items():
    m = d["measurement"]
    conc = m["value"] * to_molar[m["unit"]]
    mol = Chem.MolFromSmiles(d["smiles"])
    rows.append({
        "name": name,
        "type": m["type"],
        "value": m["value"],
        "unit": m["unit"],
        "error": m.get("error"),
        "dG_exp_kcal": RT * math.log(conc),
        "formal_charge": Chem.GetFormalCharge(mol) if mol else None,
        "heavy_atoms": mol.GetNumHeavyAtoms() if mol else None,
        "smiles": d["smiles"],
    })

df = pd.DataFrame(rows).sort_values("dG_exp_kcal")
out = f"{home}/benchmark_project/analysis/{target}_ligands.csv"
df.to_csv(out, index=False)
print(df.drop(columns="smiles").to_string(index=False))
print("\nN ligands:", len(df))
print("dG range (kcal/mol): %.2f to %.2f" % (df.dG_exp_kcal.min(), df.dG_exp_kcal.max()))
print("Formal charges:", df.formal_charge.value_counts().to_dict())
