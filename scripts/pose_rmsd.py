import os, glob, subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

home = os.path.expanduser("~")
dock = f"{home}/benchmark_project/docking"
A = f"{home}/benchmark_project/analysis"
rows = []
for molf in sorted(glob.glob(f"{dock}/ligands/*.mol")):
    name = os.path.basename(molf)[:-4]
    ref = Chem.RemoveHs(Chem.MolFromMolFile(molf, removeHs=False))
    for seed in (1, 2, 3):
        pq = f"{dock}/out/{name}_s{seed}.pdbqt"
        pdb = f"/tmp/{name}_s{seed}_top.pdb"
        rms, err = float("nan"), ""
        try:
            subprocess.run(["obabel", pq, "-f", "1", "-l", "1", "-O", pdb], check=True, stderr=subprocess.DEVNULL)
            pose = Chem.MolFromPDBFile(pdb, removeHs=True, proximityBonding=True)
            pose = AllChem.AssignBondOrdersFromTemplate(ref, pose)
            rms = rdMolAlign.CalcRMS(pose, ref)
        except Exception as e:
            err = type(e).__name__
        rows.append(dict(ligand=name, seed=seed, rmsd=rms, error=err))
df = pd.DataFrame(rows)
df.to_csv(f"{A}/vina_pose_rmsd.csv", index=False)
ok = df.dropna(subset=["rmsd"])
print("poses scored: %d of %d" % (len(ok), len(df)))
for e, g in df[df.error != ""].groupby("error"):
    print("  failed (%s): %s" % (e, ", ".join(sorted(set(g.ligand)))))
for s, g in ok.groupby("seed"):
    print("seed %d: %d of %d within 2 A | median %.2f A | max %.2f A" % (s, (g.rmsd < 2).sum(), len(g), g.rmsd.median(), g.rmsd.max()))
w = ok.pivot(index="ligand", columns="seed", values="rmsd")
w["mean"] = w[[1, 2, 3]].mean(axis=1)
w["spread"] = w[[1, 2, 3]].max(axis=1) - w[[1, 2, 3]].min(axis=1)
print("ligands with all three seeds within 2 A: %d of %d" % ((w[[1, 2, 3]] < 2).all(axis=1).sum(), len(w)))
print("\nworst 8 ligands (RMSD in A):")
print(w.sort_values("mean", ascending=False).head(8).round(2).to_string())
if "lig_p38a_3fly" in w.index:
    print("\n3FLY crystal ligand, seed 1: %.2f A (the earlier redocking check gave 0.81)" % w.loc["lig_p38a_3fly", 1])
