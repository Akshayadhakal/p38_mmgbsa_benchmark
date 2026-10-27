import os, numpy as np
from rdkit import Chem
home = os.path.expanduser("~")
sdf = f"{home}/benchmark_project/protein-ligand-benchmark/data/p38/02_ligands/ligands.sdf"
allc, ext = [], []
for m in Chem.SDMolSupplier(sdf, removeHs=False):
    if m is None: continue
    xyz = m.GetConformer().GetPositions()
    allc.append(xyz)
    ext.append(xyz.max(0) - xyz.min(0))
pts = np.vstack(allc)
center = (pts.max(0) + pts.min(0)) / 2
size = pts.max(0) - pts.min(0) + 10  # 10 A margin
print("ligands read:", len(allc))
print("center  x y z: %.2f %.2f %.2f" % tuple(center))
print("size    x y z: %.1f %.1f %.1f" % tuple(size))
print("largest single-ligand extent: %.1f %.1f %.1f" % tuple(np.max(ext, axis=0)))
