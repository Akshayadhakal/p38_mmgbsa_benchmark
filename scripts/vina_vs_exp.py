import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
home = os.path.expanduser("~")
exp = pd.read_csv(f"{home}/benchmark_project/analysis/p38_ligands.csv")[["name","dG_exp_kcal"]]
v = pd.read_csv(f"{home}/benchmark_project/docking/vina_scores.csv")
v = v[v.name != "DONE"].copy()
v["vina_top"] = pd.to_numeric(v.vina_top)
m = v.groupby("name").vina_top.mean().rename("vina_mean").reset_index().merge(exp, on="name")
rng = np.random.default_rng(0)
def stats(d):
    return spearmanr(d.vina_mean, d.dG_exp_kcal)[0], pearsonr(d.vina_mean, d.dG_exp_kcal)[0]
def boot(d, n=5000):
    s = []
    for _ in range(n):
        b = d.sample(len(d), replace=True, random_state=int(rng.integers(1e9)))
        if b.vina_mean.nunique() > 1 and b.dG_exp_kcal.nunique() > 1:
            s.append(stats(b))
    s = np.array(s)
    return np.percentile(s, [2.5, 97.5], axis=0)
for label, d in [("all 29", m), ("without 2ee, 2l (1 nM)", m[~m.name.isin(["lig_p38a_2ee","lig_p38a_2l"])])]:
    sp, pe = stats(d); ci = boot(d)
    print(f"{label}: Spearman {sp:.2f} [{ci[0,0]:.2f},{ci[1,0]:.2f}]  Pearson {pe:.2f} [{ci[0,1]:.2f},{ci[1,1]:.2f}]  N={len(d)}")
m.to_csv(f"{home}/benchmark_project/analysis/vina_vs_exp.csv", index=False)
