import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A = os.path.expanduser('~/benchmark_project/analysis')
exp = pd.read_csv(f'{A}/p38_ligands.csv')[['name', 'dG_exp_kcal', 'heavy_atoms']].rename(columns={'name': 'ligand'})
vina = pd.read_csv(f'{A}/vina_vs_exp.csv')[['name', 'vina_mean']].rename(columns={'name': 'ligand'})
mm = pd.read_csv(f'{A}/mmpbsa_table.csv')
w = mm.pivot(index='ligand', columns='rep', values='dTOTAL')
w.columns = ['rep1', 'rep2']
w['mm_avg'] = w[['rep1', 'rep2']].mean(axis=1)
df = exp.merge(vina, on='ligand').merge(w.reset_index(), on='ligand')
assert len(df) == 29, len(df)
rng = np.random.default_rng(12345)
def sp(x, y): return spearmanr(x, y)[0]
def boot_ci(x, y, n=5000):
    out, N = [], len(x)
    for _ in range(n):
        i = rng.integers(0, N, N)
        if np.unique(x[i]).size > 1 and np.unique(y[i]).size > 1: out.append(sp(x[i], y[i]))
    return np.percentile(out, [2.5, 97.5])
def report(label, d, col):
    x, yy = d[col].values, d.dG_exp_kcal.values
    r = sp(x, yy); lo, hi = boot_ci(x, yy)
    print(f'{label:<38} N={len(d):2d}  Spearman {r:5.2f} [{lo:5.2f}, {hi:5.2f}]  Pearson {pearsonr(x, yy)[0]:5.2f}')
    return r
y = df.dG_exp_kcal.values
print('PRIMARY')
rv = report('Vina (mean of 3 seeds)', df, 'vina_mean')
rm = report('MM/GBSA (mean of 2 replicas)', df, 'mm_avg')
x1, x2, N, diffs = df.vina_mean.values, df.mm_avg.values, len(df), []
for _ in range(10000):
    i = rng.integers(0, N, N)
    if min(np.unique(x1[i]).size, np.unique(x2[i]).size, np.unique(y[i]).size) > 1:
        diffs.append(sp(x2[i], y[i]) - sp(x1[i], y[i]))
diffs = np.array(diffs)
lo, hi = np.percentile(diffs, [2.5, 97.5])
print('difference (MM/GBSA - Vina): %.2f, 95%% CI [%.2f, %.2f]; resamples with MM/GBSA better: %.0f%%' % (rm - rv, lo, hi, 100 * (diffs > 0).mean()))
print('\nSECONDARY')
sub = df[~df.ligand.isin(['lig_p38a_2ee', 'lig_p38a_2l'])]
report('Vina, without 2ee and 2l', sub, 'vina_mean')
report('MM/GBSA average, without 2ee and 2l', sub, 'mm_avg')
report('MM/GBSA replica 1 only', df, 'rep1')
report('MM/GBSA replica 2 only', df, 'rep2')
print('\nEXPLORATORY: Spearman with ligand heavy-atom count')
for lab, col in (('experimental dG', 'dG_exp_kcal'), ('Vina', 'vina_mean'), ('MM/GBSA average', 'mm_avg')):
    print(f'  {lab:<16} {sp(df[col].values, df.heavy_atoms.values):5.2f}')
fig, ax = plt.subplots(1, 2, figsize=(9, 4))
ax[0].scatter(df.vina_mean, df.dG_exp_kcal); ax[0].set_xlabel('Vina score (kcal/mol)'); ax[0].set_title('Spearman %.2f' % rv)
ax[1].errorbar(df.mm_avg, df.dG_exp_kcal, xerr=(df.rep1 - df.rep2).abs() / 2, fmt='o'); ax[1].set_xlabel('MM/GBSA dTOTAL, mean of replicas (kcal/mol)'); ax[1].set_title('Spearman %.2f' % rm)
for a in ax: a.set_ylabel('Experimental dG from IC50 (kcal/mol)')
plt.tight_layout(); plt.savefig(f'{A}/scatter_primary.png', dpi=200)
df.to_csv(f'{A}/final_table.csv', index=False)
print('\nsaved final_table.csv and scatter_primary.png in', A)
