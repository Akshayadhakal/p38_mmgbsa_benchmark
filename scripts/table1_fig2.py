import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

A = os.path.expanduser('~/benchmark_project/analysis')
os.makedirs(f'{A}/figures', exist_ok=True); os.makedirs(f'{A}/tables', exist_ok=True)
df = pd.read_csv(f'{A}/final_table_d2.csv')
assert len(df) == 29
EXCL = ['lig_p38a_2ee', 'lig_p38a_2l']
rng = np.random.default_rng(2026)

def boot(x, y, n=5000):
    s, p = [], []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.unique(x[i]).size > 1 and np.unique(y[i]).size > 1:
            s.append(spearmanr(x[i], y[i])[0]); p.append(pearsonr(x[i], y[i])[0])
    return np.percentile(s, [2.5, 97.5]), np.percentile(p, [2.5, 97.5])

rows = []
def add(group, label, d, col, sign=1):
    x, y = sign * d[col].values, d.dG_exp_kcal.values
    (sl, sh), (pl, ph) = boot(x, y)
    rows.append(dict(group=group, predictor=label, N=len(d), spearman=spearmanr(x, y)[0], sp_lo=sl, sp_hi=sh,
                     pearson=pearsonr(x, y)[0], pe_lo=pl, pe_hi=ph))

sub = df[~df.ligand.isin(EXCL)]
add('Primary', 'Vina (mean of 3 seeds)', df, 'vina_mean')
add('Primary', 'MM/GBSA, dielectric 1 (mean of 2 replicas)', df, 'mm_avg')
add('Secondary', 'MM/GBSA, dielectric 1, replica 1', df, 'rep1')
add('Secondary', 'MM/GBSA, dielectric 1, replica 2', df, 'rep2')
add('Secondary', 'Vina, without 1 nM ligands', sub, 'vina_mean')
add('Secondary', 'MM/GBSA, dielectric 1, without 1 nM ligands', sub, 'mm_avg')
add('Sensitivity', 'MM/GBSA, dielectric 2 (mean of 2 replicas)', df, 'd2_avg')
add('Sensitivity', 'MM/GBSA, dielectric 2, replica 1', df, 'd2_rep1')
add('Sensitivity', 'MM/GBSA, dielectric 2, replica 2', df, 'd2_rep2')
add('Exploratory', 'Ligand size (minus heavy-atom count)', df, 'heavy_atoms', sign=-1)
t = pd.DataFrame(rows)
t.to_csv(f'{A}/tables/table1.csv', index=False)

out = pd.DataFrame({
    'Group': t.group, 'Score': t.predictor, 'N': t.N,
    'Spearman rho (95% CI)': [f'{r.spearman:.2f} ({r.sp_lo:.2f}, {r.sp_hi:.2f})' for r in t.itertuples()],
    'Pearson r (95% CI)': [f'{r.pearson:.2f} ({r.pe_lo:.2f}, {r.pe_hi:.2f})' for r in t.itertuples()]})
out.to_csv(f'{A}/tables/Table1_formatted.csv', index=False)
print(out.to_string(index=False))

# paired bootstrap, primary comparison
x1, x2, y = df.vina_mean.values, df.mm_avg.values, df.dG_exp_kcal.values
d = []
for _ in range(10000):
    i = rng.integers(0, len(df), len(df))
    if min(np.unique(x1[i]).size, np.unique(x2[i]).size, np.unique(y[i]).size) > 1:
        d.append(spearmanr(x2[i], y[i])[0] - spearmanr(x1[i], y[i])[0])
d = np.array(d)
print('\nMM/GBSA minus Vina (Spearman): %.2f, 95%% CI (%.2f, %.2f); MM/GBSA higher in %.0f%% of resamples'
      % (spearmanr(x2, y)[0] - spearmanr(x1, y)[0], *np.percentile(d, [2.5, 97.5]), 100 * (d > 0).mean()))
print('Vina score range: %.2f to %.2f kcal/mol' % (df.vina_mean.min(), df.vina_mean.max()))

# Figure 2: forest plot
col = {'Primary': '#1f4e79', 'Secondary': '#5b8fb9', 'Sensitivity': '#c0504d', 'Exploratory': '#7f7f7f'}
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ypos = np.arange(len(t))[::-1]
for yy, r in zip(ypos, t.itertuples()):
    ax.errorbar(r.spearman, yy, xerr=[[r.spearman - r.sp_lo], [r.sp_hi - r.spearman]], fmt='o', ms=4, color=col[r.group], capsize=3)
ax.axvline(0, color='k', lw=0.8)
ax.set_yticks(ypos); ax.set_yticklabels([f'{r.predictor} (N={r.N})' for r in t.itertuples()], fontsize=7.5)
ax.set_xlabel('Spearman correlation with experimental \u0394G (95% bootstrap CI)')
for g, c in col.items(): ax.plot([], [], 'o', color=c, label=g)
ax.legend(fontsize=7, loc='lower right', frameon=False)
plt.tight_layout()
plt.savefig(f'{A}/figures/Figure2_forest.png', dpi=300); plt.savefig(f'{A}/figures/Figure2_forest.pdf')
print('\nsaved Table1_formatted.csv, Figure2_forest.png/.pdf')
