import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

A = os.path.expanduser('~/benchmark_project/analysis'); M = os.path.expanduser('~/benchmark_project/md')
F = f'{A}/figures'
df = pd.read_csv(f'{A}/final_table_d2.csv'); t1 = pd.read_csv(f'{A}/tables/table1.csv')
short = lambda s: s.replace('lig_p38a_', '')
BLUE, RED = '#1f4e79', '#c0504d'
def save(name): plt.tight_layout(); plt.savefig(f'{F}/{name}.png', dpi=300); plt.savefig(f'{F}/{name}.pdf'); plt.close()
lab = lambda r: '\u03c1 = %.2f (%.2f, %.2f)' % (r.spearman, r.sp_lo, r.sp_hi)

# Figure 3: scores versus experiment
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4), sharex=True)
ax[0].scatter(df.dG_exp_kcal, df.vina_mean, s=22, color=BLUE)
ax[0].set_ylabel('Vina score (kcal/mol)'); ax[0].text(0.04, 0.95, lab(t1.iloc[0]), transform=ax[0].transAxes, va='top', fontsize=8)
ax[1].errorbar(df.dG_exp_kcal, df.mm_avg, yerr=(df.rep1 - df.rep2).abs() / 2, fmt='o', ms=4, color=BLUE, ecolor='#9db7d1', capsize=2)
ax[1].set_ylabel('MM/GBSA \u0394TOTAL (kcal/mol)'); ax[1].text(0.04, 0.95, lab(t1.iloc[1]), transform=ax[1].transAxes, va='top', fontsize=8)
for a, l in zip(ax, 'ab'):
    a.set_xlabel('Experimental \u0394G from IC$_{50}$ (kcal/mol)'); a.set_title(l, loc='left', fontweight='bold')
save('Figure3_scores_vs_experiment')

# Figure 4: replica agreement
diff = df.rep1 - df.rep2
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4))
lo, hi = min(df.rep1.min(), df.rep2.min()) - 1, max(df.rep1.max(), df.rep2.max()) + 1
ax[0].plot([lo, hi], [lo, hi], 'k--', lw=0.8); ax[0].scatter(df.rep1, df.rep2, s=22, color=BLUE)
for _, r in df[diff.abs() > 5].iterrows(): ax[0].annotate(short(r.ligand), (r.rep1, r.rep2), fontsize=6, xytext=(3, 3), textcoords='offset points')
ax[0].set_xlabel('Replica 1 \u0394TOTAL (kcal/mol)'); ax[0].set_ylabel('Replica 2 \u0394TOTAL (kcal/mol)')
ax[0].text(0.04, 0.95, '\u03c1 = %.2f' % spearmanr(df.rep1, df.rep2)[0], transform=ax[0].transAxes, va='top', fontsize=8)
ax[1].hist(diff, bins=np.arange(-10, 10.1, 2), color=BLUE, edgecolor='w')
ax[1].set_xlabel('Replica 1 \u2212 replica 2 (kcal/mol)'); ax[1].set_ylabel('Number of ligands')
ax[1].text(0.04, 0.95, 'mean |difference| = %.1f kcal/mol' % diff.abs().mean(), transform=ax[1].transAxes, va='top', fontsize=8)
for a, l in zip(ax, 'ab'): a.set_title(l, loc='left', fontweight='bold')
save('Figure4_replicas')

# Figure 5: ligand size
fig, ax = plt.subplots(figsize=(3.6, 3.4))
ax.scatter(df.heavy_atoms, df.dG_exp_kcal, s=22, color=BLUE)
ax.set_xlabel('Heavy-atom count'); ax.set_ylabel('Experimental \u0394G (kcal/mol)')
ax.text(0.04, 0.05, '\u03c1 = %.2f' % spearmanr(df.heavy_atoms, df.dG_exp_kcal)[0], transform=ax.transAxes, fontsize=8)
save('Figure5_size')

# Figure S1: Vina pose RMSD
pr = pd.read_csv(f'{A}/vina_pose_rmsd.csv').dropna(subset=['rmsd'])
w = pr.pivot(index='ligand', columns='seed', values='rmsd'); w['mean'] = w[[1, 2, 3]].mean(axis=1); w = w.sort_values('mean')
fig, ax = plt.subplots(figsize=(7.2, 3.2)); x = np.arange(len(w))
for s, m in zip((1, 2, 3), ('o', 's', '^')): ax.scatter(x, w[s], marker=m, s=18, label=f'seed {s}')
ax.axhline(2.0, ls='--', color='k', lw=0.8); ax.set_xticks(x); ax.set_xticklabels([short(l) for l in w.index], rotation=90, fontsize=7)
ax.set_ylabel('Heavy-atom RMSD to benchmark pose (\u00c5)'); ax.legend(fontsize=7, frameon=False)
save('FigureS1_pose_rmsd')

# Figure S2: ligand RMSD versus time for three ligands
def read_xvg(p):
    a = np.array([l.split() for l in open(p) if l.strip() and l[0] not in '#@'], dtype=float); return a[:, 0], a[:, 1]
fig, ax = plt.subplots(1, 3, figsize=(7.5, 2.8), sharey=True)
for a, l in zip(ax, ['2bb', '2x', '2i']):
    for r, c in ((1, BLUE), (2, RED)):
        tt, yy = read_xvg(f'{M}/lig_p38a_{l}/rep{r}/rmsd_lig3.xvg'); a.plot(tt, yy, color=c, lw=0.8, label=f'replica {r}')
    a.set_title(l); a.set_xlabel('Time (ns)')
ax[0].set_ylabel('Ligand RMSD (nm)'); ax[0].legend(fontsize=7, frameon=False)
save('FigureS2_ligand_rmsd')
print('saved figures to', F)
