import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

A = os.path.expanduser('~/benchmark_project/analysis'); F = f'{A}/figures'
df = pd.read_csv(f'{A}/final_table_d2.csv'); t = pd.read_csv(f'{A}/tables/table1.csv')
short = lambda s: s.replace('lig_p38a_', '')
BLUE = '#1f4e79'
def save(fig, name):
    fig.savefig(f'{F}/{name}.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{F}/{name}.pdf', bbox_inches='tight'); plt.close(fig)

# Figure 2: forest plot
col = {'Primary': '#1f4e79', 'Secondary': '#5b8fb9', 'Sensitivity': '#c0504d', 'Exploratory': '#7f7f7f'}
fig, ax = plt.subplots(figsize=(7.5, 4.3))
ypos = np.arange(len(t))[::-1]
for yy, r in zip(ypos, t.itertuples()):
    ax.errorbar(r.spearman, yy, xerr=[[r.spearman - r.sp_lo], [r.sp_hi - r.spearman]], fmt='o', ms=4, color=col[r.group], capsize=3)
ax.axvline(0, color='k', lw=0.8)
ax.set_yticks(ypos); ax.set_yticklabels([f'{r.predictor} (N={r.N})' for r in t.itertuples()], fontsize=7.5)
ax.set_xlim(-0.5, 1.0)
ax.set_xlabel('Spearman \u03c1 with experimental \u0394G (95% bootstrap CI)', fontsize=9)
for g, c in col.items(): ax.plot([], [], 'o', color=c, label=g)
ax.legend(fontsize=7.5, loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=4, frameon=False)
save(fig, 'Figure2_forest')

# Figure 3: scores versus experiment (statistics in the panel titles, away from the data)
lab = lambda r: 'Spearman \u03c1 = %.2f (%.2f, %.2f)' % (r.spearman, r.sp_lo, r.sp_hi)
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.5))
ax[0].scatter(df.dG_exp_kcal, df.vina_mean, s=22, color=BLUE)
ax[0].set_ylabel('Vina score (kcal/mol)')
ax[1].errorbar(df.dG_exp_kcal, df.mm_avg, yerr=(df.rep1 - df.rep2).abs() / 2, fmt='o', ms=4, color=BLUE, ecolor='#9db7d1', capsize=2)
ax[1].set_ylabel('MM/GBSA \u0394TOTAL, replica mean (kcal/mol)', fontsize=8.5)
for a, l, k in zip(ax, 'ab', (0, 1)):
    a.set_xlabel('Experimental \u0394G from IC$_{50}$ (kcal/mol)')
    a.set_title(l, loc='left', fontweight='bold'); a.set_title(lab(t.iloc[k]), loc='right', fontsize=8)
fig.tight_layout(); save(fig, 'Figure3_scores_vs_experiment')

# Figure 4: replica agreement
diff = df.rep1 - df.rep2
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.5))
lo, hi = min(df.rep1.min(), df.rep2.min()) - 1, max(df.rep1.max(), df.rep2.max()) + 1
ax[0].plot([lo, hi], [lo, hi], 'k--', lw=0.8); ax[0].scatter(df.rep1, df.rep2, s=22, color=BLUE)
off = {'2g': (-13, -9), '2f': (5, 3), '2r': (-13, 3), '2m': (5, -9)}
for _, r in df[diff.abs() > 5].iterrows():
    ax[0].annotate(short(r.ligand), (r.rep1, r.rep2), fontsize=6.5, xytext=off.get(short(r.ligand), (4, 4)), textcoords='offset points')
ax[0].set_xlabel('Replica 1 \u0394TOTAL (kcal/mol)'); ax[0].set_ylabel('Replica 2 \u0394TOTAL (kcal/mol)')
ax[0].set_title('a', loc='left', fontweight='bold'); ax[0].set_title('Spearman \u03c1 = %.2f' % spearmanr(df.rep1, df.rep2)[0], loc='right', fontsize=8)
ax[1].hist(diff, bins=np.arange(-10, 10.1, 2), color=BLUE, edgecolor='w'); ax[1].set_ylim(0, 9)
ax[1].set_xlabel('Replica 1 \u2212 replica 2 (kcal/mol)'); ax[1].set_ylabel('Number of ligands')
ax[1].set_title('b', loc='left', fontweight='bold'); ax[1].set_title('mean |difference| = %.1f kcal/mol' % diff.abs().mean(), loc='right', fontsize=8)
fig.tight_layout(); save(fig, 'Figure4_replicas')
print('saved Figures 2, 3 and 4 in', F)
