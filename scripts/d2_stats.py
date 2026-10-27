import os, re
import numpy as np, pandas as pd
from scipy.stats import spearmanr
home = os.path.expanduser('~'); A = f'{home}/benchmark_project/analysis'
runs = [l.strip() for l in open(f'{A}/run_list.txt') if l.strip()]
rows, problems = [], []
for d in runs:
    parts = d.rstrip('/').split('/')
    lig, rep = (parts[-2], int(parts[-1][3:])) if parts[-1].startswith('rep') else (parts[-1], 1)
    res, outp = f'{d}/final_results_d2.dat', f'{d}/mmpbsa_d2.out'
    if not (os.path.exists(f'{d}/mmpbsa_d2.ok') and os.path.exists(res)):
        problems.append(f'{lig} rep{rep}: no d2 result'); continue
    txt = open(res, encoding='utf-8', errors='replace').read()
    if re.search(r'nan|\*{5,}', txt, re.I): problems.append(f'{lig} rep{rep}: NaN/stars'); continue
    log = open(outp, 'rb').read().decode('utf-8', 'replace').replace('\r', '\n')
    m = re.search(r'Complex completed: (\d+) frames', log)
    if not m or int(m.group(1)) != 251: problems.append(f'{lig} rep{rep}: frames {m.group(1) if m else None}'); continue
    sec = txt[txt.find('Delta (Complex'):]
    internal = [re.search(r'^\u0394' + k + r'\s+(-?[\d.]+)', sec, re.M) for k in ('BOND', 'ANGLE', 'DIHED')]
    if any(x is None or abs(float(x.group(1))) > 0.01 for x in internal): problems.append(f'{lig} rep{rep}: internal terms'); continue
    t = re.search(r'^\u0394TOTAL\s+(-?[\d.]+)', sec, re.M)
    if not t: problems.append(f'{lig} rep{rep}: no TOTAL'); continue
    rows.append((lig, rep, float(t.group(1))))
print('runs read: %d | with problems: %d' % (len(rows) + len(problems), len(problems)))
for p in problems: print('  ', p)
d2 = pd.DataFrame(rows, columns=['ligand', 'rep', 'dTOTAL_d2']).pivot(index='ligand', columns='rep', values='dTOTAL_d2')
d2.columns = ['d2_rep1', 'd2_rep2']; d2['d2_avg'] = d2.mean(axis=1)
df = pd.read_csv(f'{A}/final_table.csv').merge(d2.reset_index(), on='ligand')
print('ligands with both replicas:', len(df))
rng = np.random.default_rng(12345)
def ci(x, y, n=5000):
    out = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.unique(x[i]).size > 1 and np.unique(y[i]).size > 1: out.append(spearmanr(x[i], y[i])[0])
    return np.percentile(out, [2.5, 97.5])
y = df.dG_exp_kcal.values
print('\nSPEARMAN vs experimental dG (N=%d)' % len(df))
for lab, col in (('Vina', 'vina_mean'), ('MM/GBSA, dielectric 1 (primary)', 'mm_avg'), ('MM/GBSA, dielectric 2 (sensitivity)', 'd2_avg'), ('  dielectric 2, replica 1 only', 'd2_rep1'), ('  dielectric 2, replica 2 only', 'd2_rep2')):
    x = df[col].values; lo, hi = ci(x, y)
    print('%-40s %5.2f [%5.2f, %5.2f]' % (lab, spearmanr(x, y)[0], lo, hi))
print('\nreplica agreement, dielectric 2: mean |rep1-rep2| = %.2f kcal/mol (dielectric 1: %.2f)' % ((df.d2_rep1 - df.d2_rep2).abs().mean(), (df.rep1 - df.rep2).abs().mean()))
print('ranking agreement between dielectric 1 and 2 (Spearman): %.2f' % spearmanr(df.mm_avg, df.d2_avg)[0])
print('mean dTOTAL: dielectric 1 %.1f, dielectric 2 %.1f kcal/mol' % (df.mm_avg.mean(), df.d2_avg.mean()))
df.to_csv(f'{A}/final_table_d2.csv', index=False)
