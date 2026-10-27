import os, re, numpy as np, pandas as pd
A = os.path.expanduser('~/benchmark_project/analysis'); T = f'{A}/tables'
short = lambda s: s.replace('lig_p38a_', '')
df = pd.read_csv(f'{A}/final_table_d2.csv')
p = pd.read_csv(f'{A}/p38_ligands.csv')[['name', 'value', 'unit']].rename(columns={'name': 'ligand', 'value': 'IC50', 'unit': 'IC50_unit'})
pr = pd.read_csv(f'{A}/vina_pose_rmsd.csv')
pm = pr.groupby('ligand').rmsd.mean().rename('vina_pose_rmsd_A').reset_index()

# Table S1: per-ligand data
s1 = df.merge(p, on='ligand').merge(pm, on='ligand').sort_values('dG_exp_kcal')
s1['ligand'] = s1.ligand.map(short)
s1 = s1.rename(columns={'dG_exp_kcal': 'dG_exp_kcal_mol', 'heavy_atoms': 'heavy_atoms', 'vina_mean': 'vina_score',
    'rep1': 'MMGBSA_eps1_rep1', 'rep2': 'MMGBSA_eps1_rep2', 'mm_avg': 'MMGBSA_eps1_mean',
    'd2_rep1': 'MMGBSA_eps2_rep1', 'd2_rep2': 'MMGBSA_eps2_rep2', 'd2_avg': 'MMGBSA_eps2_mean'})
s1 = s1[['ligand', 'IC50', 'IC50_unit', 'dG_exp_kcal_mol', 'heavy_atoms', 'vina_score', 'vina_pose_rmsd_A',
         'MMGBSA_eps1_rep1', 'MMGBSA_eps1_rep2', 'MMGBSA_eps1_mean', 'MMGBSA_eps2_rep1', 'MMGBSA_eps2_rep2', 'MMGBSA_eps2_mean']]
s1.round(2).to_csv(f'{T}/TableS1_ligands.csv', index=False)

# Table S2: per-run MD metrics
txt = open(f'{A}/check3_all.txt').read(); rows = []
def g(blk, key, k):
    m = re.search(key + r'.*?mean ([\d.]+) max ([\d.]+)', blk); return float(m.group(k)) if m else np.nan
for blk in txt.split('== ')[1:]:
    parts = blk.splitlines()[0].strip().rstrip('/').split('/')
    lig, rep = (parts[-2], int(parts[-1][3:])) if parts[-1].startswith('rep') else (parts[-1], 1)
    rows.append(dict(ligand=short(lig), replica=rep,
        backbone_rmsd_mean_nm=g(blk, 'backbone RMSD', 1), backbone_rmsd_max_nm=g(blk, 'backbone RMSD', 2),
        ligand_rmsd_mean_nm=g(blk, 'ligand RMSD', 1), ligand_rmsd_max_nm=g(blk, 'ligand RMSD', 2),
        closest_distance_mean_nm=g(blk, 'min distance', 1), closest_distance_max_nm=g(blk, 'min distance', 2)))
s2 = pd.DataFrame(rows).sort_values(['ligand', 'replica']); s2.to_csv(f'{T}/TableS2_md_runs.csv', index=False)

# Table S3: Vina pose RMSD per seed
s3 = pr.pivot(index='ligand', columns='seed', values='rmsd'); s3.columns = ['seed1_A', 'seed2_A', 'seed3_A']
s3['mean_A'] = s3.mean(axis=1); s3.index = s3.index.map(short)
s3.sort_values('mean_A', ascending=False).round(2).to_csv(f'{T}/TableS3_pose_rmsd.csv')
print('Table S1: %d rows | Table S2: %d rows (expect 58) | Table S3: %d rows' % (len(s1), len(s2), len(s3)))
print(s2.isna().sum().sum(), 'missing values in S2')
