import os, re, csv
from scipy.stats import pearsonr, spearmanr
home = os.path.expanduser('~')
runs = [l.strip() for l in open(f'{home}/benchmark_project/analysis/run_list.txt') if l.strip()]
rows = []
for d in runs:
    parts = d.rstrip('/').split('/')
    if parts[-1].startswith('rep'):
        lig, rep = parts[-2], int(parts[-1][3:])
    else:
        lig, rep = parts[-1], 1
    r = dict(ligand=lig, rep=rep, frames='', VDW='', EEL='', EGB='', ESURF='', TOTAL='', SD='', blockSEM='', flags=[])
    res = os.path.join(d, 'final_results.dat')
    if not (os.path.exists(res) and os.path.exists(os.path.join(d, 'mmpbsa.ok'))):
        r['flags'].append('no result'); rows.append(r); continue
    txt = open(res, encoding='utf-8', errors='replace').read()
    if re.search(r'nan|\*{5,}', txt, re.I): r['flags'].append('NaN/stars')
    outp = os.path.join(d, 'mmpbsa.out')
    log = open(outp, 'rb').read().decode('utf-8', 'replace').replace('\r', '\n') if os.path.exists(outp) else ''
    m = re.search(r'Complex completed: (\d+) frames', log)
    r['frames'] = int(m.group(1)) if m else ''
    if r['frames'] != 251: r['flags'].append('frames=%s' % r['frames'])
    i = txt.find('Delta (Complex')
    sec = txt[i:] if i >= 0 else ''
    vals = {}
    for key in ('BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL'):
        mm = re.search(r'^\u0394' + key + r'\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)', sec, re.M)
        if mm: vals[key] = [float(x) for x in mm.groups()]
        else: r['flags'].append('missing ' + key)
    if all(k in vals for k in ('BOND', 'ANGLE', 'DIHED')):
        if max(abs(vals[k][0]) for k in ('BOND', 'ANGLE', 'DIHED')) > 0.01: r['flags'].append('internal terms not zero')
    if 'TOTAL' in vals:
        r['TOTAL'], r['SD'], r['blockSEM'] = vals['TOTAL'][0], vals['TOTAL'][2], vals['TOTAL'][6]
        for k, c in (('VDWAALS', 'VDW'), ('EEL', 'EEL'), ('EGB', 'EGB'), ('ESURF', 'ESURF')):
            if k in vals: r[c] = vals[k][0]
    rows.append(r)
out = f'{home}/benchmark_project/analysis/mmpbsa_table.csv'
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['ligand', 'rep', 'frames', 'dVDW', 'dEEL', 'dEGB', 'dESURF', 'dTOTAL', 'SD', 'blockSEM', 'flags'])
    for r in rows:
        w.writerow([r['ligand'], r['rep'], r['frames'], r['VDW'], r['EEL'], r['EGB'], r['ESURF'], r['TOTAL'], r['SD'], r['blockSEM'], ';'.join(r['flags'])])
bad = [r for r in rows if r['flags']]
print('runs read: %d | with problems: %d' % (len(rows), len(bad)))
for r in bad: print('  ', r['ligand'], 'rep', r['rep'], '->', ';'.join(r['flags']))
by = {}
for r in rows:
    if r['TOTAL'] != '': by.setdefault(r['ligand'], {})[r['rep']] = r['TOTAL']
pair = [(l, v[1], v[2]) for l, v in by.items() if 1 in v and 2 in v]
print('\nligands with both replicas: %d' % len(pair))
diffs = sorted(pair, key=lambda t: -abs(t[1] - t[2]))
print('%-16s %9s %9s %8s' % ('ligand', 'rep1', 'rep2', 'diff'))
for l, a, b in diffs: print('%-16s %9.2f %9.2f %8.2f' % (l, a, b, a - b))
if len(pair) > 2:
    a = [p[1] for p in pair]; b = [p[2] for p in pair]
    print('\nmean |rep1-rep2| = %.2f kcal/mol; max = %.2f' % (sum(abs(x - y) for x, y in zip(a, b)) / len(a), max(abs(x - y) for x, y in zip(a, b))))
    print('replica agreement across ligands: Pearson %.2f, Spearman %.2f' % (pearsonr(a, b)[0], spearmanr(a, b)[0]))
print('\ntable saved to', out)
