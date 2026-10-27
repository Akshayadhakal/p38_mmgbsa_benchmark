import re, os
txt = open(os.path.expanduser('~/benchmark_project/analysis/check3_all.txt')).read()
flag = []
for blk in txt.split('== ')[1:]:
    name = blk.splitlines()[0].replace(os.path.expanduser('~'), '~')
    vals = {}
    for key in ('backbone RMSD', 'ligand RMSD', 'min distance'):
        m = re.search(key + r'.*?: first ([\d.]+) mean ([\d.]+) max ([\d.]+) last ([\d.]+) \((\d+) frames\)', blk)
        if m: vals[key] = tuple(float(x) for x in m.groups())
    if len(vals) < 3: flag.append((name, 'MISSING VALUES')); continue
    why = []
    if vals['ligand RMSD'][1] > 0.30: why.append('ligand RMSD mean %.2f' % vals['ligand RMSD'][1])
    if vals['min distance'][2] > 0.50: why.append('min distance max %.2f' % vals['min distance'][2])
    if vals['backbone RMSD'][1] > 0.35: why.append('backbone RMSD mean %.2f' % vals['backbone RMSD'][1])
    if int(vals['ligand RMSD'][4]) != 1001: why.append('frames %d' % vals['ligand RMSD'][4])
    if why: flag.append((name, '; '.join(why)))
print(len(txt.split('== ')) - 1, 'runs read;', len(flag), 'flagged')
for n, w in flag: print(n, '->', w)
