import sys, math
top, pdb = sys.argv[1], sys.argv[2]
atoms, bonds, sec = set(), [], None
for l in open(top):
    s = l.split(';')[0].strip()
    if not s or s.startswith('#'): continue
    if s.startswith('['): sec = s.strip('[] ').strip(); continue
    p = s.split()
    if sec == 'atoms' and len(p) >= 5 and p[0].isdigit(): atoms.add(int(p[0]))
    elif sec == 'bonds' and len(p) >= 3 and p[0].isdigit(): bonds.append((int(p[0]), int(p[1])))
models, cur = [], []
for l in open(pdb):
    if l.startswith(('ATOM', 'HETATM')):
        cur.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    elif l.startswith('ENDMDL'):
        models.append(cur); cur = []
worst, bad = 0.0, 0
for xyz in models:
    m = max(math.dist(xyz[a-1], xyz[b-1]) for a, b in bonds if a <= len(xyz) and b <= len(xyz))
    worst = max(worst, m)
    if m > 2.5: bad += 1
print('frames scanned: %d | frames with a bond > 2.5 A: %d | longest bond over all frames: %.2f A' % (len(models), bad, worst))
