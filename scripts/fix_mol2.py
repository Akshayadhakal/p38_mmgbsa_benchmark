import sys
f = sys.argv[1]
out, sec, cnt = [], None, {}
for l in open(f).read().splitlines():
    if l.startswith('@<TRIPOS>'):
        sec = l.strip()
        out.append(l)
        if sec == '@<TRIPOS>SUBSTRUCTURE':
            out.append('     1 LIG         1 TEMP              0 ****  ****    0 ROOT')
        continue
    if sec == '@<TRIPOS>SUBSTRUCTURE':
        continue
    if sec == '@<TRIPOS>ATOM' and l.strip():
        p = l.split()
        el = p[5].split('.')[0]
        cnt[el] = cnt.get(el, 0) + 1
        name = f"{el}{cnt[el]}"
        q = float(p[8]) if len(p) > 8 else 0.0
        out.append("%7s %-6s %9.4f %9.4f %9.4f %-6s %4s %-6s %9.4f" % (p[0], name, float(p[2]), float(p[3]), float(p[4]), p[5], '1', 'LIG', q))
        continue
    out.append(l)
open(f, 'w').write("\n".join(out) + "\n")
