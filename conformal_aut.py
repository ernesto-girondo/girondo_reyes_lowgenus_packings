"""Full conformal automorphism group of every orientable surface of a family, including the automorphisms that are
not induced by the triangle group (Theorem 6.3(iv) for fuchs_N9, Remark 6.7 for fuchs_N18).

The automorphisms are obtained from the hidden packings as in autfix.all_automorphisms; the order of each one is
computed in Aut(D/K). The isomorphism type is read off from the order statistics (enough for the groups that occur:
Z2, Z2 x Z2, D3, D6, ...). Classes are taken up to Delta^pm-conjugacy (one representative per surface).
Usage:  python3 conformal_aut.py fuchs_N9      -> groups/conformal_aut_fuchs_N9.txt   (about 1 minute)
        python3 conformal_aut.py fuchs_N18     -> groups/conformal_aut_fuchs_N18.txt  (about 1 minute)
"""
import json, sys
from collections import Counter
from geom import Geometry
from surface import Surface
from autfix import all_automorphisms, order_of


def group_type(orders):
    """isomorphism type from the multiset of element orders, for the small groups that occur"""
    n = sum(orders.values())
    o = dict(orders)
    if n == 1:
        return '1'
    if max(o) == n:
        return 'Z%d' % n
    if n == 4:
        return 'Z2xZ2'
    if n == 6 and o.get(2) == 3:
        return 'D3'
    if n == 8 and o.get(2) == 7:
        return 'Z2^3'
    if n == 8 and o.get(2) == 3 and o.get(4) == 4:
        return 'Z2xZ4'
    if n == 8 and o.get(2) == 5:
        return 'D4'
    if n == 12 and o.get(2) == 7 and o.get(6) == 2:
        return 'D6'
    return 'order %d, element orders %s' % (n, o)


fam = sys.argv[1] if len(sys.argv) > 1 else 'fuchs_N9'
groups = [json.loads(l) for l in open('groups/%s.jsonl' % fam) if l.strip()]
res = json.load(open('groups/%s_results.json' % fam))
G = Geometry(groups[0]['N'])
lines = []
types = Counter()
seen = set()
for i, g in enumerate(groups):
    k = tuple(map(tuple, res[i]['key']))
    if k in seen:
        continue
    seen.add(k)
    S = Surface(G, {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]},
                label='%s#%d' % (fam, i + 1))
    S.analyse(verbose=False)
    conf, ind = Counter(), Counter()
    for gm, src, par in all_automorphisms(S):
        if par == 0:
            o = order_of(S, gm)
            conf[o] += 1
            if src == 'induced':
                ind[o] += 1
    t = group_type(conf)
    types[t] += 1
    line = '%s#%d packings=%d |Aut+|=%d  type %s  element orders %s  (induced by the triangle group: %s)' % (
        fam, i + 1, len(S.packings), sum(conf.values()), t, dict(sorted(conf.items())), dict(sorted(ind.items())))
    print(line, flush=True)
    lines.append(line)
lines.append('summary (one line per Delta^pm-class): %s' % dict(sorted(types.items())))
print(lines[-1])
open('groups/conformal_aut_%s.txt' % fam, 'w').write('\n'.join(lines) + '\n')
