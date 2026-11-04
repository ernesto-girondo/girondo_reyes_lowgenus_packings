"""Hidden automorphisms, real structures and identifications between classes (Theorem 6.3(i),(ii)).

For every hidden packing of every orientable surface D/K, the group K_j = alpha^-1 K alpha < Delta is identified by
its Delta-conjugacy class (delta_key): K_j ~_Delta K gives a hidden conformal automorphism, K_j ~_Delta aKa^-1 a
hidden anticonformal one, and any other class means that two different classes give the same Riemann surface.
Output: groups/surface_identifications.json (per class) and groups/surface_identifications.txt (summary lines).
About 7 minutes."""
import json, glob, sys
from collections import deque, Counter
from geom import Geometry
from surface import Surface
from gapiface import coset_action_from_words, reduce_word
from batch import parity_and_words

def delta_key(action):
    d = action['d']
    par, _ = parity_and_words({'a': action['a'], 'b': action['b'], 'c': action['c']}, d)
    x = [action['b'][action['a'][p]] for p in range(d)]   # ab : apply a then b (right action)
    y = [action['c'][action['b'][p]] for p in range(d)]
    best = None
    for s in [p for p in range(d) if par[p] == 0]:
        lab = {s: 0}; order = [s]; q = deque([s])
        while q:
            p = q.popleft()
            for perm in (x, y):
                r = perm[p]
                if r not in lab:
                    lab[r] = len(order); order.append(r); q.append(r)
        key = tuple(tuple(lab[perm[p]] for p in order) for perm in (x, y))
        if best is None or key < best: best = key
    return best

if __name__ == "__main__":
    summary = {}
    lines = []
    for fn in sorted(glob.glob('groups/fuchs_*.jsonl')):
        name = fn.split('/')[-1].replace('.jsonl', '')
        groups = [json.loads(l) for l in open(fn) if l.strip()]
        res = json.load(open(fn.replace('.jsonl', '_results.json')))
        G = Geometry(groups[0]['N'])
        actions = {}
        dkeys = {}
        for idx, g in enumerate(groups):
            label = '%s#%d' % (name, idx + 1)
            actions[label] = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
            dkeys[label] = delta_key(actions[label])
        inv = {}
        for l, k in dkeys.items(): inv.setdefault(k, []).append(l)
        assert all(len(v) == 1 for v in inv.values())
        inv = {k: v[0] for k, v in inv.items()}
        pmkey = {r['label']: tuple(map(tuple, r['key'])) for r in res}
        mirror = {}
        for l in dkeys:
            mates = [m for m in dkeys if m != l and pmkey[m] == pmkey[l]]
            mirror[l] = mates[0] if mates else None
        info = {}
        for label, action in actions.items():
            S = Surface(G, action, label=label); S.analyse(verbose=False)
            hid_conf = 0; hid_anti = 0; others = []
            for P in S.packings:
                if G.reduce(P['alpha']) is not None: continue
                a2 = coset_action_from_words(S.N, [reduce_word(w) for w in P['words']])
                k = delta_key(a2)
                tgt = inv.get(k)
                if tgt == label: hid_conf += 1
                elif tgt is not None and tgt == mirror[label]: hid_anti += 1
                elif tgt is None: others.append('UNKNOWN')
                else: others.append(tgt)
            info[label] = {'packings': len(S.packings), 'induced_aut': sum(1 for a in res[int(label.split('#')[1]) - 1]['automorphisms']),
                           'induced_anti': res[int(label.split('#')[1]) - 1]['anticonformal_aut'],
                           'hidden_conf_links': hid_conf, 'hidden_anti_links': hid_anti, 'other_surfaces': sorted(set(others)), 'mirror': mirror[label]}
        summary[name] = info
        n_real_hidden = sum(1 for l, v in info.items() if not v['induced_anti'] and v['hidden_anti_links'] > 0)
        n_nonreal = sum(1 for l, v in info.items() if not v['induced_anti'] and v['hidden_anti_links'] == 0)
        ident = sorted(set(tuple(sorted([l] + v['other_surfaces'])) for l, v in info.items() if v['other_surfaces']))
        line = ('%-10s classes in Delta=%3d | with induced anticonformal automorphism=%3d | with only hidden anticonformal '
                'automorphisms=%3d | without anticonformal automorphisms=%3d | distinct Delta-classes giving the same Riemann '
                'surface: %s' % (name, len(info), sum(1 for v in info.values() if v['induced_anti']), n_real_hidden,
                                 n_nonreal, ident))
        print(line)
        lines.append(line)
        sys.stdout.flush()
    json.dump(summary, open('groups/surface_identifications.json', 'w'), indent=1)
    open('groups/surface_identifications.txt', 'w').write('\n'.join(lines) + '\n')
