"""Run the full analysis over the groups enumerated by GAP.  Usage: python3 batch.py groups/nec_N9.jsonl ..."""
import sys, json, time, traceback
import mpmath as mp
from collections import deque
from geom import Geometry, to_disc, hdist, act, mdet
from surface import Surface
from field import pari

mp.mp.dps = 40


def canonical_key(perm, d):
    """canonical form of the transitive <a,b,c>-set (isomorphism class)"""
    best = None
    for s in range(d):
        lab = {s: 0}
        order = [s]
        q = deque([s])
        while q:
            p = q.popleft()
            for ch in 'abc':
                r = perm[ch][p]
                if r not in lab:
                    lab[r] = len(order)
                    order.append(r)
                    q.append(r)
        key = tuple(tuple(lab[perm[ch][p]] for p in order) for ch in 'abc')
        if best is None or key < best:
            best = key
    return best


def parity_and_words(perm, d):
    par = [None] * d
    word = [None] * d
    par[0] = 0
    word[0] = ''
    q = deque([0])
    while q:
        p = q.popleft()
        for ch in 'abc':
            r = perm[ch][p]
            if par[r] is None:
                par[r] = 1 - par[p]
                word[r] = word[p] + ch
                q.append(r)
    return par, word


def automorphisms(S):
    """points j whose stabiliser is K (all side pairings fix j) -> automorphisms of D/K induced by Delta^pm"""
    d = S.d
    auts = []
    par, word = parity_and_words(S.perm, d)
    for j in range(d):
        if all(S.apply(j, s.word) == j for s in S.sp):
            # phi_j(0^w) = j^w
            phi = [None] * d
            for p in range(d):
                phi[p] = S.apply(j, word[p])
            # order
            o = 1
            cur = phi[:]
            while cur != list(range(d)):
                cur = [phi[x] for x in cur]
                o += 1
                if o > 10000:
                    break
            fixes = any(S.perm[ch][p] == phi[p] for p in range(d) for ch in 'abc')
            auts.append({'j': j, 'parity': par[j], 'order': o, 'has_fixed_points': fixes, 'word': word[j]})
    return auts


def double_cover_action(action):
    """action of Delta^pm on cosets of K^+ = K cap Delta, from the action on cosets of K"""
    d = action['d']
    out = {'d': 2 * d}
    for ch in 'abc':
        p = action[ch]
        out[ch] = [None] * (2 * d)
        for j in range(d):
            for e in (0, 1):
                out[ch][2 * j + e] = 2 * p[j] + (1 - e)
    return out


def analyse_one(G, action, label):
    t = time.time()
    S = Surface(G, action, label=label)
    S.analyse(verbose=False)
    auts = automorphisms(S)
    i0 = mp.mpc(0, 1)
    rec = {
        'label': label, 'N': S.N, 'd': S.d, 'k': S.k, 'orientable': S.orientable,
        'n_sidepairings': len(S.sp), 'sidepairings': [s.word for s in S.sp],
        'n_candidates': len(S.candidates),
        'n_uncertified_alphas': S.n_uncertified,
        'n_rejected_candidates': len(S.rejected),
        'rejected_confirmed_independently': all(c for _, c in S.rejected),
        'n_packings': len(S.packings),
        'packings': [{
            'centre_disc': [mp.nstr(to_disc(P['z']).real, 20), mp.nstr(to_disc(P['z']).imag, 20)],
            'dist_to_origin': mp.nstr(hdist(i0, P['z']), 20),
            'alpha_det_sign': P['det_sign'],
            'alpha_exact': str(P['alpha_exact']),
            'words': P['words'],
        } for P in S.packings],
        'automorphisms': auts,
        'n_aut': len(auts),
        'anticonformal_aut': any(a['parity'] == 1 for a in auts),
        'key': canonical_key(S.perm, S.d),
        'time': round(time.time() - t, 2),
    }
    return rec, S


if __name__ == '__main__':
    files = sys.argv[1:]
    for fn in files:
        out = fn.replace('.jsonl', '_results.json')
        logf = open(fn.replace('.jsonl', '_run.log'), 'w')
        groups = [json.loads(l) for l in open(fn) if l.strip()]
        N = groups[0]['N']
        G = Geometry(N)
        results = []
        for idx, g in enumerate(groups):
            action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
            label = '%s#%d' % (fn.split('/')[-1].replace('.jsonl', ''), idx + 1)
            try:
                rec, S = analyse_one(G, action, label)
                rec['gens_gap'] = g.get('gens')
                if not rec['orientable']:
                    rec2, S2 = analyse_one(G, double_cover_action(action), label + '+')
                    rec['double_cover'] = rec2
                msg = '%s: k=%d ori=%s packings=%d aut=%d unc=%d rej=%d/%s%s' % (
                    label, rec['k'], rec['orientable'], rec['n_packings'], rec['n_aut'], rec['n_uncertified_alphas'], rec['n_rejected_candidates'], rec['rejected_confirmed_independently'],
                    ('  cover: packings=%d aut=%d' % (rec['double_cover']['n_packings'], rec['double_cover']['n_aut'])) if 'double_cover' in rec else '')
            except Exception as e:
                rec = {'label': label, 'error': traceback.format_exc()}
                msg = '%s: ERROR %s' % (label, e)
            results.append(rec)
            logf.write(msg + '\n')
            logf.flush()
            json.dump(results, open(out, 'w'))
        logf.write('DONE\n')
        logf.close()
