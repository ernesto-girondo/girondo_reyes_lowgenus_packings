"""Summary tables of Section 6 of the paper and consistency checks, from groups/*_results.json.
Output: groups/report.txt (also printed). Only the families of the paper are used (not the external test nec_N18).
"""
import json
from collections import Counter

NEC = ['nec_N12', 'nec_N9', 'nec_N8', 'nec_N7']
FUCHS = ['fuchs_N18', 'fuchs_N12', 'fuchs_N10', 'fuchs_N9', 'fuchs_N8', 'fuchs_N7']
res = {name: json.load(open('groups/%s_results.json' % name)) for name in NEC + FUCHS}

lines = []
P = lines.append


def dist(counter):
    return '  '.join('%d(%d)' % (n, counter[n]) for n in sorted(counter))


P('=== 1. Robustness of the computation (Section 4) ===')
tot = 0
for name in NEC + FUCHS:
    L = res[name]
    errs = [r for r in L if 'error' in r]
    unc = sum(r.get('n_uncertified_alphas', 0) for r in L if 'error' not in r)
    rej_ok = all(r.get('rejected_confirmed_independently', True) for r in L if 'error' not in r)
    covers = [r['double_cover'] for r in L if 'double_cover' in r]
    unc_c = sum(r.get('n_uncertified_alphas', 0) for r in covers)
    rej_c = all(r.get('rejected_confirmed_independently', True) for r in covers)
    tot += len(L) + len(covers)
    P('%-10s surfaces=%3d double covers=%2d errors=%d alphas not verified=%d (covers: %d) '
      'all rejections confirmed=%s' % (name, len(L), len(covers), len(errs), unc, unc_c, rej_ok and rej_c))
P('total: %d surfaces (373 of Table 2 and 50 double covers)' % tot)

P('')
P('=== 2. Non-orientable surfaces of genus 3: number of extremal k-packings (Theorem 6.1) ===')
P('notation m(n): n surfaces have exactly m packings')
for name in ('nec_N7', 'nec_N8', 'nec_N9', 'nec_N12'):
    L = res[name]
    P('  N=%2d k=%d  surfaces=%2d  packings: %s' % (L[0]['N'], L[0]['k'], len(L), dist(Counter(r['n_packings'] for r in L))))

P('')
P('=== 3. Canonical double covers (Theorem 6.2) ===')
for name in ('nec_N7', 'nec_N8', 'nec_N9', 'nec_N12'):
    L = res[name]
    bad = [r['label'] for r in L if r['n_packings'] != r['double_cover']['n_packings']]
    P('  N=%2d: double cover with the same number of packings: %s' % (
        L[0]['N'], 'all %d' % len(L) if not bad else 'NOT for ' + ', '.join(bad)))
    twice = sum(1 for r in L if r['double_cover']['n_aut'] == 2 * r['n_aut'])
    P('        induced group of the cover = twice the induced group of X: %d of %d' % (twice, len(L)))
    for r in L:
        if r['double_cover']['n_aut'] != 2 * r['n_aut']:
            assert r['double_cover']['n_aut'] % (2 * r['n_aut']) == 0
            P('        exception: %s  |induced on X| = %d, |induced on the cover %s| = %d' % (
                r['label'], r['n_aut'], r['double_cover']['label'], r['double_cover']['n_aut']))

P('')
P('=== 4. Orientable surfaces of genus 2: conjugacy classes (Table 2) ===')
fuchs_keys = {}
for name in FUCHS:
    L = res[name]
    keys = {}
    for r in L:
        keys.setdefault(tuple(map(tuple, r['key'])), []).append(r)
    fuchs_keys[name] = keys
    nsym = sum(1 for v in keys.values() if v[0]['anticonformal_aut'])
    # a Delta^pm-class without induced anticonformal automorphism splits into two Delta-classes (K and aKa^-1)
    ok = all((len(v) == 1) == v[0]['anticonformal_aut'] for v in keys.values())
    P('  N=%2d k=%2d  classes in Delta: %3d   classes in Delta^pm: %3d   with induced anticonformal symmetry: %2d   '
      '(splitting consistent: %s)' % (L[0]['N'], L[0]['k'], len(L), len(keys), nsym, ok))

P('')
P('=== 5. Orientable surfaces of genus 2: number of extremal k-packings (Theorem 6.3) ===')
P('notation m(n): n classes have exactly m packings')
for name in FUCHS:
    keys = fuchs_keys[name]
    c = Counter(v[0]['n_packings'] for v in keys.values())
    c2 = Counter(r['n_packings'] for r in res[name])
    P('  N=%2d k=%2d  per Delta^pm-class: %s    per Delta-class: %s' % (
        res[name][0]['N'], res[name][0]['k'], dist(c), dist(c2)))
P('  (for N=9 the two classes with 5 packings, fuchs_N9#36 and #37, are the same surface X6: see identify2.py)')

P('')
P('=== 6. Fuchsian classes with induced anticonformal symmetry that are not double covers (Section 6.2) ===')
for name in ('nec_N12', 'nec_N9', 'nec_N8', 'nec_N7'):
    N = res[name][0]['N']
    L = res[name]
    cover_keys = set(tuple(map(tuple, r['double_cover']['key'])) for r in L)
    keys = fuchs_keys['fuchs_N%d' % N]
    sym = [k for k, v in keys.items() if v[0]['anticonformal_aut']]
    notcov = [k for k in sym if k not in cover_keys]
    missing = [k for k in cover_keys if k not in keys]
    P('  N=%2d: %d non-orientable surfaces -> %d distinct Fuchsian classes (all in the list: %s); '
      'classes with induced anticonformal symmetry: %d, of which double covers: %d, not double covers: %d' % (
          N, len(L), len(cover_keys), not missing, len(sym), len(sym) - len(notcov), len(notcov)))
    for k in notcov:
        r = keys[k][0]
        odd = [a for a in r['automorphisms'] if a['parity'] == 1]
        desc = Counter((a['order'], a['has_fixed_points']) for a in odd)
        P('      %s: packings=%d |induced|=%d; induced anticonformal automorphisms (order, has fixed points): %s' % (
            r['label'], r['n_packings'], r['n_aut'], dict(sorted(desc.items()))))

P('')
P('=== 7. k=1, genus 2 (N=18), surface by surface (Remark 6.7) ===')
for r in res['fuchs_N18']:
    hidden = sorted(round(float(p['dist_to_origin']), 6) for p in r['packings'] if float(p['dist_to_origin']) > 1e-9)
    P('  %s: discs=%d |induced|=%d anticonformal=%s  distances from O of the other disc centres: %s' % (
        r['label'], r['n_packings'], r['n_aut'], r['anticonformal_aut'], hidden))

open('groups/report.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
