"""Numerical safety margins of the search (Steps 3-4 of the paper).

Every decision of the search that uses floating point numbers (40 significant digits) has the form
"is the displacement rho(z, s z) equal to an element of the finite list D_N or not?". This script re-runs the
search on all the surfaces and double covers of the paper and records, over all the decisions taken:
  * accepted equalities: the largest |rho(z,sz) - d| that was treated as an equality (should be ~ rounding error);
  * rejections: for each rejected point, the gap |rho(z,sz) - D_N| that caused the rejection; we report the
    smallest such gap over all rejected points (the safety margin of the rejections);
  * the same for the confirmation of the candidates rejected in Step 4: an element of K other than the identity
    (a side pairing or a product of two side pairings) that displaces the candidate by a non-admissible distance;
    when its displacement is beyond the list of admissible distances of the search, a longer (complete) list is
    used, as in Surface.confirm_not_centre.
If the first number is many orders of magnitude below the precision threshold (1e-8) and the second many orders
above it, no rounding error can have changed any decision.
Output: groups/margins.txt. About 6 minutes.
"""
import json
import mpmath as mp
from geom import Geometry, hdist, act, mmul, minv
from surface import Surface, TOL
from batch import double_cover_action
from outfile import tee

FAMILIES = ['nec_N12', 'nec_N9', 'nec_N8', 'nec_N7', 'fuchs_N18', 'fuchs_N12', 'fuchs_N10', 'fuchs_N9',
            'fuchs_N8', 'fuchs_N7']


def gap(x, dists):
    return min(abs(x - d) for d in dists)


tee('groups/margins.txt')
worst_eq, worst_rej, worst_conf = mp.mpf(0), mp.mpf(10), mp.mpf(10)
where_rej = where_conf = where_eq = None
n_points = n_rej = n_conf = n_single = n_ext = 0
for fam in FAMILIES:
    groups = [json.loads(l) for l in open('groups/%s.jsonl' % fam) if l.strip()]
    G = Geometry(groups[0]['N'])
    for i, g in enumerate(groups):
        action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
        acts = [(action, '%s#%d' % (fam, i + 1))]
        if fam.startswith('nec'):
            acts.append((double_cover_action(action), '%s#%d+' % (fam, i + 1)))
        for act_, label in acts:
            S = Surface(G, act_, label=label)
            S.analyse(verbose=False)
            dists = S.dists
            top = dists[-1] - TOL
            good = set(id(z) for z in S.candidates)
            for z in S.intersection_points:
                n_points += 1
                gaps = [gap(hdist(z, act(s.Mn, z)), dists) for s in S.sp if hdist(z, act(s.Mn, z)) <= top]
                if id(z) in good:
                    m = max(gaps) if gaps else mp.mpf(0)
                    if m > worst_eq:
                        worst_eq, where_eq = m, label
                else:
                    n_rej += 1
                    m = max(gaps)                    # the decisive (largest) gap of this point
                    if m < worst_rej:
                        worst_rej, where_rej = m, label
            for z, _ in S.rejected:                  # candidates without alpha: independent confirmation
                n_conf += 1
                r = S.confirm_not_centre(z, dists)   # the element of K with the largest gap (identity excluded)
                best = r[4] if r is not None else mp.mpf(0)
                if r is not None and r[1] is None:
                    n_single += 1
                if r is not None and r[5] > dists[-1]:
                    n_ext += 1
                if best < worst_conf:
                    worst_conf, where_conf = best, label
    print('%-10s done' % fam, flush=True)
print()
print('points examined in Step 3: %d (of which rejected by the displacement test: %d)' % (n_points, n_rej))
print('largest deviation treated as an equality (accepted candidates): %s  (%s)' % (mp.nstr(worst_eq, 3), where_eq))
print('smallest gap that caused a rejection in Step 3:                  %s  (%s)' % (mp.nstr(worst_rej, 3), where_rej))
print('candidates rejected in Step 4: %d; smallest gap of their confirmation by an element of K (a side pairing or a '
      'product of two): %s  (%s)' % (n_conf, mp.nstr(worst_conf, 3) if n_conf else '-', where_conf))
print('   confirmed by a product of two side pairings: %d, by a single side pairing: %d; with a list of admissible '
      'distances longer than the one of the search: %d' % (n_conf - n_single, n_single, n_ext))
print('threshold used in the program: %s; working precision: %d significant digits' % (mp.nstr(TOL, 1), mp.mp.dps))
