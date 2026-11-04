"""Run the search of Section 4 of the paper again for one surface, from its permutation representation.

Usage:  python3 check_surface.py fuchs_N9 21            (the surface X_3 of the paper; a few seconds)
        python3 check_surface.py fuchs_N9 21 --details  (also every ED-cycle and every point examined, see below)
        python3 check_surface.py nec_N9 2 --cover      (the canonical double cover of a non-orientable surface)

Without --details it prints the steps of the search in summary: the bound D and the number of admissible
distances, the number of intersections of the ED-cycles inside the search disc, the number of candidates that pass
the displacement test, the packings found and verified exactly (one centre of each, in the disc model), and, for
the candidates without a conjugator, whether some element of K other than the identity (a side pairing or a
product of two) displaces them by a non-admissible distance (so that they are not centres).

With --details it prints, in addition, everything needed to follow Steps 2 and 3 by hand:
  * the admissible distances d_0 < d_1 < ... up to D (Step 2);
  * the two side pairings s_1, s_2 used in Step 3 (as words in a, b, c), their translation lengths, and the
    admissible distances d_i for which their ED-cycles b_s(d_i) = {z : d(z, s z) = d_i} are used;
  * every intersection point of an ED-cycle of s_1 with one of s_2 inside the search disc (the circumscribed disc of
    the central polygon, enlarged by 0.1), numbered #1, #2, ..., with its position, the pairs (i, j) such that it
    lies on b_{s_1}(d_i) and on b_{s_2}(d_j), and what happens to it: either a side pairing moves it by a distance
    that is not admissible (it is discarded in Step 3, and the side pairing and the distance are printed), or it is
    a candidate, and then it is the centre of a disc of the packing(s) printed (Steps 4 and 5), or it has no
    conjugator and a product of two side pairings moves it by a non-admissible distance (it is not a centre).
The packings are numbered as in draw_surface.py: packing 1 is the one of the uniformization. The same report is
written by `python3 draw_surface.py <family> <n> --edcycles` next to the picture of the ED-cycles.
"""
import json
import sys
import mpmath as mp
import numpy as np
from geom import Geometry, to_disc, hdist
from surface import Surface
from batch import double_cover_action
from wordfmt import compress


def disc_point(z):
    """position of a point of the upper half plane in the disc model, with distance from O and angle"""
    w = complex(to_disc(z))
    d = float(hdist(mp.mpc(0, 1), z))
    if d < 1e-9:
        return '%+.8f %+.8fi  (this is O)' % (0.0, 0.0)
    return '%+.8f %+.8fi  (distance %.8f from O, angle %+9.4f deg)' % (w.real, w.imag, d, np.degrees(np.angle(w)))


def details_lines(G, S):
    """the detailed report of Steps 2 and 3 (S must have been analysed, packings sorted as wanted)"""
    se = S.search
    s1, s2 = se['s1'], se['s2']
    dists = S.dists
    out = ['Step 2. Admissible distances up to D = %s (%d values):' % (mp.nstr(S.D, 8), len(dists))]
    for i, d in enumerate(dists):
        out.append('   d_%-3d = %s' % (i, mp.nstr(d, 12)))
    out.append('')
    out.append('Step 3. The side pairings used (the two that move O the least among those with different axes):')
    for name, s, Dm in (('s_1', s1, se['D1']), ('s_2', s2, se['D2'])):
        used = [i for i, d in enumerate(dists) if d <= Dm]
        nonempty = [i for i in used if s.ed_cycle(dists[i])]
        out.append('   %s = %s   (%s, translation length %s, moves O by %s)' % (
            name, compress(s.word), 'glide reflection' if s.rev else 'hyperbolic', mp.nstr(s.T, 10),
            mp.nstr(s.d0, 10)))
        out.append('      ED-cycles b_%s(d_i) for d_i <= %s: i = %s (empty for i = %s, because d_i is smaller than '
                   'what %s moves any point)' % (name[-1], mp.nstr(Dm, 8), ', '.join(map(str, nonempty)) or '-',
                                                 ', '.join(str(i) for i in used if i not in nonempty) or '-', name)
                   if len(nonempty) < len(used) else
                   '      ED-cycles b_%s(d_i) for d_i <= %s: i = %s' % (name[-1], mp.nstr(Dm, 8),
                                                                       ', '.join(map(str, nonempty))))
    out.append('   search disc: centred at O, radius %s (circumradius of the polygon + 0.1)' % mp.nstr(se['radius'], 10))
    out.append('')
    pts = se['points']
    ncand = sum(1 for p in pts if p['killed_by'] is None)
    out.append('Intersection points of an ED-cycle of s_1 and one of s_2 in the search disc: %d, of which %d pass '
               'the displacement test (candidates).' % (len(pts), ncand))
    for n, p in enumerate(pts, 1):
        out.append('#%-3d %s' % (n, disc_point(p['z'])))
        out.append('      on b_1(d_i) and b_2(d_j) for (i, j) = %s' % ', '.join('(%d, %d)' % ij for ij in p['pairs']))
        if p['killed_by'] is not None:
            w, dz = p['killed_by']
            out.append('      discarded: the side pairing %s moves it by %s, which is not an admissible distance'
                       % (compress(w), mp.nstr(dz, 10)))
        elif p['packings']:
            idx = sorted(S.packings.index(P) + 1 for P in p['packings'])
            out.append('      candidate: centre of a disc of packing%s %s (verified exactly)' % (
                's' if len(idx) > 1 else '', ', '.join(map(str, idx))))
        else:
            r = p['rejection']
            if r is None:
                out.append('      candidate without a conjugator (no product of two side pairings excludes it)')
            else:
                out.append('      candidate without a conjugator: %s %s%s moves it by %s, not admissible, so it is not '
                           'a centre' % (compress(r[0]), compress(r[1]), '^-1' if r[2] else '', mp.nstr(r[3], 10)))
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        print(__doc__)
        sys.exit(1)
    family, n = args[0], int(args[1])
    groups = [json.loads(l) for l in open('groups/%s.jsonl' % family) if l.strip()]
    g = groups[n - 1]
    action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
    label = '%s#%d' % (family, n)
    if '--cover' in sys.argv:
        action = double_cover_action(action)
        label += '+'
    G = Geometry(g['N'])
    S = Surface(G, action, label=label)
    S.analyse(verbose=True)
    if S.rejected:
        ok = sum(1 for _, conf in S.rejected if conf)
        print('   candidates without a conjugator: %d; shown not to be centres by an element of K (a product of two side '
              'pairings or a side pairing): %d' % (len(S.rejected), ok))
    else:
        print('   every candidate has a conjugator (several candidates may give the same packing)')
    if '--details' in sys.argv:
        S.packings.sort(key=lambda P: G.reduce(P['alpha']) is None)   # packing 1 = the one of the uniformization
        print()
        print('\n'.join(details_lines(G, S)))


if __name__ == '__main__':
    main()
