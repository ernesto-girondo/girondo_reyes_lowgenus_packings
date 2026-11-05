"""Section 6.5 of the paper (Proposition 6.10): the canonical double covers of the non-orientable extremal surfaces
of genus 4.

Usage:  python3 genus4_covers.py
Output: groups/genus4_covers.txt
Time:   about 10 minutes.

The 144 non-orientable 1-extremal surfaces X of genus 4 are the NEC classes of groups/nec_N18.jsonl (torsion-free
subgroups of index 36 of the extended triangle group Delta^pm(2,3,18), the external test of Nakamura 2009). Their
canonical double covers X+ = D/K+, K+ = K cap Delta, are 2-extremal Riemann surfaces of genus 3. For each X:

  * the extremal packings of X and of X+ are computed as for the surfaces of the paper (Surface.analyse);
  * an extremal 2-packing of X+, given by alpha (K+ < alpha Delta^pm alpha^-1), descends to X if and only if it is
    invariant under the covering involution tau, i.e. if alpha^-1 t alpha lies in Delta^pm for an orientation
    reversing element t of K (a lift of tau); this is decided by the reduction algorithm of Step 4 (Geometry.reduce);
  * for each packing of X+ the group K_j = alpha^-1 K+ alpha is identified by its Delta^pm-conjugacy class (the
    canonical form `key` of the data files, abbreviated as the class_id of surfaces_summary.csv): K_j in the class of
    K+ means that the packing is the image of the original one under an isometry of X+, another class means that X+
    is also uniformized by that class (Lemma 4.2);
  * when X+ has more than one packing, the orbits of the isometry group of X+ on its packings are computed from the
    explicit isometries (autfix.all_automorphisms).

Finally the double covers are grouped into isometry classes (same Delta^pm-class of K+, or linked through the groups
K_j). The summary at the end gives the statements of Proposition 6.10.
"""
import hashlib
import json
from geom import Geometry, mmul, minv, mdet, hdist
import mpmath as mp
from surface import Surface
from batch import double_cover_action, canonical_key
from gapiface import coset_action_from_words, reduce_word
from autfix import all_automorphisms


def class_id(S):
    return hashlib.sha1(json.dumps(canonical_key(S.perm, S.d)).encode()).hexdigest()[:10]


def isometry_orbits(G, S):
    """orbits of Isom(D/K) on the extremal packings of S (lists of indices)"""
    P = S.packings
    n = len(P)
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for M, _, _ in all_automorphisms(S):
        for i in range(n):
            for j in range(n):
                if root(i) != root(j) and G.reduce(mmul(mmul(minv(P[j]['alpha']), M), P[i]['alpha'])) is not None:
                    parent[root(i)] = root(j)
    orb = {}
    for i in range(n):
        orb.setdefault(root(i), []).append(i + 1)
    return sorted(orb.values())


def main():
    G = Geometry(18)
    groups = [json.loads(l) for l in open('groups/nec_N18.jsonl') if l.strip()]
    rows = []
    for n, g in enumerate(groups, 1):
        action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
        S = Surface(G, action, label='nec_N18#%d' % n)
        S.analyse(verbose=False)
        if S.orientable:
            continue                                   # the 8 Fuchsian classes of the file are not part of this test
        # the non-orientable surfaces are pairwise non-isometric: every K_j of X is in the class of K (Lemma 4.2)
        selfcls = class_id(S)
        x_distinct = all(class_id(Surface(G, coset_action_from_words(18, [reduce_word(w) for w in P['words']]),
                                          label='K_j')) == selfcls for P in S.packings)
        S2 = Surface(G, double_cover_action(action), label='nec_N18#%d+' % n)
        S2.analyse(verbose=False)
        t = [s for s in S.sp if mdet(s.Mn) < 0][0].Mn          # an orientation reversing element of K
        packs = []
        for P in S2.packings:
            lift = G.reduce(mmul(mmul(minv(P['alpha']), t), P['alpha'])) is not None
            Kj = Surface(G, coset_action_from_words(18, [reduce_word(w) for w in P['words']]), label='K_j')
            packs.append((lift, class_id(Kj)))
        packs.sort(key=lambda p: (not p[0], p[1]))             # lifts first
        orbits = isometry_orbits(G, S2) if len(S2.packings) > 1 else [[1]]
        ok = all(c for _, c in S.rejected) and all(c for _, c in S2.rejected) and not S.n_uncertified and not S2.n_uncertified
        dX = sorted(float(hdist(mp.mpc(0, 1), P['z'])) for P in S.packings)
        rows.append({'n': n, 'dX': dX, 'X': len(S.packings), 'Xplus': len(S2.packings), 'cls': class_id(S2), 'packs': packs,
                     'orbits': orbits, 'ok': ok, 'x_distinct': x_distinct})
        print('nec_N18#%d done' % n, flush=True)
    # isometry classes of the double covers
    parent = {}

    def root(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            a = parent[a]
        return a
    for r in rows:
        for _, c in r['packs']:
            parent[root(c)] = root(r['cls'])
    surf = {}
    for r in rows:
        surf.setdefault(root(r['cls']), []).append(r['n'])
    out = ['Canonical double covers of the 144 non-orientable 1-extremal surfaces of genus 4 (nec_N18, Nakamura 2009):',
           'Section 6.5 and Proposition 6.10 of the paper. Columns: number of extremal discs of X; number of extremal',
           '2-packings of X+; class of K+ (class_id); for each packing of X+, "L" if it is a lift of a packing of X',
           '(invariant under the covering involution) or "-" if not, with the class of K_j; orbits of Isom(X+) on the',
           'packings of X+ (packings numbered as listed).', '']
    for r in rows:
        out.append('nec_N18#%-4d X: %d  X+: %d  K+: %s  packings: %s  orbits: %s%s' % (
            r['n'], r['X'], r['Xplus'], r['cls'], ' '.join('%s:%s' % ('L' if l else '-', c) for l, c in r['packs']),
            r['orbits'], '' if r['ok'] else '  CHECK FAILED'))
    out.append('')
    from collections import Counter
    cnt = Counter((r['X'], r['Xplus']) for r in rows)
    out.append('Summary')
    out.append('  surfaces: %d; all packings certified and all rejections confirmed: %s' % (len(rows), all(r['ok'] for r in rows)))
    for (a, b), m in sorted(cnt.items()):
        out.append('  X with %d extremal disc(s), X+ with %d extremal 2-packing(s): %d surfaces' % (a, b, m))
    exc = [r for r in rows if r['Xplus'] != r['X']]
    out.append('  surfaces whose double cover has more packings: %s' % ', '.join('nec_N18#%d' % r['n'] for r in exc))
    non_lift_everywhere_else = [r['n'] for r in rows if r['Xplus'] == r['X'] and not all(l for l, _ in r['packs'])]
    out.append('  packings of X+ that are not lifts, outside these surfaces: %s' % (non_lift_everywhere_else or 'none'))
    for r in exc:
        others = sorted(set(c for l, c in r['packs'] if not l) - {r['cls']})
        if others:
            mates = [s['n'] for s in rows if s['cls'] in others]
            how = 'the packings that are not lifts have K_j in the class of the double cover of %s' % ', '.join(
                'nec_N18#%d' % m for m in mates)
        else:
            how = 'the packings that are not lifts have K_j in the class of K+ (images of lifts under isometries of X+)'
        out.append('  nec_N18#%d: centres of the extremal discs of X at distances %s from O (r_c(18) = %.6f: a vertex '
                   'of the 18-gon); orbits of Isom(X+) %s; %s' % (r['n'], ', '.join('%.6f' % x for x in r['dX']),
                                                                  float(G.F.H), r['orbits'], how))
    out.append('  the %d non-orientable surfaces are pairwise non-isometric (each K_j of X in the class of K): %s' % (
        len(rows), all(r['x_distinct'] for r in rows)))
    out.append('  double covers: %d Delta^pm-classes, %d isometry classes of surfaces' % (
        len(set(r['cls'] for r in rows)), len(surf)))
    out.append('  surfaces of genus 3 that are double covers of two non-orientable surfaces: %s' % '; '.join(
        ' and '.join('nec_N18#%d' % m for m in v) for v in surf.values() if len(v) > 1))
    open('groups/genus4_covers.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out[-12:]))


if __name__ == '__main__':
    main()
