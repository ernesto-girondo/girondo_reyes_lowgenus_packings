"""Identify a surface given by its monodromy, and find it in the catalogue of the paper.

The monodromy can be given in two ways (points are numbered 1, ..., d; a permutation is written either in cycle
notation, "(1,2)(3,4)...", as GAP prints it, or as the list of images, "[2,1,4,3,...]"):

  1. Three permutations a, b, c: the action of the three reflections of the extended triangle group
     Delta^pm(2,3,N) on the cosets of the surface group K (the sheets of the covering D/K -> D/Delta^pm), with the
     conventions of the paper: a, b, c are involutions and ab, bc, ca have orders 2, 3, N. This works for orientable
     and non-orientable surfaces, and it is the form in which the classes are stored in groups/<family>.jsonl.

         python3 identify_class.py 9 abc "<a>" "<b>" "<c>"

  2. Two permutations x, y (orientable surfaces only): the monodromy of the Belyi function D/K -> D/Delta(2,3,N),
     i.e. of the dessin of type (2,3,N), with x of order 2 and y of order 3 (in either order), whose product has
     order N. In the notation of the paper x = ab and y = bc.

         python3 identify_class.py 9 xy "<x>" "<y>"

What it does:
  * It checks that the permutations define a compact surface: they generate a transitive group, they satisfy the
    relations of the triangle group, and no point is fixed by a, b, c, ab, bc, ca (resp. x, y, xy) or a power of
    them other than the identity: in other words K has no elliptic elements and no reflections. It reports the
    index, the Euler characteristic, the orientability, the genus and k (the number of N-gons, i.e. of discs).
  * It computes the canonical form of the monodromy (field `key` of the results files) and the fingerprint
    `class_id` of the table groups/surfaces_summary.csv. They do not depend on the numbering of the points: two
    monodromies give the same canonical form exactly when they differ by a renumbering of the sheets, i.e. when the
    two groups K are conjugate in Delta^pm(2,3,N).
  * It looks for the class in the catalogue (the ten families of the paper and the external test nec_N18) and
    prints its line of groups/surfaces_summary.csv (number of packings, induced isometries, ...) and of
    groups/automorphism_groups.txt (isometry group and conformal automorphism group).

Mirror images. Two Fuchsian groups K and sigma K sigma^-1 (sigma a reflection of the triangle group) uniformize
mirror-image surfaces; they are conjugate in Delta^pm, so they have the same canonical form and the same class_id,
but they can be two different classes of the GAP enumeration (for instance fuchs_N18#6 and fuchs_N18#7). In that case
the program says which of the two the monodromy is, but this depends on the conventions: with the conventions of
the paper (permutations act on the right, products are read from left to right, x = ab, y = bc), and with the
opposite convention (left actions, or y replaced by its inverse), which gives the other one. Both are reported.

The program uses only the files in groups/ and plain Python (no PARI, no GAP). `python3 identify_class.py --selftest`
recomputes the canonical form of every stored class and checks it against the stored one, and identifies every
class from its own monodromy, after a random renumbering of its points.
"""
import csv
import hashlib
import json
import os
import random
import re
import sys
from collections import deque

HERE = os.path.dirname(os.path.abspath(__file__))
FAMILIES = ['nec_N12', 'nec_N9', 'nec_N8', 'nec_N7', 'fuchs_N18', 'fuchs_N12', 'fuchs_N10', 'fuchs_N9', 'fuchs_N8',
            'fuchs_N7', 'nec_N18']

NAMES = {   # the surfaces that the paper calls by name (table of Section 4.1 of VERIFYING.md)
    'fuchs_N18#9': 'X_6 of Girondo-Gonzalez-Diez (Theorem 6.5), as a 1-extremal surface',
    'fuchs_N9#36': 'X_6 of Theorem 6.5, as a 4-extremal surface (the group K_orb)',
    'fuchs_N9#37': 'X_6 of Theorem 6.5, as a 4-extremal surface (the group K_inv)',
    'fuchs_N18#2': 'X_2 of Girondo-Gonzalez-Diez, as a 1-extremal surface',
    'fuchs_N9#29': 'X_2 of Girondo-Gonzalez-Diez, as a 4-extremal surface',
    'fuchs_N18#3': 'X_3 of Remark 6.6, as a 1-extremal surface',
    'fuchs_N9#21': 'X_3 of Remark 6.6, as a 4-extremal surface',
    'fuchs_N18#1': 'X_1 of Girondo-Gonzalez-Diez (Remark 6.7)',
    'fuchs_N18#4': 'X_4 of Girondo-Gonzalez-Diez (Remark 6.7)',
    'fuchs_N18#8': 'X_5 of Girondo-Gonzalez-Diez (Remark 6.7)',
    'fuchs_N18#5': 'X_7 of Girondo-Gonzalez-Diez (Remark 6.7)',
    'fuchs_N18#6': 'X_8 of Girondo-Gonzalez-Diez (Remark 6.7)',
    'fuchs_N18#7': 'the mirror image of X_8 of Girondo-Gonzalez-Diez (Remark 6.7)',
    'fuchs_N8#77': "Bolza's curve y^2 = x^5 - x (Remark 6.4)",
    'fuchs_N10#18': 'the curve y^2 = x^5 - 1 (Remark 6.4)',
}


# ---------------------------------------------------------------- permutations (0-based lists of images)
def parse_perm(text, d=None):
    """'(1,2)(3,4)' or '[2,1,4,3]' (points numbered from 1) -> list of images, 0-based"""
    t = text.strip()
    if t.startswith('['):
        imgs = [int(x) - 1 for x in re.findall(r'-?\d+', t)]
        if sorted(imgs) != list(range(len(imgs))):
            raise ValueError('not a permutation of 1..%d: %s' % (len(imgs), text[:60]))
        return imgs
    cycles = re.findall(r'\(([^()]*)\)', t)
    if t and not cycles and t != '()':
        raise ValueError('cannot read the permutation: %s' % text[:60])
    pts = [int(x) for c in cycles for x in re.findall(r'\d+', c)]
    n = max([d or 0] + pts)
    p = list(range(n))
    for c in cycles:
        xs = [int(x) - 1 for x in re.findall(r'\d+', c)]
        for i, x in enumerate(xs):
            p[x] = xs[(i + 1) % len(xs)]
    if sorted(p) != list(range(n)):
        raise ValueError('the cycles overlap: %s' % text[:60])
    return p


def pad(p, d):
    return p + list(range(len(p), d))


def then(p, q):
    """first p, then q (right actions: point^(pq) = (point^p)^q)"""
    return [q[p[i]] for i in range(len(p))]


def inverse(p):
    r = [0] * len(p)
    for i, x in enumerate(p):
        r[x] = i
    return r


def cycle_lengths(p):
    seen, out = set(), []
    for i in range(len(p)):
        if i not in seen:
            j, n = i, 0
            while j not in seen:
                seen.add(j)
                j = p[j]
                n += 1
            out.append(n)
    return out


def orbits(gens, d):
    seen, out = set(), []
    for s in range(d):
        if s in seen:
            continue
        orb, q = [s], deque([s])
        seen.add(s)
        while q:
            p = q.popleft()
            for g in gens:
                r = g[p]
                if r not in seen:
                    seen.add(r)
                    orb.append(r)
                    q.append(r)
        out.append(orb)
    return out


# ---------------------------------------------------------------- canonical forms
def canonical_key(perm, d):
    """canonical form of the transitive <a,b,c>-set (the same function as batch.canonical_key)"""
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


def canonical_xy(x, y):
    """canonical form of a transitive <x,y>-set (classes of subgroups of Delta up to conjugation in Delta)"""
    d = len(x)
    best = None
    for s in range(d):
        lab = {s: 0}
        order = [s]
        q = deque([s])
        while q:
            p = q.popleft()
            for g in (x, y):
                r = g[p]
                if r not in lab:
                    lab[r] = len(order)
                    order.append(r)
                    q.append(r)
        key = (tuple(lab[x[p]] for p in order), tuple(lab[y[p]] for p in order))
        if best is None or key < best:
            best = key
    return best


def class_id(key):
    """the fingerprint of the table groups/surfaces_summary.csv (same as summary_tables.class_id)"""
    return hashlib.sha1(json.dumps([list(k) for k in key]).encode()).hexdigest()[:10]


def abc_from_xy(x, y):
    """the action of a, b, c on the 2n cosets of K in Delta^pm, from the action of x = ab, y = bc on the n cosets
    of K in Delta. Point i is the coset K h (h in Delta), point n + i the coset K h a."""
    n = len(x)
    z = inverse(then(x, y))                      # z = ca, since (ab)(bc)(ca) = 1
    zi = inverse(z)
    a = [i + n for i in range(n)] + list(range(n))
    b = [None] * (2 * n)
    c = [None] * (2 * n)
    for i in range(n):
        # b = a.x: (Kh)b = (Kh a)x = (K h (a x a)) a = (K h x^-1) a ; (Kha)b = (Kh) x
        b[i] = n + x[i]                          # x is an involution, x^-1 = x
        b[n + i] = x[i]
        # c = z.a: (Kh)c = (Kh z) a ; (Kha)c = (K h (a z a)) = (K h z^-1)
        c[i] = n + z[i]
        c[n + i] = zi[i]
    return {'a': a, 'b': b, 'c': c}


def xy_of_class(perm, d):
    """the action of x = ab, y = bc on the orbit of the point 0 (the cosets of K in Delta, for a Fuchsian K)"""
    x = then(perm['a'], perm['b'])
    y = then(perm['b'], perm['c'])
    orb = orbits([x, y], d)
    o = [o for o in orb if 0 in o][0]
    pos = {p: i for i, p in enumerate(o)}
    return [pos[x[p]] for p in o], [pos[y[p]] for p in o]


# ---------------------------------------------------------------- checks
def check_abc(N, a, b, c):
    """returns (problems, info)"""
    d = len(a)
    problems = []
    for name, p, m in (('a', a, 2), ('b', b, 2), ('c', c, 2), ('ab', then(a, b), 2), ('bc', then(b, c), 3),
                       ('ca', then(c, a), N)):
        ls = cycle_lengths(p)
        if any(m % l for l in ls):
            problems.append('%s does not have order dividing %d (cycle lengths %s)' % (name, m, sorted(set(ls))))
        elif any(l != m for l in ls):
            problems.append('%s has cycles of length %s < %d: K contains %s' % (
                name, sorted(set(l for l in ls if l != m)), m,
                'a reflection' if m == 2 and len(name) == 1 else 'an elliptic element'))
    if len(orbits([a, b, c], d)) != 1:
        problems.append('the permutations do not generate a transitive group (the covering is not connected)')
    two = len(orbits([then(a, b), then(b, c)], d)) == 2
    chi = d // (2 * N) + d // 6 + d // 4 - 3 * d // 2 + d
    k = d // (2 * N)
    info = {'d': d, 'chi': chi, 'k': k, 'orientable': two,
            'genus': (2 - chi) // 2 if two else 2 - chi}
    return problems, info


def check_xy(N, x, y):
    n = len(x)
    problems = []
    for name, p, m in (('x', x, 2), ('y', y, 3), ('xy', then(x, y), N)):
        ls = cycle_lengths(p)
        if any(m % l for l in ls):
            problems.append('%s does not have order dividing %d (cycle lengths %s)' % (name, m, sorted(set(ls))))
        elif any(l != m for l in ls):
            problems.append('%s has cycles of length %s < %d: K contains an elliptic element' % (
                name, sorted(set(l for l in ls if l != m)), m))
    if len(orbits([x, y], n)) != 1:
        problems.append('the permutations do not generate a transitive group (the covering is not connected)')
    return problems


# ---------------------------------------------------------------- the catalogue
def load_catalogue():
    cat = {}
    for fam in FAMILIES:
        path = os.path.join(HERE, 'groups', fam + '.jsonl')
        if not os.path.exists(path):
            continue
        groups = [json.loads(l) for l in open(path) if l.strip()]
        res = json.load(open(os.path.join(HERE, 'groups', fam + '_results.json')))
        for n, (g, r) in enumerate(zip(groups, res), 1):
            if fam == 'nec_N18' and r['orientable']:
                continue          # the 8 orientable classes of that test file are the classes of fuchs_N18
            perm = {ch: [x - 1 for x in g[ch]] for ch in 'abc'}
            key = tuple(tuple(k) for k in r['key'])
            cat.setdefault((g['N'], key), []).append({'label': '%s#%d' % (fam, n), 'perm': perm, 'd': g['d'],
                                                      'family': fam})
    return cat


def table_lines(label):
    out = []
    if label in NAMES:
        out.append('  in the paper: ' + NAMES[label])
    p = os.path.join(HERE, 'groups', 'surfaces_summary.csv')
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            if r['surface'] == label:
                centres = [r[c] for c in r if c.startswith('centre_packing') and r[c]]
                anti = r['orientable'] == 'True' and r['induced_anticonformal'] == 'True'
                out.append('  packings: %s;  isometry group induced by the triangle group: order %s%s' % (
                    r['n_packings'], r['induced_isometries'], ', with an anticonformal isometry' if anti else ''))
                if centres:
                    out.append('  disc centres of the other packings (disc model): %s' % ' ; '.join(centres))
    p = os.path.join(HERE, 'groups', 'automorphism_groups.txt')
    if os.path.exists(p):
        for l in open(p):
            if l.startswith(label + ' '):
                out.append('  ' + ' '.join(l.split()[4:]).replace('Isom:', 'isometry group:').replace(
                    'Aut+:', ' conformal automorphisms:'))
    return out


def report(N, perm, xy=None, cat=None):
    """perm: action of a, b, c (0-based). xy: the (x, y) given by the user, if any. Returns the printed lines."""
    cat = cat if cat is not None else load_catalogue()
    d = len(perm['a'])
    key = canonical_key(perm, d)
    out = ['canonical form computed; class_id = %s' % class_id(key)]
    found = cat.get((N, key), [])
    if not found:
        out.append('This surface is not in the catalogue (the catalogue contains the k-extremal surfaces of genus 2 '
                   '(orientable) and 3 (non-orientable), and the external test nec_N18).')
        return out
    if len(found) == 1:
        e = found[0]
        out.append('It is %s.' % e['label'])
        out += table_lines(e['label'])
        return out
    # a mirror pair of Fuchsian classes
    labels = [e['label'] for e in found]
    out.append('It is one of the classes %s, which uniformize mirror-image surfaces (same canonical form in '
               'Delta^pm).' % ' and '.join(labels))
    if xy is not None:
        mine = canonical_xy(*xy)
        mirror = canonical_xy(xy[0], inverse(xy[1]))
        same = [e['label'] for e in found if canonical_xy(*xy_of_class(e['perm'], e['d'])) == mine]
        other = [e['label'] for e in found if canonical_xy(*xy_of_class(e['perm'], e['d'])) == mirror]
        out.append('  With the conventions of the paper (right actions, x = ab, y = bc) it is %s;' % ', '.join(same))
        out.append('  with the opposite convention (left actions, or y^-1 instead of y) it is %s.' % ', '.join(other))
    else:
        out.append('  (They are mirror images of each other; to tell them apart give the monodromy as x, y.)')
    for e in found:
        out.append('%s:' % e['label'])
        out += table_lines(e['label'])
    return out


# ---------------------------------------------------------------- command line
def main(argv):
    if len(argv) >= 1 and argv[0] == '--selftest':
        return selftest()
    if len(argv) < 2 or argv[1] not in ('abc', 'xy'):
        print(__doc__)
        return 1
    N = int(argv[0])
    if argv[1] == 'abc':
        if len(argv) != 5:
            sys.exit('give N, the word abc and three permutations')
        ps = [parse_perm(t) for t in argv[2:5]]
        d = max(len(p) for p in ps)
        a, b, c = (pad(p, d) for p in ps)
        problems, info = check_abc(N, a, b, c)
        xy = None
    else:
        if len(argv) != 4:
            sys.exit('give N, the word xy and two permutations')
        ps = [parse_perm(t) for t in argv[2:4]]
        n = max(len(p) for p in ps)
        x, y = (pad(p, n) for p in ps)
        if all(l == 3 for l in cycle_lengths(x)) and all(l == 2 for l in cycle_lengths(y)):
            x, y = y, x
            print('(the first permutation has order 3 and the second order 2: they are used as y and x)')
        problems = check_xy(N, x, y)
        abc = abc_from_xy(x, y)
        a, b, c = abc['a'], abc['b'], abc['c']
        p2, info = check_abc(N, a, b, c)
        problems += [q for q in p2 if q not in problems and not problems]
        xy = (x, y)
    print('N = %d: %d sheets over the orbifold of Delta^pm(2,3,%d)' % (N, info['d'], N))
    if problems:
        print('The permutations do not define a compact surface uniformized by a subgroup of Delta^pm(2,3,%d):' % N)
        for q in problems:
            print('  * ' + q)
        return 2
    print('%s surface of genus %d, Euler characteristic %d, index %d in Delta%s(2,3,%d), made of %d regular %d-gon%s:'
          ' a %d-extremal surface' % (
              'orientable' if info['orientable'] else 'non-orientable', info['genus'], info['chi'],
              info['d'] // 2 if info['orientable'] else info['d'], '' if info['orientable'] else '^pm', N,
              info['k'], N, 's' if info['k'] > 1 else '', info['k']))
    for l in report(N, {'a': a, 'b': b, 'c': c}, xy):
        print(l)
    return 0


def selftest():
    """every stored class: (1) the stored canonical form is recomputed; (2) the class is identified from its own
    monodromy after a random renumbering of the points, given as a, b, c and, for Fuchsian classes, as x, y"""
    random.seed(1)
    cat = load_catalogue()
    n = bad = 0
    for (N, key), entries in cat.items():
        for e in entries:
            n += 1
            perm, d = e['perm'], e['d']
            if canonical_key(perm, d) != key:
                bad += 1
                print('stored key differs:', e['label'])
            s = list(range(d))
            random.shuffle(s)
            si = inverse(s)
            rp = {ch: [s[perm[ch][si[i]]] for i in range(d)] for ch in 'abc'}
            if canonical_key(rp, d) != key:
                bad += 1
                print('renumbering changes the key:', e['label'])
            if e['family'].startswith('fuchs'):
                x, y = xy_of_class(perm, d)
                m = len(x)
                t = list(range(m))
                random.shuffle(t)
                ti = inverse(t)
                rx = [t[x[ti[i]]] for i in range(m)]
                ry = [t[y[ti[i]]] for i in range(m)]
                lines = report(N, abc_from_xy(rx, ry), (rx, ry), cat)
                if len(entries) == 1:
                    ok = ('It is %s.' % e['label']) in lines
                else:
                    ok = any(l.startswith('  With the conventions of the paper') and e['label'] in l for l in lines)
                if not ok:
                    bad += 1
                    print('not identified from x, y:', e['label'])
    print('%d classes checked, %d problems' % (n, bad))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
