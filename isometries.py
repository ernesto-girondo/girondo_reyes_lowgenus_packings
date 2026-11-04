"""The isometry group of a surface D/K and the fixed points of its elements, for the pictures of draw_surface.py.

For every isometry gamma of X = D/K (induced by the triangle group or hidden, from autfix.all_automorphisms) we
determine
  * its isolated fixed points on X (points where gamma acts locally as a rotation), and
  * its mirrors (curves of fixed points, where gamma acts locally as a reflection).
A point of X fixed by gamma is the projection of the fixed set of an element h of the coset K gamma.

Induced isometries (gamma in the triangle group) are treated exactly: the elements of K gamma with fixed points are
conjugates w r w^-1 of the rotations r about the vertices O, A, C of the triangle T_0 and of the reflections a, b, c
in its sides, w running over the triangles of the fundamental domain, and the condition w r w^-1 in K gamma is a
condition on the permutation representation.
Hidden isometries are treated numerically: the elements h of K gamma moving some point of the fundamental domain
F back into F are found by sampling F, and the list is closed under multiplication by the side pairings.

As a check, the formula chi(X/G) = (chi(X) + sum_{gamma != 1} #isolated fixed points of gamma) / |G| (mirrors are
circles and do not contribute) must give an integer, both for the whole group and, for orientable X, for the group
of conformal automorphisms, where it must moreover be even and at most 2 (the quotient is a closed surface).
"""
import random
import cmath
import mpmath as mp
from geom import act, mmul, minv, mnorm, mdet
from autfix import all_automorphisms, in_K, triangle_word, reduce_mod_K

I0 = mp.mpc(0, 1)


# ---------------------------------------------------------------- floating point helpers (pictures only)
def fmat(M):
    return [[complex(M[0][0]), complex(M[0][1])], [complex(M[1][0]), complex(M[1][1])]]


def fdet(M):
    return (M[0][0] * M[1][1] - M[0][1] * M[1][0]).real


def fact(M, z):
    """action on the upper half plane (orientation reversing matrices act on conj(z))"""
    if fdet(M) < 0:
        z = z.conjugate()
    return (M[0][0] * z + M[0][1]) / (M[1][0] * z + M[1][1])


def fmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def finv(A):
    d = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return [[A[1][1] / d, -A[0][1] / d], [-A[1][0] / d, A[0][0] / d]]


def fdist(z, w):
    return 2 * cmath.atanh(abs((z - w) / (z - w.conjugate()))).real


def todisc(z):
    return 1j * (z - 1j) / (z + 1j)


def fromdisc(w):
    """inverse of todisc"""
    return (1 - 1j * w) / (w - 1j)


class Domain:
    """the fundamental domain F (union of k polygons) in floating point"""

    def __init__(self, S):
        G = S.G
        self.S = S
        self.Ms = [fmat(G.word_num(S.gword[P])) for P in sorted(S.gword)]
        self.Minv = [finv(M) for M in self.Ms]
        self.nb = [complex(act(R, I0)) for R in G.Rota]           # centres of the neighbours of P_0
        self.sp = []
        for s in S.sp:
            M = fmat(s.Mn)
            self.sp += [M, finv(M)]

    def contains(self, z, tol=1e-9):
        for Mi in self.Minv:
            q = fact(Mi, z)
            d0 = fdist(q, 1j)
            if all(d0 <= fdist(q, c) + tol for c in self.nb):
                return True
        return False

    def copies(self, z, tol=1e-7):
        """all the copies in F of the point z of F (closure under the side pairings)"""
        out = [z]
        i = 0
        while i < len(out):
            for M in self.sp:
                w = fact(M, out[i])
                if self.contains(w) and not any(abs(w - u) < tol * (1 + abs(u)) for u in out):
                    out.append(w)
            i += 1
            if len(out) > 60:
                break
        return out

    def canonical(self, z):
        cs = self.copies(z)
        ds = [todisc(c) for c in cs]
        return min(ds, key=lambda w: (round(w.real, 7), round(w.imag, 7))), [todisc(c) for c in cs]


# ---------------------------------------------------------------- elements of K gamma with fixed points
def classify(h):
    """('rot', fixed point) for an elliptic element, ('mirror', two points of the axis) for a reflection, else None"""
    det = fdet(h)
    a, b, c, d = h[0][0], h[0][1], h[1][0], h[1][1]
    s = cmath.sqrt(abs(det))
    a, b, c, d = a / s, b / s, c / s, d / s
    tr = (a + d).real
    if det > 0:
        if abs(tr) >= 2 - 1e-9:
            return None
        disc = (d - a) ** 2 + 4 * b * c
        if abs(c) < 1e-14:
            return None
        z = (-(d - a) + cmath.sqrt(disc)) / (2 * c)
        if z.imag < 0:
            z = (-(d - a) - cmath.sqrt(disc)) / (2 * c)
        return ('rot', z)
    else:
        if abs(tr) > 1e-7:
            return None                                  # glide reflection: no fixed points
        pts = []
        for w in (1j, 2 + 1j, -1 + 0.5j, 0.3 + 3j):
            hw = fact(h, w)
            # midpoint of the segment [w, h(w)] lies on the mirror
            pts.append(geo_midpoint(w, hw))
        p, q = pts[0], max(pts[1:], key=lambda u: abs(todisc(u) - todisc(pts[0])))
        return ('mirror', (p, q))


def geo_midpoint(z, w):
    # move z to i, take the midpoint on the imaginary axis direction, move back
    zd, wd = todisc(z), todisc(w)
    u = (wd - zd) / (1 - zd.conjugate() * wd)
    r = abs(u)
    if r < 1e-15:
        return z
    rm = cmath.tanh(cmath.atanh(r) / 2).real
    m = rm * u / r
    return fromdisc((m + zd) / (1 + zd.conjugate() * m))


def induced_elements(S, word_g):
    """exact: elements w r w^-1 in K gamma (gamma = word_g), with r a rotation about O, A, C or a reflection a, b, c.
    Returns a list of (kind, data) as in classify, with exact special points."""
    G = S.G
    target = S.apply(0, word_g)
    out = []
    rots = [('O', 'ca' * m) for m in range(1, G.N)] + [('A', 'bc'), ('A', 'bcbc'), ('C', 'ab')]
    pt = {'O': I0, 'A': G.vA, 'C': G.vC}
    for p in range(S.d):
        w = triangle_word(S, p)
        winv = w[::-1]
        W = G.word_num(w)
        for P, r in rots:
            if S.apply(0, w + r + winv) == target:
                out.append(('rot', complex(act(W, pt[P]))))
        for r, (P, Q) in (('a', ('O', 'C')), ('b', ('C', 'A')), ('c', ('A', 'O'))):
            if S.apply(0, w + r + winv) == target:
                out.append(('mirror', (complex(act(W, pt[P])), complex(act(W, pt[Q])))))
    return out


def hidden_elements(S, g, dom, nsamp=8, seed=0):
    """numerical: elements kappa*g (kappa in K) moving sample points of F into F, closed under side pairings"""
    G = S.G
    rnd = random.Random(seed)
    A, B, C = G.vA, G.vB, G.vC
    samples = []
    while len(samples) < nsamp:
        r1, r2 = rnd.random(), rnd.random()
        if r1 + r2 > 1:
            r1, r2 = 1 - r1, 1 - r2
        q = A + r1 * (B - A) + r2 * (C - A)
        if not G.wall_violations(q):
            samples.append(q)
    hs = []                                   # pairs (float matrix, high precision matrix)

    def add(hm):
        hn = fmat(hm)
        for h2, _ in hs:
            if all(abs(hn[i][j] * h2[1][1] - h2[i][j] * hn[1][1]) < 1e-7 * (1 + abs(hn[1][1])) for i in range(2)
                   for j in range(2)) and abs(hn[1][0] * h2[0][0] - h2[1][0] * hn[0][0]) < 1e-7:
                return False
        hs.append((hn, hm))
        return True

    for p in range(S.d):
        gp = G.word_num(triangle_word(S, p))
        for q in samples:
            z = act(gp, q)
            w = act(g, z)
            qq, word = w, []
            while True:
                v = G.wall_violations(qq)
                if not v:
                    break
                qq = act(G.gens[v[0]], qq)
                word.append(v[0])
            vword = ''.join(word)
            pp = S.apply(0, vword)
            kappa = mmul(G.word_num(triangle_word(S, pp)), minv(G.word_num(vword)))
            add(mnorm(mmul(kappa, g)))
    sides = [s.Mn for s in S.sp] + [minv(s.Mn) for s in S.sp]
    base = list(hs)
    for _, hm in base:
        for M in sides:
            add(mnorm(mmul(M, hm)))
    out = []
    for hf, hm in hs:
        c = classify(hf)
        if c is None:
            continue
        if c[0] == 'rot':
            # fixed point in high precision, then brought into F (reduction modulo K)
            a_, b_, c_, d_ = hm[0][0], hm[0][1], hm[1][0], hm[1][1]
            disc = (d_ - a_) ** 2 + 4 * b_ * c_
            z = (-(d_ - a_) + mp.sqrt(disc)) / (2 * c_)
            if mp.im(z) < 0:
                z = (-(d_ - a_) - mp.sqrt(disc)) / (2 * c_)
            out.append(('rot', complex(reduce_mod_K(S, z))))
        else:
            out.append(c)
    return out


# ---------------------------------------------------------------- order of an isometry of X
def reduce_mod_K_matrix(S, M):
    """kappa*M with kappa in K such that kappa*M(i) lies in F (keeps the entries small)"""
    G = S.G
    q, word = act(M, I0), []
    while True:
        v = G.wall_violations(q)
        if not v:
            break
        q = act(G.gens[v[0]], q)
        word.append(v[0])
    vword = ''.join(word)
    kappa = mmul(G.word_num(triangle_word(S, S.apply(0, vword))), minv(G.word_num(vword)))
    return mnorm(mmul(kappa, M))


def order_on_X(S, g, maxo=200):
    M = reduce_mod_K_matrix(S, g)
    for n in range(1, maxo + 1):
        if in_K(S, M):
            return n
        M = reduce_mod_K_matrix(S, mmul(M, g))
    raise RuntimeError('order not found')


# ---------------------------------------------------------------- the whole analysis
def analyse_isometries(S):
    """returns a dict with the group data and the fixed point data, ready to be drawn"""
    G = S.G
    dom = Domain(S)
    auts = all_automorphisms(S)                       # list of (matrix, source, parity)
    elems = []
    for gm, src, parity in auts:
        o = order_on_X(S, gm)
        if o == 1:
            continue                                  # the identity
        if src == 'induced':
            items = induced_elements(S, G.reduce(gm))
        else:
            items = hidden_elements(S, gm, dom)
        isolated, mirrors = {}, []
        for kind, data in items:
            if kind == 'rot':
                key, allc = dom.canonical(data)
                isolated[(round(key.real, 6), round(key.imag, 6))] = (key, allc)
            else:
                mirrors.append(data)
        elems.append({'matrix': gm, 'source': src, 'parity': parity, 'order': o,
                      'isolated': isolated, 'mirrors': mirrors})
    n = len(elems) + 1
    # chi(X) from the triangulation of X by the d triangles (vertices: the points O, A, C)
    d, N = S.d, G.N
    chiX = d // (2 * N) + d // 6 + d // 4 - 3 * d // 2 + d
    tot = chiX + sum(len(e['isolated']) for e in elems)
    res = {'order': n, 'elements': elems, 'chiX': chiX, 'chi_quotient': tot / n}
    if S.orientable:
        conf = [e for e in elems if e['parity'] == 0]
        res['n_conformal'] = len(conf) + 1
        res['chi_quotient_conformal'] = (chiX + sum(len(e['isolated']) for e in conf)) / (len(conf) + 1)
    # points with non trivial rotation stabiliser
    pts = {}
    for e in elems:
        for key, (z, allc) in e['isolated'].items():
            pts.setdefault(key, [z, allc, 1])[2] += 1
    res['points'] = list(pts.values())          # [canonical point (disc), all copies (disc), stabiliser order]
    return res, dom
