"""All automorphisms of D/K (hidden ones included) as explicit isometries, their orders, fixed points."""
import json, mpmath as mp
from collections import deque
from geom import Geometry, to_disc, hdist, act, mmul, minv, mnorm, mdet, proj_equal
from surface import Surface
from gapiface import coset_action_from_words, reduce_word
from batch import parity_and_words, canonical_key
mp.mp.dps = 40
i0 = mp.mpc(0, 1)

def triangle_word(S, p):
    for Q, poly in enumerate(S.polys):
        if p in poly['members']:
            m, e = poly['members'][p]
            return S.gword[Q] + 'ca' * m + 'a' * e

def reduce_mod_K(S, z):
    G = S.G
    q = z; word = []
    while True:
        v = G.wall_violations(q)
        if not v: break
        q = act(G.gens[v[0]], q); word.append(v[0])
    p = S.apply(0, ''.join(word))
    return act(G.word_num(triangle_word(S, p)), q)

def in_K(S, M):
    """is the isometry M (numeric) in K?  reduce in Delta^pm then check the word fixes point 0"""
    w = S.G.reduce(M)
    return w is not None and S.apply(0, w) == 0

def all_automorphisms(S):
    """list of (matrix, source packing index, parity)"""
    G = S.G
    key0 = canonical_key(S.perm, S.d)
    par, word = parity_and_words(S.perm, S.d)
    auts = []
    # induced
    for j in range(S.d):
        if all(S.apply(j, s.word) == j for s in S.sp):
            auts.append((G.word_num(word[j]), 'induced', par[j]))
    for pi, P in enumerate(S.packings):
        if G.reduce(P['alpha']) is not None: continue
        a2 = coset_action_from_words(G.N, [reduce_word(w) for w in P['words']])
        if a2['d'] != S.d or canonical_key({'a': a2['a'], 'b': a2['b'], 'c': a2['c']}, a2['d']) != key0: continue
        perm2 = {'a': a2['a'], 'b': a2['b'], 'c': a2['c']}
        def ap(p, w):
            for ch in w: p = perm2[ch][p]
            return p
        par2, word2 = parity_and_words(perm2, a2['d'])
        for q in range(a2['d']):
            if all(ap(q, s.word) == q for s in S.sp):
                g = mmul(P['alpha'], G.word_num(word2[q]))
                auts.append((g, 'packing %d' % pi, par2[q]))
    return auts

def order_of(S, g, maxo=50):
    M = mnorm(g)
    for n in range(1, maxo + 1):
        if in_K(S, M): return n
        M = mnorm(mmul(M, g))
    return None

def local_elements(S, g, nsamp=60):
    """elements h = kappa*g (kappa in K) that map sample points of the fundamental domain back into it"""
    G = S.G
    import random
    random.seed(0)
    A, Bv, C = G.vA, G.vB, G.vC
    samples = []
    while len(samples) < nsamp:
        r1, r2 = random.random(), random.random()
        if r1 + r2 > 1: r1, r2 = 1 - r1, 1 - r2
        q = A + r1 * (Bv - A) + r2 * (C - A)
        if not G.wall_violations(q): samples.append(q)
    hs = []
    for p in range(S.d):
        gp = G.word_num(triangle_word(S, p))
        for q in samples:
            z = act(gp, q)
            w = act(g, z)
            # reduce w
            qq = w; word = []
            while True:
                v = G.wall_violations(qq)
                if not v: break
                qq = act(G.gens[v[0]], qq); word.append(v[0])
            vword = ''.join(word)
            pp = S.apply(0, vword)
            kappa = mmul(G.word_num(triangle_word(S, pp)), minv(G.word_num(vword)))
            h = mnorm(mmul(kappa, g))
            if not any(proj_equal(h, h2) for h2 in hs):
                hs.append(h)
    return hs

def fixed_points(S, g, hs):
    pts = []
    for h in hs:
        if mdet(h) < 0: continue
        tr = h[0][0] + h[1][1]
        if abs(tr) >= 2 - 1e-12: continue
        a, b, c, d = h[0][0], h[0][1], h[1][0], h[1][1]
        disc = (d - a) ** 2 + 4 * b * c
        z = (-(d - a) + mp.sqrt(disc)) / (2 * c)
        if z.imag < 0: z = (-(d - a) - mp.sqrt(disc)) / (2 * c)
        zF = reduce_mod_K(S, z)
        if not any(abs(zF - w) < 1e-9 for w in pts):
            pts.append(zF)
    return pts

if __name__ == '__main__':
    import sys
    fn, idx = sys.argv[1], int(sys.argv[2])
    g = [json.loads(l) for l in open(fn) if l.strip()][idx - 1]
    G = Geometry(g['N'])
    S = Surface(G, {'d': g['d'], 'a': [x-1 for x in g['a']], 'b': [x-1 for x in g['b']], 'c': [x-1 for x in g['c']]})
    S.analyse(verbose=False)
    auts = all_automorphisms(S)
    print('%s#%d: |Aut| = %d' % (fn, idx, len(auts)))
    for gm, src, parity in auts:
        o = order_of(S, gm)
        if parity == 0 and o == 2:
            fp = fixed_points(S, gm, local_elements(S, gm))
            print('  conformal involution from %-10s: %d fixed points: %s' % (src, len(fp), [(mp.nstr(to_disc(z), 6), mp.nstr(hdist(i0,z),6), mp.nstr(mp.degrees(mp.arg(to_disc(z))),5)) for z in fp]))
        else:
            print('  %s from %-10s order %s' % ('conformal' if parity == 0 else 'anticonformal', src, o))
