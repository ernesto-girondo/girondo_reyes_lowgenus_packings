"""Weierstrass points of X0 = fuchs_N18#9 = fuchs_N9#36 = fuchs_N9#37 (genus 2, D6 curve).
In the frame of K37 the hyperelliptic involution is induced by ab; its fixed points are 6 edge midpoints.
We transport them to the frames of K36 and K9 (N=18) through the hidden packings."""
import json, mpmath as mp
from field import pari
from geom import Geometry, to_disc, hdist, act, mmul, minv
from surface import Surface
from gapiface import coset_action_from_words, reduce_word
from batch import parity_and_words, automorphisms
from wordfmt import compress
from outfile import tee
tee('groups/weierstrass_X0.txt')
mp.mp.dps = 40
i0 = mp.mpc(0, 1)

def load(fn, idx, G):
    g = [json.loads(l) for l in open(fn) if l.strip()][idx - 1]
    S = Surface(G, {'d': g['d'], 'a': [x-1 for x in g['a']], 'b': [x-1 for x in g['b']], 'c': [x-1 for x in g['c']]}, label='%s#%d' % (fn, idx))
    S.analyse(verbose=False)
    return S

def triangle_word(S, p):
    """word g_p with triangle g_p T0 representing coset point p inside the fundamental domain"""
    for Q, poly in enumerate(S.polys):
        if p in poly['members']:
            m, e = poly['members'][p]
            return S.gword[Q] + 'ca' * m + 'a' * e
    raise KeyError

def reduce_mod_K(S, z):
    """K-equivalent point inside the fundamental domain of S (frame of S.G)"""
    G = S.G
    # Delta^pm reduction of z: find word w with z = w(z0)
    q = z
    word = []
    while True:
        v = G.wall_violations(q)
        if not v: break
        q = act(G.gens[v[0]], q); word.append(v[0])
    w = ''.join(word)            # z = w(q)  (w = r1 r2 ... rn applied right to left)
    p = S.apply(0, w)
    gp = triangle_word(S, p)
    return act(G.word_num(gp), q), gp

G9 = Geometry(9); G18 = Geometry(18)
S37 = load('groups/fuchs_N9.jsonl', 37, G9)
S36 = load('groups/fuchs_N9.jsonl', 36, G9)
# J in the frame of K37: automorphism with word 'ab'
par, word = parity_and_words(S37.perm, S37.d)
j = S37.apply(0, 'ab')
phi = [S37.apply(j, word[p]) for p in range(S37.d)]
fixed_tri = [p for p in range(S37.d) if S37.apply(p, 'ab') == phi[p]]
mids = []
for p in fixed_tri:
    gp = triangle_word(S37, p)
    m = act(G9.word_num(gp), mp.mpc(0, 1 / G9.u))
    if not any(abs(m - q) < 1e-12 for q in mids):
        mids.append(m)
print('Frame K37 (4-packing invariant under Aut): Weierstrass points = edge midpoints of the 9-gons:')
for m in mids:
    print('   disc %s   d(0,w)=%s  angle=%s deg' % (mp.nstr(to_disc(m), 8), mp.nstr(hdist(i0, m), 8), mp.nstr(mp.degrees(mp.arg(to_disc(m))), 6)))

def transport(S_from, S_to, N_to, label):
    """find in S_from a hidden packing whose K_j is the class of S_to; return the isometry T = alpha*delta
    with T K_to T^-1 = K_from (so T maps the frame of S_to onto the frame of S_from)"""
    from batch import canonical_key
    key_to = canonical_key(S_to.perm, S_to.d)
    for P in S_from.packings:
        if S_from.G.reduce(P['alpha']) is not None: continue
        a2 = coset_action_from_words(N_to, [reduce_word(w) for w in P['words']])
        if canonical_key({'a': a2['a'], 'b': a2['b'], 'c': a2['c']}, a2['d']) != key_to: continue
        perm2 = {'a': a2['a'], 'b': a2['b'], 'c': a2['c']}
        def ap(p, w):
            for ch in w: p = perm2[ch][p]
            return p
        qs = [q for q in range(a2['d']) if all(ap(q, s.word) == q for s in S_to.sp)]
        par2, word2 = parity_and_words(perm2, a2['d'])
        delta = S_to.G.word_num(word2[qs[0]])
        T = mmul(P['alpha'], delta)
        return T
    return None

for S_to_name, S_from, S_to, G_to in [('K36 (frame of a 4-packing of the orbit of 4)', S36, S37, G9)]:
    pass

# transport Weierstrass points into the frames of K36 and K9
S9 = load('groups/fuchs_N18.jsonl', 9, G18)
# S9's packings in the N=9 geometry: redo analysis with G9 (as in cross_check)
exec(open('cross_check_9_18.py').read().split("from gapiface import")[0].replace("S.analyse(verbose=True)", "S.analyse(verbose=False)").replace("G18 = Geometry(18); G9 = Geometry(9)", ""))
S9x = S   # fundamental domain of K9 (N=18) but packings/alphas in the (2,3,9) geometry
for name, S_from, G_from in [('K36', S36, G9), ('K9 (N=18, extremal disc)', S9x, G18)]:
    T = transport(S_from, S37, 9, name)
    print('Frame %s: Weierstrass points (reduced into the fundamental domain):' % name)
    S_from.G = G_from
    for m in mids:
        w = act(T, m)
        wF, gp = reduce_mod_K(S_from, w)
        print('   disc %s   d(0,w)=%s  angle=%s deg   triangle word %s' % (mp.nstr(to_disc(wF), 8), mp.nstr(hdist(i0, wF), 8), mp.nstr(mp.degrees(mp.arg(to_disc(wF))), 6), compress(gp)))
    S_from.G = G9
