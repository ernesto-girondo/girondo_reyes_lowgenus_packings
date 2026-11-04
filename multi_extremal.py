"""Theorem 6.3(v): which 1-extremal surfaces of genus 2 (classes of fuchs_N18) are also 4-extremal?

A surface D/K with K < Delta(2,3,18) of index 18 is 4-extremal iff K is contained in some PSL2(R)-conjugate
alpha Delta(2,3,9) alpha^-1. The two triangle groups are commensurable and L_18 = L_9 (same defining polynomial,
different embeddings), so the side pairings of K can be rewritten over L_9 and the 4-packing search of surface.py
can be run with the geometry of Delta(2,3,9). For every 4-packing found, K_j = alpha^-1 K alpha < Delta(2,3,9) is
identified with its class in groups/fuchs_N9.jsonl (Delta^pm-key and Delta-key).
Generalizes cross_check_9_18.py (which only treats the 9th class, X0).
Usage: python3 multi_extremal.py          (all 9 classes of fuchs_N18, about one minute)
"""
import json
import mpmath as mp
from field import pari
from geom import Geometry, to_disc
from surface import Surface
from gapiface import coset_action_from_words, reduce_word
from batch import canonical_key
from identify2 import delta_key
from outfile import tee
tee('groups/multi_extremal.txt')

G18 = Geometry(18); G9 = Geometry(9)
F18, F9 = G18.F, G9.F
assert F18.P == F9.P
nf = pari.nfinit(F9.P)
r = None
for rr in pari.nfroots(nf, pari.subst(F9.P, 'y', 'x')):
    cand = pari('Mod(%s, %s)' % (pari.lift(rr), F9.P))
    if abs(F9.num(cand) - F18.t_num) < mp.mpf(10) ** -25:
        r = cand
assert r is not None, 'field isomorphism L_18 -> L_9 not found'


def conv(e):
    return pari('Mod(%s, %s)' % (pari.subst(pari.lift(e), 'y', pari.lift(r)), F9.P))


def conv_mat(M):
    return pari.matrix(2, 2, [conv(M[0, 0]), conv(M[0, 1]), conv(M[1, 0]), conv(M[1, 1])])


def action_of(g):
    return {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}


res9 = json.load(open('groups/fuchs_N9_results.json'))
keys9 = {tuple(map(tuple, x['key'])): x['label'] for x in res9}
g9 = [json.loads(l) for l in open('groups/fuchs_N9.jsonl') if l.strip()]
dk9 = {delta_key(action_of(gg)): 'fuchs_N9#%d' % (i + 1) for i, gg in enumerate(g9)}
res18 = json.load(open('groups/fuchs_N18_results.json'))

g18 = [json.loads(l) for l in open('groups/fuchs_N18.jsonl') if l.strip()]
print('class        discs  4-packings  classes of K_j in Delta(2,3,9) (Delta-key)')
for i, g in enumerate(g18):
    S = Surface(G18, action_of(g), label='fuchs_N18#%d' % (i + 1))
    for s in S.sp:
        s.Mx = conv_mat(s.Mx)
        assert all(abs(F9.num(s.Mx[a, b]) - s.Mn[a][b]) < 1e-20 for a in range(2) for b in range(2))
    S.G = G9
    S.analyse(verbose=False)
    found = []
    for P in S.packings:
        a2 = coset_action_from_words(9, [reduce_word(w) for w in P['words']])
        assert a2['d'] == 72, 'K_j should have index 72 in Delta^pm(2,3,9)'
        found.append(dk9.get(delta_key(a2)) or keys9.get(canonical_key({'a': a2['a'], 'b': a2['b'], 'c': a2['c']}, a2['d'])))
    print('%-12s %5d  %10d  %s' % (S.label, res18[i]['n_packings'], len(S.packings), sorted(set(found)) if found else '-'))
