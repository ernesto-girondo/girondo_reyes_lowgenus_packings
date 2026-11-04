"""Does the 1-extremal surface fuchs_N18#9 also admit extremal 4-packings (subgroup of a conjugate of Delta(2,3,9))?"""
import json, mpmath as mp
from field import pari, pmpf
from geom import Geometry, to_disc, hdist
from surface import Surface
if globals().get('__file__', '').endswith('cross_check_9_18.py'):     # not when this code is reused by other scripts
    from outfile import tee; tee('groups/cross_check_9_18.txt')
G18 = Geometry(18); G9 = Geometry(9)
F18, F9 = G18.F, G9.F
assert F18.P == F9.P
# field isomorphism matching the embeddings: root r of P in Q[y]/P with r(t9) = t18
nf = pari.nfinit(F9.P)
roots = pari.nfroots(nf, pari.subst(F9.P, 'y', 'x'))
r = None
for rr in roots:
    cand = pari('Mod(%s, %s)' % (pari.lift(rr), F9.P))
    if abs(F9.num(cand) - F18.t_num) < mp.mpf(10) ** -25:
        r = cand
print('isomorphism found:', r is not None)
def conv(e):
    return pari('Mod(%s, %s)' % (pari.subst(pari.lift(e), 'y', pari.lift(r)), F9.P))
def conv_mat(M):
    return pari.matrix(2, 2, [conv(M[0, 0]), conv(M[0, 1]), conv(M[1, 0]), conv(M[1, 1])])
g = [json.loads(l) for l in open('groups/fuchs_N18.jsonl') if l.strip()][8]
action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
S = Surface(G18, action, label='N18#9 as a subgroup of conjugates of Delta(2,3,9)')
for s in S.sp:
    s.Mx = conv_mat(s.Mx)
    assert all(abs(F9.num(s.Mx[i, j]) - s.Mn[i][j]) < 1e-20 for i in range(2) for j in range(2))
S.G = G9
S.analyse(verbose=True)

from gapiface import coset_action_from_words, reduce_word
from batch import canonical_key
from identify2 import delta_key
res9 = json.load(open('groups/fuchs_N9_results.json'))
keys9 = {tuple(map(tuple, r['key'])): r['label'] for r in res9}
g9 = [json.loads(l) for l in open('groups/fuchs_N9.jsonl') if l.strip()]
dk = {}
for i, gg in enumerate(g9):
    dk[delta_key({'d': gg['d'], 'a': [x-1 for x in gg['a']], 'b': [x-1 for x in gg['b']], 'c': [x-1 for x in gg['c']]})] = 'fuchs_N9#%d' % (i+1)
for P in S.packings:
    if G9.reduce(P['alpha']) is not None:
        print('alpha in Delta(2,3,9)?? unexpected'); continue
    a2 = coset_action_from_words(9, [reduce_word(w) for w in P['words']])
    k = canonical_key({'a': a2['a'], 'b': a2['b'], 'c': a2['c']}, a2['d'])
    print('4-packing centre', mp.nstr(to_disc(P['z']), 6), '-> K_j index', a2['d'], 'Delta+- class:', keys9.get(k), ' Delta class:', dk.get(delta_key(a2)))
