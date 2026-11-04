"""Interface to GAP: coset action of the extended triangle group on the cosets of a subgroup."""
import subprocess, json, os, tempfile


class W:
    """words in the involutions a,b,c (inverse = reversed string)"""
    def __init__(self, s=''):
        self.s = s

    def __mul__(self, o):
        return W(self.s + o.s)

    def inv(self):
        return W(self.s[::-1])

    def __pow__(self, n):
        if n < 0:
            return self.inv() ** (-n)
        return W(self.s * n)

    def __repr__(self):
        return self.s


def reduce_word(s):
    """free reduction (aa -> 1)"""
    out = []
    for ch in s:
        if out and out[-1] == ch:
            out.pop()
        else:
            out.append(ch)
    return ''.join(out)


GAP_TEMPLATE = r"""
N := %d;;
F := FreeGroup("a","b","c");;
a := F.1;; b := F.2;; c := F.3;;
G := F / [a^2, b^2, c^2, (a*b)^2, (b*c)^3, (c*a)^N];;
a := G.1;; b := G.2;; c := G.3;;
wordfromstring := function(s)
  local w, ch;
  w := One(G);
  for ch in s do
    if ch = 'a' then w := w*a; elif ch = 'b' then w := w*b; else w := w*c; fi;
  od;
  return w;
end;;
%s
hom := FactorCosetAction(G, K);;
d := Index(G, K);;
pa := ListPerm(Image(hom, a), d);; pb := ListPerm(Image(hom, b), d);; pc := ListPerm(Image(hom, c), d);;
PrintTo("%s", "{\"d\":", d, ",\"a\":", pa, ",\"b\":", pb, ",\"c\":", pc, "}");
QUIT;
"""


def coset_action_from_words(N, words):
    """words: list of letter strings generating K. Returns dict with index and perms (0-based lists)."""
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, 'out.json')
        gens = 'K := Subgroup(G, [%s]);;' % ', '.join('wordfromstring("%s")' % w for w in words)
        script = GAP_TEMPLATE % (N, gens, out)
        sp = os.path.join(td, 's.g')
        open(sp, 'w').write(script)
        subprocess.run(['gap', '-q', '-b', sp], check=True, capture_output=True, timeout=3600,
                       stdin=subprocess.DEVNULL)      # no input: on an error GAP stops instead of waiting
        res = json.load(open(out))
    return {'d': res['d'], 'a': [x - 1 for x in res['a']], 'b': [x - 1 for x in res['b']], 'c': [x - 1 for x in res['c']]}


def coset_action_from_perms(N, pa, pb, pc):
    """perms given as 0-based lists (already the coset action). Just validates relations."""
    d = len(pa)
    def comp(p, q):  # x -> q[p[x]]  (apply p then q) : right action convention
        return [q[p[x]] for x in range(d)]
    def power(p, n):
        r = list(range(d))
        for _ in range(n):
            r = comp(r, p)
        return r
    idp = list(range(d))
    assert power(pa, 2) == idp and power(pb, 2) == idp and power(pc, 2) == idp
    assert power(comp(pa, pb), 2) == idp
    assert power(comp(pb, pc), 3) == idp
    assert power(comp(pc, pa), N) == idp
    return {'d': d, 'a': pa, 'b': pb, 'c': pc}
