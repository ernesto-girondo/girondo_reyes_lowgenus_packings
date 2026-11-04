"""Analysis of one extremal surface D/K: fundamental domain, side pairings, search of all extremal
packings and exact verification  K < alpha Delta alpha^{-1}."""
import mpmath as mp
from geom import Geometry, mmul, mdet, minv, mnorm, act, hdist, to_disc, proj_equal
from field import pari
from gapiface import reduce_word

mp.mp.dps = 40
TOL = mp.mpf(10) ** -8
MARGIN = mp.mpf('0.1')


def circle3(p, q, r):
    """general circle A(x^2+y^2)+Bx+Cy+D=0 through three finite points (complex)"""
    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    pts = [p, q, r]
    s = [abs(t) ** 2 for t in pts]
    x = [t.real for t in pts]
    y = [t.imag for t in pts]
    A = det3([[x[0], y[0], 1], [x[1], y[1], 1], [x[2], y[2], 1]])
    B = -det3([[s[0], y[0], 1], [s[1], y[1], 1], [s[2], y[2], 1]])
    C = det3([[s[0], x[0], 1], [s[1], x[1], 1], [s[2], x[2], 1]])
    D = -det3([[s[0], x[0], y[0]], [s[1], x[1], y[1]], [s[2], x[2], y[2]]])
    return (A, B, C, D)


def circ_norm(c):
    A, B, C, D = c
    if abs(A) > mp.mpf(10) ** -20 * max(abs(B), abs(C), 1):
        return (mp.mpf(1), B / A, C / A, D / A)
    m = max(abs(B), abs(C))
    return (mp.mpf(0), B / m, C / m, D / m)


def line_circle(line, circ):
    _, B, C, D = line
    _, b, c, d = circ          # x^2+y^2+bx+cy+d = 0
    sols = []
    if abs(C) >= abs(B):
        # y = -(Bx+D)/C
        # x^2 + (Bx+D)^2/C^2 + b x - c (Bx+D)/C + d = 0
        k1, k0 = -B / C, -D / C
        qa = 1 + k1 ** 2
        qb = 2 * k1 * k0 + b + c * k1
        qc = k0 ** 2 + c * k0 + d
        disc = qb ** 2 - 4 * qa * qc
        if disc < -mp.mpf(10) ** -24:
            return []
        disc = mp.sqrt(max(disc, 0))
        for sgn in (1, -1):
            x = (-qb + sgn * disc) / (2 * qa)
            sols.append(mp.mpc(x, k1 * x + k0))
    else:
        k1, k0 = -C / B, -D / B
        qa = 1 + k1 ** 2
        qb = 2 * k1 * k0 + c + b * k1
        qc = k0 ** 2 + b * k0 + d
        disc = qb ** 2 - 4 * qa * qc
        if disc < -mp.mpf(10) ** -24:
            return []
        disc = mp.sqrt(max(disc, 0))
        for sgn in (1, -1):
            y = (-qb + sgn * disc) / (2 * qa)
            sols.append(mp.mpc(k1 * y + k0, y))
    return sols


def circle_intersections(c1, c2):
    c1, c2 = circ_norm(c1), circ_norm(c2)
    if c1[0] == 0 and c2[0] == 0:
        _, B1, C1, D1 = c1
        _, B2, C2, D2 = c2
        det = B1 * C2 - B2 * C1
        if abs(det) < mp.mpf(10) ** -20:
            return []
        x = (-D1 * C2 + D2 * C1) / det
        y = (-B1 * D2 + B2 * D1) / det
        return [mp.mpc(x, y)]
    if c1[0] == 0:
        return line_circle(c1, c2)
    if c2[0] == 0:
        return line_circle(c2, c1)
    line = (mp.mpf(0), c1[1] - c2[1], c1[2] - c2[2], c1[3] - c2[3])
    if max(abs(line[1]), abs(line[2])) < mp.mpf(10) ** -20:
        return []
    return line_circle(line, c1)


class Isom:
    """a side pairing: numeric and exact matrices, word"""
    def __init__(self, G, word):
        self.word = word
        self.Mn = G.word_num(word)
        self.Mx = G.F.word_matrix(word)
        self.rev = mdet(self.Mn) < 0
        M2 = mmul(self.Mn, self.Mn) if self.rev else self.Mn
        M2 = mnorm(M2)
        tr = abs(M2[0][0] + M2[1][1])
        T2 = 2 * mp.acosh(tr / 2)          # translation length of M2 (conformal, hyperbolic)
        self.T = T2 / 2 if self.rev else T2
        # axis endpoints: fixed points of M2
        A, B, C, D = M2[0][0], M2[0][1], M2[1][0], M2[1][1]
        if abs(C) < mp.mpf(10) ** -25:
            self.axis = (B / (D - A), None)      # second endpoint infinity
        else:
            disc = mp.sqrt((D - A) ** 2 + 4 * B * C)
            t1 = (-(D - A) - disc) / (2 * C)
            t2 = (-(D - A) + disc) / (2 * C)
            self.axis = (min(t1, t2), max(t1, t2))

    def ed_cycle(self, d):
        """list of general circles forming {z : d(z, s z) = d}"""
        if self.rev:
            cr = mp.cosh(d / 2) / mp.cosh(self.T / 2)
        else:
            if d < self.T - mp.mpf(10) ** -12:
                return []
            cr = mp.sinh(d / 2) / mp.sinh(self.T / 2)
        if cr < 1 - mp.mpf(10) ** -12:
            return []
        rho = mp.acosh(max(cr, 1))
        sh = mp.sinh(rho)
        p1, p2 = self.axis
        out = []
        if sh < mp.mpf(10) ** -12:
            # degenerate ED-cycle: d equals the translation length, the cycle is the axis itself.
            # (a residual sh ~ 1e-20 would displace the geodesic and lose tangential intersections)
            sh = mp.mpf(0)
            signs = (1,)
        else:
            signs = (1, -1)
        if p2 is None:
            for sg in signs:
                out.append((mp.mpf(0), mp.mpf(1), -sg * sh, -p1))
        else:
            for sg in signs:
                w = mp.mpc(sg * sh, 1)
                q = (p2 * w + p1) / (w + 1)
                out.append(circle3(mp.mpc(p1, 0), mp.mpc(p2, 0), q))
        return out


class Surface:
    def __init__(self, G, action, label=''):
        self.G = G
        self.N = G.N
        self.label = label
        self.d = action['d']
        self.perm = {'a': action['a'], 'b': action['b'], 'c': action['c']}
        self.check_torsion_free()
        self.build_domain()

    # ---- coset action ----
    def apply(self, pt, word):
        for ch in word:
            pt = self.perm[ch][pt]
        return pt

    def check_torsion_free(self):
        N = self.N
        for w in ['a', 'b', 'c', 'ab'] + ['bc' * k for k in (1, 2)] + ['ca' * k for k in range(1, N)]:
            for p in range(self.d):
                assert self.apply(p, w) != p, "torsion: %s fixes a coset" % w
        self.orientable = all(self.apply(0, w) != 0 or len(w) % 2 == 0 for w in ['a', 'b', 'c'])
        # orientable iff K < Delta: no odd word fixes point 0 -> parity well defined on cosets
        par = [None] * self.d
        par[0] = 0
        from collections import deque
        q = deque([0])
        ok = True
        while q:
            p = q.popleft()
            for ch in 'abc':
                r = self.perm[ch][p]
                if par[r] is None:
                    par[r] = 1 - par[p]
                    q.append(r)
                elif par[r] != 1 - par[p]:
                    ok = False
        self.orientable = ok

    # ---- fundamental domain ----
    def build_domain(self):
        N = self.N
        # polygons = orbits of <a,c>
        orbit = [None] * self.d
        polys = []
        for p in range(self.d):
            if orbit[p] is None:
                idx = len(polys)
                members = {}
                for m in range(N):
                    for e in (0, 1):
                        pt = self.apply(p, 'ca' * m + 'a' * e)
                        orbit[pt] = idx
                        members[pt] = (m, e)
                polys.append({'root': p, 'members': members})
        self.k = len(polys)
        assert self.k * 2 * N == self.d
        # BFS tree from the polygon of point 0, g = word; triangle labels relative to g
        gword = {0: ''}
        tree = []
        from collections import deque
        q = deque([0])
        sp_raw = []
        visited_edges = set()
        while q:
            P = q.popleft()
            gP = gword[P]
            for kk in range(N):
                rot = 'ca' * kk + 'ab'
                pt = self.apply(self.apply(0, gP), rot)
                Q = orbit[pt]
                if Q not in gword:
                    gword[Q] = reduce_word(gP + rot)
                    # relabel members of Q relative to the new root
                    root = self.apply(0, gword[Q])
                    members = {}
                    for m in range(N):
                        for e in (0, 1):
                            members[self.apply(root, 'ca' * m + 'a' * e)] = (m, e)
                    polys[Q]['members'] = members
                    tree.append((P, kk, Q))
                    q.append(Q)
                else:
                    m, e = polys[Q]['members'][pt]
                    # s = gP rot a^e (ca)^-m gQ^-1
                    w = gP + rot + 'a' * e + ('ac' * m) + gword[Q][::-1]
                    w = reduce_word(w)
                    if proj_equal(self.G.word_num(w), self.G.I):
                        continue   # the tree edge seen from the other side
                    sp_raw.append((P, kk, Q, m, e, w))
        self.sp_raw = sp_raw
        self.polys = polys
        self.gword = gword
        self.tree = tree
        # every s must lie in K: fixes point 0
        for (P, kk, Q, m, e, w) in sp_raw:
            assert self.apply(0, w) == 0
        # remove inverse duplicates
        isos = []
        for (P, kk, Q, m, e, w) in sp_raw:
            I = Isom(self.G, w)
            dup = False
            for J in isos:
                if proj_equal(mmul(I.Mn, J.Mn), self.G.I):
                    dup = True
                    break
            if not dup:
                isos.append(I)
        assert len(isos) <= (N * self.k - 2 * (self.k - 1)) // 2, (len(isos), N, self.k)
        self.n_edge_pairs = (N * self.k - 2 * (self.k - 1)) // 2
        self.sp = isos
        i0 = mp.mpc(0, 1)
        for s in self.sp:
            s.d0 = hdist(i0, act(s.Mn, i0))

    # ---- search ----
    def analyse(self, verbose=True):
        G = self.G
        H = G.F.H
        i0 = mp.mpc(0, 1)
        # distance list long enough
        sps = sorted(self.sp, key=lambda s: s.d0)
        s1 = sps[0]
        s2 = None
        for s in sps[1:]:
            if not (abs(s.axis[0] - s1.axis[0]) < 1e-10 and (
                    (s.axis[1] is None and s1.axis[1] is None) or
                    (s.axis[1] is not None and s1.axis[1] is not None and abs(s.axis[1] - s1.axis[1]) < 1e-10))):
                s2 = s
                break
        assert s2 is not None
        Hm = H + MARGIN      # search disc slightly larger than the circumdisc (tangential intersections on the boundary)
        D = max(2 * Hm + s2.d0, 2 * Hm + 2 * G.F.R + mp.mpf('0.5')) + mp.mpf('0.1')
        self.D = D
        if not hasattr(G, '_centres') or G._centres[0] < D:
            G._centres = (D, *G.centres(D))
        _, centres, dl = G._centres
        dists = [d for d, _ in dl if d <= D]
        self.dists = dists
        self.centres = centres
        if verbose:
            print('[%s] N=%d k=%d orientable=%s  sidepairings=%d  D=%s  #dist=%d' % (
                self.label, self.N, self.k, self.orientable, len(self.sp), mp.nstr(D, 5), len(dists)))
        D1 = 2 * Hm + s1.d0 + mp.mpf('1e-6')
        D2 = 2 * Hm + s2.d0 + mp.mpf('1e-6')
        b1 = [(i, c) for i, d in enumerate(dists) if d <= D1 for c in s1.ed_cycle(d)]
        b2 = [(j, c) for j, d in enumerate(dists) if d <= D2 for c in s2.ed_cycle(d)]
        cands = []
        pairs = []                # for each point of cands, the pairs (i, j): it lies on b_{s1}(d_i) and b_{s2}(d_j)
        for i, c1 in b1:
            for j, c2 in b2:
                for z in circle_intersections(c1, c2):
                    if z.imag <= 0:
                        continue
                    if hdist(i0, z) > Hm:
                        continue
                    old = [m for m, w in enumerate(cands) if abs(z - w) < mp.mpf(10) ** -9]
                    if old:
                        if (i, j) not in pairs[old[0]]:
                            pairs[old[0]].append((i, j))
                        continue
                    cands.append(z)
                    pairs.append([(i, j)])
        if verbose:
            print('   ED-cycle intersections in the circumdisc:', len(cands))
        # filter by all generators
        good = []
        status = []               # for each point of cands: None (candidate) or (side pairing, its displacement)
        for z in cands:
            ok = True
            for s in self.sp:
                dz = hdist(z, act(s.Mn, z))
                if dz > dists[-1] - TOL:
                    continue          # beyond the computed list: undecided, keep
                if min(abs(dz - d) for d in dists) > TOL:
                    ok = False
                    status.append((s.word, dz))
                    break
            if ok:
                good.append(z)
                status.append(None)
        if verbose:
            print('   candidates passing the displacement test:', len(good))
        self.candidates = good
        self.intersection_points = cands      # all the points examined in Step 3 (used by margins.py)
        # what is needed to show the search (check_surface.py --details, draw_surface.py --edcycles)
        self.search = {'s1': s1, 's2': s2, 'radius': Hm, 'D1': D1, 'D2': D2,
                       'points': [{'z': z, 'pairs': pr, 'killed_by': st, 'packings': [], 'rejection': None}
                                  for z, pr, st in zip(cands, pairs, status)]}
        point_of = {id(z): pt for z, pt in zip(cands, self.search['points'])}
        # alpha search
        packings = []
        self.n_uncertified = 0
        self.rejected = []          # candidates with no alpha: (z, independently_confirmed)
        for z in good:
            alphas = self.find_alphas(z, s1)
            if not alphas:
                # independent confirmation that z is not a centre: an element kappa of K, other than the identity,
                # that displaces z a distance which is not admissible (a centre is displaced by every kappa != 1 an
                # admissible distance). We try the side pairings and the products of two of them; a displacement
                # beyond the list of admissible distances is decided with a longer list (Lemma 3.1 gives the
                # complete list up to any bound).
                best = self.confirm_not_centre(z, dists)
                conf = best is not None
                point_of[id(z)]['rejection'] = best
                self.rejected.append((z, conf))
            for alpha, words in alphas:
                res = self.certify(alpha, words)
                if res is None:
                    self.n_uncertified += 1
                    if verbose:
                        print('   !! numerical alpha could not be verified exactly for z=', mp.nstr(to_disc(z), 8))
                    continue
                alpha_x, lam = res
                key = self.packing_key(alpha)
                dup = None
                for P in packings:
                    if self.same_key(P['key'], key):
                        dup = P
                        break
                if dup is None:
                    dup = {'z': z, 'alpha': alpha, 'alpha_exact': alpha_x, 'words': words,
                           'key': key, 'det_sign': 1 if mdet(alpha) > 0 else -1}
                    packings.append(dup)
                if all(Q is not dup for Q in point_of[id(z)]['packings']):
                    point_of[id(z)]['packings'].append(dup)
        self.packings = packings
        if verbose:
            print('   extremal packings found (verified exactly): %d' % len(packings))
            for P in packings:
                print('      centre in disc:', mp.nstr(to_disc(P['z']), 10), ' d(0,z)=', mp.nstr(hdist(i0, P['z']), 10))
        return packings

    def confirm_not_centre(self, z, dists, cap=10):
        """an element of K (a side pairing or a product of two) that displaces z by a non admissible distance,
        with the largest margin: (w1, w2, inv, displacement, gap, bound of the list used), meaning the element
        w1 w2^(+-1) (w2 = None for a single side pairing), or None if there is none"""
        G = self.G
        elems = [(s.word, None, False, s.Mn) for s in self.sp]
        for s in self.sp:
            for t in self.sp:
                elems.append((s.word, t.word, False, mmul(s.Mn, t.Mn)))
                elems.append((s.word, t.word, True, mmul(s.Mn, minv(t.Mn))))
        disp = []
        for w1, w2, inv, M in elems:
            dz = hdist(z, act(M, z))
            if dz > mp.mpf(10) ** -6:                 # the identity (e.g. s s^-1) is excluded
                disp.append((w1, w2, inv, dz))
        bound, dl = dists[-1], dists
        while True:
            best = None
            for w1, w2, inv, dz in disp:
                if dz < bound - TOL:
                    gap = min(abs(dz - d) for d in dl)
                    if gap > TOL and (best is None or gap > best[4]):
                        best = (w1, w2, inv, dz, gap, bound)
            if best is not None or bound >= cap:
                return best
            # extend the list of admissible distances just beyond the smallest undecided displacement
            todo = [dz for _, _, _, dz in disp if dz >= bound - TOL]
            if not todo:
                return None
            newb = min(min(todo) + mp.mpf('0.05'), mp.mpf(cap))
            if not hasattr(G, '_ext_dists') or G._ext_dists[0] < newb:
                G._ext_dists = (newb, [d for d, _ in G.centres(newb)[1]])
            bound = newb
            dl = [d for d in G._ext_dists[1] if d <= bound]

    def find_alphas(self, z, s):
        """numeric alphas (with the reduced words of all generators) such that alpha^-1 K alpha < Delta"""
        G = self.G
        i0 = mp.mpc(0, 1)
        x, y = z.real, z.imag
        sy = mp.sqrt(y)
        beta = [[sy, x / sy], [mp.mpf(0), 1 / sy]]
        betai = minv(beta)
        p = act(betai, act(s.Mn, z))
        dd = hdist(i0, p)
        targets = [c for (d, c, w) in self.centres if abs(d - dd) < TOL]
        thetas = []
        argp = mp.arg(to_disc(p))
        for c in targets:
            th = argp - mp.arg(to_disc(c))
            step = 2 * mp.pi / G.N
            th = th - mp.floor(th / step) * step
            if not any(abs(th - t) < 1e-9 or abs(abs(th - t) - step) < 1e-9 for t in thetas):
                thetas.append(th)
        out = []
        for th in thetas:
            for sg in (1, -1):
                R = [[mp.cos(sg * th / 2), mp.sin(sg * th / 2)], [-mp.sin(sg * th / 2), mp.cos(sg * th / 2)]]
                alpha = mmul(beta, R)
                ai = minv(alpha)
                words = []
                okall = True
                for t in self.sp:
                    h = mmul(ai, mmul(t.Mn, alpha))
                    w = G.reduce(h)
                    if w is None:
                        okall = False
                        break
                    words.append(w)
                if okall:
                    out.append((alpha, words))
                    # The two signs are the rotations by +theta and -theta about O (the sign convention of R is
                    # not used anywhere else). Only one of them maps the target centre c to p; if the rotation by
                    # -theta also defines a packing, that packing is found from its own target centre, whose angle
                    # -theta is also in the list thetas. So we can stop at the first sign that works.
                    break
        return out

    def certify(self, alpha_num, words):
        """exact: find alpha in M_2(L) with s_j alpha = +-alpha w_j for all j. Returns (alpha_exact, signs)"""
        F = self.G.F
        W = [F.word_matrix(w) for w in words]
        S = [s.Mx for s in self.sp]
        # linear map X -> S0 X W0^-1 on the 4-dim space
        def lin(Sm, Wm):
            Wi = Wm ** (-1)
            cols = []
            for (i, j) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                E = pari.matrix(2, 2, [F.one if (r, c) == (i, j) else F.zero for r in range(2) for c in range(2)])
                X = Sm * E * Wi
                cols.append([X[0, 0], X[0, 1], X[1, 0], X[1, 1]])
            # matrix whose columns are the images
            return pari.matrix(4, 4, [cols[c][r] for r in range(4) for c in range(4)])
        Id = pari.matrix(4, 4, [F.one if r == c else F.zero for r in range(4) for c in range(4)])
        # signs from the numeric alpha
        blocks = []
        signs = []
        for Sm, Wm, s, w in zip(S, W, self.sp, words):
            Xn = mmul(s.Mn, mmul(alpha_num, minv(self.G.word_num(w))))
            # Xn should be +- alpha_num up to the scaling of word_num (|det| = 1 for both)
            Xn = mnorm(Xn)
            An = mnorm(alpha_num)
            e1 = max(abs(Xn[i][j] - An[i][j]) for i in range(2) for j in range(2))
            e2 = max(abs(Xn[i][j] + An[i][j]) for i in range(2) for j in range(2))
            if min(e1, e2) > mp.mpf(10) ** -8:
                return None
            lam = F.one if e1 < e2 else -F.one
            signs.append(1 if e1 < e2 else -1)
            blocks.append(lin(Sm, Wm) - lam * Id)
        n = len(blocks)
        big = pari.matrix(4 * n, 4, [blocks[b][r, c] for b in range(n) for r in range(4) for c in range(4)])
        K = pari.matker(big)
        ncol = int(pari.matsize(K)[1])
        if ncol != 1:
            print('   !! kernel dimension', ncol)
            return None
        v = [K[r, 0] for r in range(4)]
        alpha = pari.matrix(2, 2, v)
        an = F.mat_num(alpha)
        if not proj_equal(an, alpha_num, tol=mp.mpf(10) ** -8):
            print('   !! exact alpha does not match the numeric one')
            return None
        for Sm, Wm, sg in zip(S, W, signs):
            X = Sm * alpha * Wm ** (-1)
            assert X == alpha if sg == 1 else X == -alpha
        if F.det(alpha) == 0:
            return None
        return alpha, signs

    def packing_key(self, alpha):
        """set of centres of the packing alpha Delta(i) within distance H + d0 + 0.3 of i (disc coords)"""
        G = self.G
        i0 = mp.mpc(0, 1)
        rad = G.F.H + 2 * G.F.R + mp.mpf('0.3')
        w0 = act(minv(alpha), i0)
        pts = [to_disc(act(alpha, c)) for (d, c, w) in self.centres if hdist(w0, c) <= rad]
        return pts

    @staticmethod
    def same_key(k1, k2):
        if len(k1) != len(k2):
            return False
        for p in k1:
            if not any(abs(p - q) < mp.mpf(10) ** -8 for q in k2):
                return False
        return True
