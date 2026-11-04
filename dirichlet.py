"""Dirichlet domains of the surface groups, centred at the centre of a disc of an extremal packing.

For a torsion-free K < Delta^pm(2,3,N) (given by its permutation action on the cosets) the Dirichlet domain
    D_O(K) = { z : rho(z, O) <= rho(z, kappa O) for every kappa in K }
is computed as follows.
  * The orbit K.O consists of centres g(O) of the tessellation {N,3}; g(O) belongs to it if and only if g r lies in K
    for some r in the stabiliser <a, c> of O, a condition on the permutation representation. This gives the orbit
    points together with explicit words kappa = g r in K (see orbit_points).
  * In the Klein model the bisector of O and kappa(O) is a straight line, so D_O(K) is the intersection of half
    planes. It is first computed in floating point with all the orbit points at distance at most D from O; the lines
    that are active or nearly active (within 1e-6) are then recomputed from their words in 40 digit arithmetic, and
    the final polygon, its side pairings and the checks below are obtained in 40 digits.
  * Completeness: if rho_max is the largest distance from O to a vertex of the computed polygon, every orbit point
    whose bisector meets the polygon is at distance at most 2 rho_max from O; the computation is repeated with
    D >= 2 rho_max, so no side is missing.
  * Side pairings: the side on the bisector of O and kappa(O) is mapped by kappa^-1 onto the side on the bisector of
    O and kappa^-1(O).
  * Checks: every side is paired; every vertex cycle has angle sum 2 pi; the area equals 2 pi |chi(X)| (Gauss-Bonnet).
    By Poincare's theorem the side pairings then generate K and the polygon is a fundamental domain.
For the other extremal packings the same is done for K_j = alpha_j^-1 K alpha_j (whose coset action is computed by
GAP from the words w_s), and the result is mapped by alpha_j: this is the Dirichlet domain of K centred at
alpha_j(O), the centre of a disc of the j-th packing.
"""
import mpmath as mp
from geom import act, mmul, minv, mdet

I0 = mp.mpc(0, 1)


def to_klein(z):
    """upper half plane -> Klein disc (through the Poincare disc of geom.to_disc)"""
    w = mp.mpc(0, 1) * (z - I0) / (z + I0)
    return 2 * w / (1 + abs(w) ** 2)


def klein_to_poincare(k):
    return k / (1 + mp.sqrt(1 - abs(k) ** 2))


def lift(k):
    """point of the hyperboloid (B(x,x) = 1, B = x0 y0 - x1 y1 - x2 y2) over the Klein point k"""
    s = 1 / mp.sqrt(1 - abs(k) ** 2)
    return (s, s * k.real, s * k.imag)


def B(x, y):
    return x[0] * y[0] - x[1] * y[1] - x[2] * y[2]


def kdist(k1, k2):
    c = B(lift(k1), lift(k2))
    return mp.acosh(max(c, mp.mpf(1)))


def orbit_points(G, perm, d, D):
    """centres kappa(O), kappa in K, with rho(O, kappa O) <= D, in floating point: list of (Klein point, word).

    The orbit K.O consists of centres g(O) of the tessellation {N,3}; g(O) belongs to it if and only if g r is in K
    for some r in the stabiliser <a, c> of O, i.e. if the coset K g is one of the cosets K r^-1 (a condition on the
    permutation representation, followed along the search). The centres are enumerated by a breadth first search
    from polygon to adjacent polygon, keeping those at distance at most D + r_c from O (as in Lemma 3.1 of the
    paper this reaches every centre at distance at most D)."""
    import math
    stab = ['ca' * m + e for m in range(G.N) for e in ('', 'a')]

    def apply(p, w):
        for ch in w:
            p = perm[ch][p]
        return p
    target = {}
    for r in stab:
        target.setdefault(apply(0, r[::-1]), r)          # K g r = K  iff  K g = K r^-1
    steps = []
    for kk in range(G.N):
        R = G.Rota[kk]
        steps.append(('ca' * kk + 'ab', (complex(R[0][0]), complex(R[0][1]), complex(R[1][0]), complex(R[1][1]))))
    lim = float(D + G.F.H) + 0.01
    Df = float(D)
    # queue entries: (matrix a, b, c, d (float), coset point, parent index, step index)
    queue = [((1 + 0j, 0j, 0j, 1 + 0j), 0, -1, -1)]
    seen = {(0, 0)}
    found = []
    head = 0
    while head < len(queue):
        (a, b, c, dd), p, par, st = queue[head]
        if head and p in target:
            det = (a * dd - b * c).real
            z = -1j if det < 0 else 1j
            zz = (a * z + b) / (c * z + dd)
            dist = math.acosh(1 + abs(zz - 1j) ** 2 / (2 * zz.imag))
            if dist <= Df + 1e-9:
                wd = 1j * (zz - 1j) / (zz + 1j)
                found.append((2 * wd / (1 + abs(wd) ** 2), head, target[p]))
        for si, (ws, (ra, rb, rc, rd)) in enumerate(steps):
            m2 = (a * ra + b * rc, a * rb + b * rd, c * ra + dd * rc, c * rb + dd * rd)
            det = (m2[0] * m2[3] - m2[1] * m2[2]).real
            z = -1j if det < 0 else 1j
            z2 = (m2[0] * z + m2[1]) / (m2[2] * z + m2[3])
            wd = 1j * (z2 - 1j) / (z2 + 1j)
            key = (round(wd.real * 1e11), round(wd.imag * 1e11))
            if key in seen:
                continue
            seen.add(key)
            if math.acosh(1 + abs(z2 - 1j) ** 2 / (2 * z2.imag)) > lim:
                continue
            queue.append((m2, apply(p, ws), head, si))
        head += 1

    def word(i):
        parts = []
        while queue[i][2] >= 0:
            parts.append(steps[queue[i][3]][0])
            i = queue[i][2]
        return ''.join(reversed(parts))
    return [(k, word(i) + r) for k, i, r in found]


def fclip(poly, line):
    """floating point version of clip (vertices complex, labels arbitrary)"""
    c0, c1, c2, lab = line
    out = []
    n = len(poly)
    for i in range(n):
        v, e = poly[i]
        w, _ = poly[(i + 1) % n]
        fv = c0 + c1 * v.real + c2 * v.imag
        fw = c0 + c1 * w.real + c2 * w.imag
        if fv >= 0:
            out.append((v, e))
            if fw < 0:
                out.append((v + fv / (fv - fw) * (w - v), lab))
        elif fw >= 0:
            out.append((v + fv / (fv - fw) * (w - v), e))
    return out


def clip(poly, line):
    """poly: list of (vertex, label of the edge from this vertex to the next); keep {(1, k) . line >= 0}"""
    c0, c1, c2, lab = line
    f = lambda v: c0 + c1 * v.real + c2 * v.imag
    out = []
    n = len(poly)
    for i in range(n):
        v, e = poly[i]
        w, _ = poly[(i + 1) % n]
        fv, fw = f(v), f(w)
        if fv >= 0:
            out.append((v, e))
            if fw < 0:
                t = fv / (fv - fw)
                out.append((v + t * (w - v), lab))
        elif fw >= 0:
            t = fv / (fv - fw)
            out.append((v + t * (w - v), e))
    # remove (numerically) repeated vertices
    clean = []
    for v, e in out:
        if clean and abs(v - clean[-1][0]) < mp.mpf(10) ** -25:
            clean[-1] = (clean[-1][0], e)
            continue
        clean.append((v, e))
    if len(clean) > 1 and abs(clean[0][0] - clean[-1][0]) < mp.mpf(10) ** -25:
        clean.pop()
    return clean


def angle_at(u, v, w):
    """interior angle at v of the hyperbolic triangle u v w (Klein points)"""
    a, b, c = kdist(v, u), kdist(v, w), kdist(u, w)
    x = (mp.cosh(a) * mp.cosh(b) - mp.cosh(c)) / (mp.sinh(a) * mp.sinh(b))
    return mp.acos(max(-1, min(1, x)))


def kact(M, k):
    """isometry M (matrix, upper half plane) applied to a Klein point"""
    p = klein_to_poincare(k)
    z = (1 - mp.mpc(0, 1) * p) / (p - mp.mpc(0, 1))          # Poincare disc -> upper half plane
    return to_klein(act(M, z))


def dirichlet_at_O(G, perm, d, chi, D=mp.mpf(4)):
    """Dirichlet domain of the group with permutation action perm (point 0 = K) centred at O"""
    import math
    while True:
        pts = orbit_points(G, perm, d, D)
        flines = []
        for k, w in pts:
            s_ = 1 / math.sqrt(1 - abs(k) ** 2)
            flines.append((s_ - 1, -s_ * k.real, -s_ * k.imag, w))
        poly = [(complex(-1, -1), None), (complex(1, -1), None), (complex(1, 1), None), (complex(-1, 1), None)]
        for ln in sorted(flines, key=lambda l: l[0]):
            poly = fclip(poly, ln)
        bounded = poly and all(e is not None for _, e in poly) and all(abs(v) < 1 - 1e-9 for v, _ in poly)
        if bounded:
            rho = max(math.acosh(1 / math.sqrt(1 - abs(v) ** 2)) for v, _ in poly)
            if 2 * rho + 1e-6 <= D:
                break
            D = mp.mpf(2 * rho + 0.2)
        else:
            D = D + 2
    # exact stage: the lines that are active or nearly active on the floating point polygon, in 40 digits.
    # A line with c0 + c1 x + c2 y >= 1e-6 * |(c1, c2)| at every vertex of that polygon cannot cut the domain
    # (the floating point errors are far below 1e-6), so only the others are recomputed.
    keep = []
    for c0, c1, c2, w in flines:
        sc = math.hypot(c1, c2) + abs(c0)
        if min(c0 + c1 * v.real + c2 * v.imag for v, _ in poly) <= 1e-6 * sc:
            keep.append(w)
    Ok = lift(mp.mpc(0))
    lines = []
    for w in keep:
        M = G.word_num(w)
        q = lift(to_klein(act(M, I0)))
        # B(x, q - O) >= 0  with x = (1, k): (q0 - 1) - q1 k1 - q2 k2 >= 0
        lines.append((q[0] - Ok[0], -q[1], -q[2], (w, M)))
    poly = [(mp.mpc(-1, -1), None), (mp.mpc(1, -1), None), (mp.mpc(1, 1), None), (mp.mpc(-1, 1), None)]
    for ln in sorted(lines, key=lambda l: l[0]):
        poly = clip(poly, ln)
    rho = max(kdist(mp.mpc(0), v) for v, _ in poly)
    m = len(poly)
    verts = [v for v, _ in poly]
    words = [e[0] for _, e in poly]                    # side i: from verts[i] to verts[i+1], on bisector of O, kappa O
    mats = [e[1] for _, e in poly]
    # pairing: side i (kappa) <-> side j whose kappa_j = kappa_i^-1 (same orbit point kappa^-1 O)
    target = [to_klein(act(minv(M), I0)) for M in mats]
    own = [to_klein(act(M, I0)) for M in mats]
    partner = []
    for i in range(m):
        js = [j for j in range(m) if abs(own[j] - target[i]) < mp.mpf(10) ** -20]
        partner.append(js[0] if len(js) == 1 else None)
    angles = [angle_at(verts[i - 1], verts[i], verts[(i + 1) % m]) for i in range(m)]
    # vertex cycles: vertex i lies on sides i-1 and i; the pairing of side i maps it to a vertex of the partner side
    cycles, used = [], set()
    for i0 in range(m):
        if i0 in used:
            continue
        cyc, i, side = [], i0, i0
        for _ in range(4 * m):
            if i in used and i == i0 and cyc:
                break
            used.add(i)
            cyc.append(i)
            j = partner[side]
            if j is None:
                break
            img = kact(minv(mats[side]), verts[i])     # kappa^-1 maps side 'side' onto side j
            cand = [t for t in (j, (j + 1) % m) if abs(verts[t] - img) < mp.mpf(10) ** -15]
            if not cand:
                break
            i = cand[0]
            side = (i - 1) % m if j == i else i        # the other side at the new vertex
            if i == i0:
                break
        cycles.append(cyc)
    cycle_sums = [sum(angles[i] for i in c) for c in cycles]
    area = (m - 2) * mp.pi - sum(angles)
    checks = {
        'sides': m,
        'all_paired': all(p is not None for p in partner) and all(partner[partner[i]] == i for i in range(m)),
        'in_K': True,                                   # by construction (words fixing the coset K)
        'cycle_sums_ok': all(abs(s - 2 * mp.pi) < mp.mpf(10) ** -20 for s in cycle_sums),
        'area_ok': abs(area - 2 * mp.pi * abs(chi)) < mp.mpf(10) ** -20,
        'euler_ok': len(cycles) - m // 2 + 1 == chi,
        'cycles': [len(c) for c in cycles],
        'rho_max': rho,
        'D': D,
    }
    checks['ok'] = checks['all_paired'] and checks['cycle_sums_ok'] and checks['area_ok'] and checks['euler_ok']
    return {'verts': verts, 'words': words, 'mats': mats, 'partner': partner, 'angles': angles,
            'cycles': cycles, 'area': area, 'checks': checks}


def transport(dom, alpha):
    """the Dirichlet domain of K_j at O mapped by alpha (Dirichlet domain of K at alpha(O))"""
    verts = [kact(alpha, v) for v in dom['verts']]
    mats = [mmul(mmul(alpha, M), minv(alpha)) for M in dom['mats']]
    out = dict(dom)
    out['verts'] = verts
    out['mats'] = mats
    out['centre'] = to_klein(act(alpha, I0))
    return out
