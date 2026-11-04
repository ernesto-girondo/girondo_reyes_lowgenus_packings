"""Figures for the article (Poincare disc, matplotlib)."""
import json, sys
import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from geom import Geometry, to_disc, hdist, act, mmul, minv
from surface import Surface, Isom
from gapiface import coset_action_from_words, reduce_word

mp.mp.dps = 40
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures') + os.sep   # the figures of the paper
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'cm'})


def cplx(z):
    return complex(float(z.real), float(z.imag))


def disc(z):
    return cplx(to_disc(z))


def geodesic_pts(z1, z2, n=40):
    """points of the hyperbolic segment between two points of the disc"""
    z1, z2 = complex(z1), complex(z2)
    # move z1 to 0
    w = (z2 - z1) / (1 - np.conj(z1) * z2)
    ts = np.linspace(0, 1, n)
    pts = []
    for t in ts:
        u = t * w
        pts.append((u + z1) / (1 + np.conj(z1) * u))
    return np.array(pts)


def draw_segment(ax, z1, z2, **kw):
    p = geodesic_pts(z1, z2)
    ax.plot(p.real, p.imag, **kw)


def polygon_vertices(G, M=None):
    """vertices (disc) of the image of the central polygon under the numeric matrix M"""
    N = G.N
    vs = []
    for j in range(N):
        R = G.word_num('ca' * j)
        v = act(R, G.vA)
        if M is not None:
            v = act(M, v)
        vs.append(disc(v))
    # order by argument around the centre
    return vs


def draw_polygon(ax, G, M=None, **kw):
    vs = polygon_vertices(G, M)
    c = disc(act(M, mp.mpc(0, 1))) if M is not None else 0j
    # sort vertices by angle around centre c (hyperbolic polygon is convex, star-shaped w.r.t. c)
    vs = sorted(vs, key=lambda v: np.angle(v - c))
    for i in range(len(vs)):
        draw_segment(ax, vs[i], vs[(i + 1) % len(vs)], **kw)


def setup(ax, r=1.0):
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=0.6, color='k'))
    ax.set_xlim(-r, r); ax.set_ylim(-r, r)
    ax.set_aspect('equal'); ax.axis('off')


# ---------------------------------------------------------------- data of X0 (the 18-gon of Figure 1)
def X0_data():
    G18 = Geometry(18)
    g = [json.loads(l) for l in open('groups/fuchs_N18.jsonl') if l.strip()][8]
    action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
    S = Surface(G18, action, label='X0')
    S.analyse(verbose=False)
    return G18, S


def fig_X0_domain():
    G, S = X0_data()
    N = 18
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    setup(ax, 0.95)
    O = mp.mpc(0, 1)
    # inscribed (extremal) disc: Euclidean circle of radius tanh(R/2)
    ax.add_patch(Circle((0, 0), float(mp.tanh(G.F.R / 2)), fill=True, color='#dbe9f6', lw=0))
    draw_polygon(ax, G, None, color='k', lw=1.1)
    # side pairings: edge kk of the polygon paired with edge m; label both with the same number
    pairs = {}
    label = 0
    seen = set()
    for (P, kk, Q, m, e, w) in S.sp_raw:
        key = tuple(sorted((kk, m)))
        if key in seen:
            continue
        seen.add(key)
        label += 1
        pairs[kk] = label; pairs[m] = label
    # midpoints of edges: (ca)^j (C) ; edge j of the central polygon has midpoint (ca)^j(i/u)
    C = mp.mpc(0, 1 / G.u)
    for j in range(N):
        z = disc(act(G.word_num('ca' * j), C))
        zz = z * 1.09
        ax.text(zz.real, zz.imag, str(pairs.get(j, '?')), ha='center', va='center', fontsize=7.5)
    # centres of the extremal discs
    names = iter(['u', 'v', 'w'])
    for P in S.packings:
        z = disc(P['z'])
        ax.plot([z.real], [z.imag], 'o', color='C3', ms=5)
        if abs(z) > 0.1:
            ax.text(z.real * 1.12, z.imag * 1.12, '$%s$' % next(names), color='C3', ha='center', va='center')
    # Weierstrass points (from the earlier computation: distances/angles in this frame)
    wp = [(0.927539, 60), (0.927539, 180), (0.927539, -60), (1.445482, 0), (1.445482, 120), (1.445482, -120)]
    for d, ang in wp:
        r = float(mp.tanh(mp.mpf(d) / 2))
        z = r * np.exp(1j * np.radians(ang))
        ax.plot([z.real], [z.imag], 's', color='C0', ms=4)
    ax.plot([0], [0], 'o', color='C3', ms=5)
    ax.text(0.03, 0.03, r'$O$', color='C3')
    fig.savefig(OUT + 'fig_X0_domain.pdf', bbox_inches='tight')
    plt.close(fig)
    return pairs


def fig_edcycles():
    G, S = X0_data()
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    rc = float(mp.tanh(G.F.H / 2))
    setup(ax, rc * 1.15)
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=0))
    draw_polygon(ax, G, None, color='0.6', lw=0.8)
    ax.add_patch(Circle((0, 0), rc, fill=False, lw=0.8, ls='--', color='k'))
    sps = sorted(S.sp, key=lambda s: s.d0)
    s1 = sps[0]
    s2 = None
    for s in sps[1:]:
        if abs(s.axis[0] - s1.axis[0]) > 1e-8:
            s2 = s; break
    dists = S.dists
    H = G.F.H
    cols = ['C2', 'C1']
    def hypercycles(s, d):
        """points (disc) of the two hypercycles {rho(z, s z) = d}, by moving the axis to standard position"""
        if s.rev:
            cr = mp.cosh(d / 2) / mp.cosh(s.T / 2)
        else:
            if d < s.T - mp.mpf(10) ** -12:
                return []
            cr = mp.sinh(d / 2) / mp.sinh(s.T / 2)
        if cr < 1 - mp.mpf(10) ** -12:
            return []
        rho = mp.acosh(max(cr, 1))
        p1, p2 = s.axis
        p1 = float(p1)
        t = np.linspace(-9, 9, 1500)
        out = []
        for sg in (1, -1):
            xs = sg * float(mp.sinh(rho)) * np.exp(t); ys = np.exp(t)
            z = xs + 1j * ys
            if p2 is None:
                w = z + p1
            else:
                p2f = float(p2)
                w = (p2f * z + p1) / (z + 1)      # phi: 0 -> p1, inf -> p2
            out.append(1j * (w - 1j) / (w + 1j))
        return out

    for s, col in ((s1, cols[0]), (s2, cols[1])):
        Dmax = 2 * (H + mp.mpf('0.1')) + s.d0
        for d in dists:
            if d > Dmax: break
            for w in hypercycles(s, d):
                if np.min(np.abs(w)) > rc + 0.02:
                    continue
                ax.plot(w.real, w.imag, color=col, lw=0.7, alpha=0.9)
    for z in S.candidates:
        zz = disc(z)
        ax.plot([zz.real], [zz.imag], 'o', color='k', ms=4, mfc='w')
    for P in S.packings:
        zz = disc(P['z'])
        ax.plot([zz.real], [zz.imag], 'o', color='C3', ms=5)
    fig.savefig(OUT + 'fig_edcycles.pdf', bbox_inches='tight')
    plt.close(fig)
    return s1.word, s2.word


def fig_two_tessellations():
    """the 18-gon of X0 with the nonagons of the Aut-invariant 4-packing and of one of the other four"""
    src = open('cross_check_9_18.py').read().split("from gapiface import")[0].replace("S.analyse(verbose=True)", "S.analyse(verbose=False)")
    ns = {}
    exec(src, ns)
    S = ns['S']; G9 = ns['G9']; G18 = ns['G18']
    res9, _ = G9.centres(mp.mpf(4.5))
    # identify which packing is the invariant one: its K_j class is #37 -> use canonical keys
    from batch import canonical_key
    res = json.load(open('groups/fuchs_N9_results.json'))
    key37 = tuple(map(tuple, res[36]['key']))
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2))
    chosen = []
    for P in S.packings:
        if G9.reduce(P['alpha']) is not None:
            continue
        a2 = coset_action_from_words(9, [reduce_word(w) for w in P['words']])
        k = canonical_key({'a': a2['a'], 'b': a2['b'], 'c': a2['c']}, a2['d'])
        chosen.append((k == key37, P))
    inv = [P for f, P in chosen if f][0]
    oth = [P for f, P in chosen if not f and hdist(mp.mpc(0, 1), P['z']) < 1e-6]
    oth = oth[0] if oth else [P for f, P in chosen if not f][0]
    for ax, P, title in ((axes[0], inv, 'invariant 4-packing'), (axes[1], oth, 'one of the other four')):
        setup(ax, 0.97)
        ax.add_patch(Circle((0, 0), float(mp.tanh(G18.F.R / 2)), fill=True, color='#dbe9f6', lw=0))
        draw_polygon(ax, G18, None, color='k', lw=1.1)
        al = P['alpha']
        for d, z, w in res9:
            M = mmul(al, G9.word_num(w))
            c = disc(act(M, mp.mpc(0, 1)))
            if abs(c) > 0.93:
                continue
            draw_polygon(ax, G9, M, color='C3', lw=0.6)
            ax.plot([c.real], [c.imag], 'o', color='C3', ms=2.5)
        ax.set_title(title, fontsize=9)
    fig.savefig(OUT + 'fig_two_tessellations.pdf', bbox_inches='tight')
    plt.close(fig)


def fig_dirichlet_X0():
    """Dirichlet domains of X_0 (frame of the 37th class of index 36 in Delta(2,3,9)) centred at O, a disc of the
    invariant 4-packing, and at a disc of one of the other four packings (Section 7 of the paper)"""
    import draw_surface as ds
    G, S, label = ds.load('fuchs_N9', 37, False)
    S.packings.sort(key=lambda P: G.reduce(P['alpha']) is None)
    polyM = {P: G.word_num(S.gword[P]) for P in sorted(S.gword)}
    doms, chi = ds.dirichlet_domains(G, S, 'all')
    # the invariant packing (index 0) and the first packing of the orbit whose chosen disc is not centred at O
    other = [(i, d) for i, d in doms[1:] if abs(complex(d['centre'])) > 1e-6][0]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6))
    exts = []
    for (idx, dom), ax in zip([doms[0], other], axes):
        exts.append(ds.draw_dirichlet_panel(ax, G, S, polyM, idx, dom))
    r = min(1.0, max(exts) + 0.05)
    for ax, t in zip(axes, ['centred at a disc of the invariant packing', 'centred at a disc of another packing']):
        ax.set_xlim(-r, r); ax.set_ylim(-r, r)
        ax.set_title(t, fontsize=9)
    fig.savefig(OUT + 'fig_dirichlet_X0.pdf', bbox_inches='tight')
    plt.close(fig)
    return [(i, d['checks']['sides'], d['checks']['cycles']) for i, d in (doms[0], other)]


def fig_classical_curves():
    """Remark 6.4: Bolza's curve (fuchs_N8#77, 6-extremal) and y^2 = x^5 - 1 (fuchs_N10#18, 3-extremal). For each,
    the Dirichlet domain centred at a disc of the extremal packing, with the tessellation by regular N-gons and the
    discs of the packing inscribed in them, clipped to the domain."""
    import draw_surface as ds
    from dirichlet import klein_to_poincare
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.7))
    out = []
    for ax, (fam, n, title) in zip(axes, [('fuchs_N8', 77, r'$y^2=x^5-x$ (Bolza): six octagons'),
                                          ('fuchs_N10', 18, r'$y^2=x^5-1$: three decagons')]):
        G, S, label = ds.load(fam, n, False)
        doms, chi = ds.dirichlet_domains(G, S, 'O')
        dom = doms[0][1]
        ax.set_aspect('equal'); ax.axis('off')
        ax.add_patch(Circle((0, 0), 1, fill=False, lw=0.6, color='k'))
        V = dom['verts']; m = len(V)
        bnd, ext = [], 0
        for i in range(m):
            k1, k2 = V[i], V[(i + 1) % m]
            bnd += [complex(klein_to_poincare(k1 + t * (k2 - k1))) for t in np.linspace(0, 1, 60)[:-1]]
        ext = max(abs(z) for z in bnd)
        clip = PathPatch(Path([(z.real, z.imag) for z in bnd] + [(bnd[0].real, bnd[0].imag)], closed=True),
                         transform=ax.transData, fc='none', ec='none')
        ax.add_patch(clip)
        s = float(mp.tanh(G.F.R / 2))
        cents = G.centres(2 * float(mp.atanh(ext)) + 2 * float(G.F.H))[0]
        for d, z, w in cents:
            wd = disc(z)
            # inscribed disc of the polygon centred at wd: Euclidean circle
            a = abs(wd) ** 2
            c = wd * (1 - s * s) / (1 - a * s * s)
            r = s * (1 - a) / (1 - a * s * s)
            circ = Circle((c.real, c.imag), r, fc='#dbe9f6', ec='#8fb3d9', lw=0.5)
            ax.add_patch(circ); circ.set_clip_path(clip)
            M = G.word_num(w)
            vs = ds.polygon_disc_vertices(G, M)
            vs.sort(key=lambda v: np.angle((v - wd) / (1 - np.conj(wd) * v)))
            for i in range(len(vs)):
                pts = geodesic_pts(vs[i], vs[(i + 1) % len(vs)])
                ln, = ax.plot(pts.real, pts.imag, color='0.55', lw=0.5)
                ln.set_clip_path(clip)
        lab, num = {}, 0
        for i in range(m):
            if i in lab:
                continue
            num += 1
            lab[i] = num
            if dom['partner'][i] is not None:
                lab[dom['partner'][i]] = num
        for i in range(m):
            k1, k2 = V[i], V[(i + 1) % m]
            pts = [complex(klein_to_poincare(k1 + t * (k2 - k1))) for t in np.linspace(0, 1, 40)]
            ax.plot([p.real for p in pts], [p.imag for p in pts], color='C2', lw=1.4)
            mid = pts[20]
            t = mid * (1 + 0.06)
            ax.text(t.real, t.imag, str(lab[i]), ha='center', va='center', fontsize=8, color='C2')
        ax.plot([0], [0], 'o', color='k', ms=2.5)
        r0 = min(1.0, ext * 1.12)
        ax.set_xlim(-r0, r0); ax.set_ylim(-r0, r0)
        ax.set_title(title, fontsize=9)
        out.append((label, dom['checks']['sides'], dom['checks']['cycles'], [(i, dom['partner'][i]) for i in range(m)]))
    fig.savefig(OUT + 'fig_classical_curves.pdf', bbox_inches='tight')
    plt.close(fig)
    return out


def fig_X3_packings():
    """Remark 6.6: the four extremal 4-packings of X_3 (fuchs_N9#21) drawn in one fixed Dirichlet domain (centred
    at a disc of the packing of the uniformization). Discs shared with another packing are drawn darker."""
    import draw_surface as ds
    from dirichlet import klein_to_poincare
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    G, S, label = ds.load('fuchs_N9', 21, False)
    S.packings.sort(key=lambda P: G.reduce(P['alpha']) is None)       # the packing of the uniformization first
    doms, chi = ds.dirichlet_domains(G, S, 'O')
    dom = doms[0][1]
    V = dom['verts']; m = len(V)
    bnd = []
    for i in range(m):
        k1, k2 = V[i], V[(i + 1) % m]
        bnd += [complex(klein_to_poincare(k1 + t * (k2 - k1))) for t in np.linspace(0, 1, 60)[:-1]]
    ext = max(abs(z) for z in bnd)
    s = float(mp.tanh(G.F.R / 2))
    cents = G.centres(2 * float(mp.atanh(ext)) + 4 * float(G.F.H))[0]
    # the centres of each packing that can meet the domain
    pk = []
    for P in S.packings:
        A = P['alpha']
        lst = []
        for d, z, w in cents:
            wd = disc(act(A, z))
            if abs(wd) < 0.995 and float(2 * mp.atanh(abs(wd))) < 2 * float(mp.atanh(ext)) + 2 * float(G.F.H):
                lst.append((wd, mmul(A, G.word_num(w))))
        pk.append(lst)
    near = lambda z, L: any(abs(z - w) < 1e-7 for w, _ in L)
    shared = [[any(near(wd, pk[j]) for j in range(len(pk)) if j != i) for wd, _ in pk[i]] for i in range(len(pk))]
    cols = [('#dbe9f6', '#9cbde0'), ('#f8dcdc', '#e3a1a1'), ('#dcefd9', '#9fd09a'), ('#efe2f5', '#c9a6dc')]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 7.6))
    for idx, ax in enumerate(axes.flat):
        ax.set_aspect('equal'); ax.axis('off')
        clip = PathPatch(Path([(z.real, z.imag) for z in bnd] + [(bnd[0].real, bnd[0].imag)], closed=True),
                         transform=ax.transData, fc='none', ec='none')
        ax.add_patch(clip)
        light, dark = cols[idx]
        for (wd, M), sh in zip(pk[idx], shared[idx]):
            a = abs(wd) ** 2
            c = wd * (1 - s * s) / (1 - a * s * s)
            r = s * (1 - a) / (1 - a * s * s)
            circ = Circle((c.real, c.imag), r, fc=dark if sh else light, ec='0.5', lw=0.4)
            ax.add_patch(circ); circ.set_clip_path(clip)
            vs = ds.polygon_disc_vertices(G, M)
            vs.sort(key=lambda v: np.angle((v - wd) / (1 - np.conj(wd) * v)))
            for i in range(len(vs)):
                pts = geodesic_pts(vs[i], vs[(i + 1) % len(vs)])
                ln, = ax.plot(pts.real, pts.imag, color='0.45', lw=0.5)
                ln.set_clip_path(clip)
        ax.plot([z.real for z in bnd + bnd[:1]], [z.imag for z in bnd + bnd[:1]], color='C2', lw=1.2)
        ax.plot([0], [0], 'o', color='k', ms=2)
        r0 = ext * 1.06
        ax.set_xlim(-r0, r0); ax.set_ylim(-r0, r0)
        ax.set_title('packing %d' % (idx + 1), fontsize=9, pad=10)
    fig.savefig(OUT + 'fig_X3_packings.pdf', bbox_inches='tight')
    plt.close(fig)
    return [sum(x) for x in shared], [len(x) for x in pk]


if __name__ == '__main__':
    pairs = fig_X0_domain(); print('domain done', pairs)
    print('ED-cycles', fig_edcycles())
    fig_two_tessellations(); print('tessellations done')
    print('Dirichlet domains of X0', fig_dirichlet_X0())
    print('Bolza and y^2=x^5-1', fig_classical_curves())
    print('X3', fig_X3_packings())
