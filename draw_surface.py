"""Picture of one surface: fundamental domain, side pairings and the centres of all its extremal packings.

Usage (from this folder):
    python3 draw_surface.py fuchs_N9 21            # surface fuchs_N9#21 (the surface X_3 of the paper)
    python3 draw_surface.py nec_N9 2               # a non-orientable surface
    python3 draw_surface.py nec_N9 2 --cover       # its canonical double cover
    python3 draw_surface.py fuchs_N9 37 --tess     # also the polygons of the other packings
    python3 draw_surface.py fuchs_N9 21 --dirichlet       # Dirichlet domain centred at O (--dirichlet-all: one
                                                          # centred at a disc of each packing)
    python3 draw_surface.py fuchs_N9 21 --edcycles        # Step 3 of the search (see the end of this text)
Output: surface_figures/<family>_<n>[_cover].pdf (and .png), with the coordinates of the points drawn in a .txt
file of the same name; with --dirichlet, <...>_dirichlet.pdf/.txt; with --edcycles, <...>_edcycles.pdf/.txt.
Takes a few seconds.

What is drawn (Poincare disc, same position as in the paper: polygon P_0 centred at O, an edge midpoint on the
negative imaginary axis):
  * the k polygons of the fundamental domain F built in Step 1 (thin grey lines: edges along which two
    polygons of F are glued; thick black lines: edges of the boundary of F);
  * the side pairings: the two edges identified by a side pairing carry the same number; a number with a
    tilde (e.g. 3~) marks an orientation reversing pairing (glide reflection);
  * the centres of the discs of every extremal packing, one colour per packing, at the points of F where
    they lie. Packing 1 is the one given by the uniformization (centres = centres of the polygons of F).
    A centre lying on the boundary of F appears at each of its copies on the boundary;
  * (with --tess) for the other packings, the regular N-gons of their tessellation centred at those
    points, dashed, in the colour of the packing.
The centres are those of the packings found and verified exactly by the search (surface.py); the picture
itself is drawn with floating point numbers.

With --edcycles the picture shows Step 3 of the search (Section 4 of the paper, Figure 2 for X_6): the central
polygon, its circumscribed disc (dashed) and the search disc (dotted, 0.1 larger), the ED-cycles
b_s(d) = {z : d(z, s z) = d} of the two side pairings s_1 (green) and s_2 (orange) for the admissible distances d
used, all their intersection points in the search disc (small grey dots: discarded by the displacement test) and
the candidates, numbered #n: with the marker of packing n if the candidate is the centre of a disc of that
packing (as in the main picture), with a red cross if it is shown not to be a centre. The .txt file contains the
report of `check_surface.py --details`: every admissible distance, every point #n with its coordinates, the
ED-cycles it lies on and the reason of its fate. When many ED-cycles cross the search disc (typically for N = 7
and 8) the picture can be crowded, with labels on top of each other; the .txt file always gives the complete
information.
"""
import json
import os
import sys
import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from geom import Geometry, to_disc, hdist, act, mmul, minv, mdet
from surface import Surface
from batch import double_cover_action
from wordfmt import compress

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'surface_figures')
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
COLOURS = ['k', 'C3', 'C0', 'C2', 'C1', 'C4', 'C5', 'C6']
MARKERS = ['o', 'D', 's', '^', 'v', 'p', 'h', '*']


def cplx(z):
    return complex(float(z.real), float(z.imag))


def disc(z):
    return cplx(to_disc(z))


def geodesic_pts(z1, z2, n=40):
    z1, z2 = complex(z1), complex(z2)
    w = (z2 - z1) / (1 - np.conj(z1) * z2)
    pts = [(t * w + z1) / (1 + np.conj(z1) * t * w) for t in np.linspace(0, 1, n)]
    return np.array(pts)


def polygon_disc_vertices(G, M):
    """vertices (disc coordinates) of M(P_0), in cyclic order"""
    vs = [disc(act(mmul(M, G.word_num('ca' * j)), G.vA)) for j in range(G.N)]
    c = disc(act(M, mp.mpc(0, 1)))
    return sorted(vs, key=lambda v: np.angle((v - c) / (1 - np.conj(c) * v)))


def draw_polygon(ax, G, M, **kw):
    vs = polygon_disc_vertices(G, M)
    for i in range(len(vs)):
        p = geodesic_pts(vs[i], vs[(i + 1) % len(vs)])
        ax.plot(p.real, p.imag, **kw)


def edge_segment(G, M, j):
    """end points (disc) of edge j of the polygon M(P_0): the edge through M (ca)^j (C)"""
    R = mmul(M, G.word_num('ca' * j))
    # the edge of P_0 through C has end points A and a(A), a being the reflection in the line OC
    v1 = act(R, G.vA)
    v2 = act(mmul(R, G.word_num('a')), G.vA)
    return disc(v1), disc(v2)


def load(family, n, cover):
    groups = [json.loads(l) for l in open(os.path.join(HERE, 'groups', family + '.jsonl')) if l.strip()]
    g = groups[n - 1]
    action = {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']], 'c': [x - 1 for x in g['c']]}
    label = '%s#%d' % (family, n)
    if cover:
        if not family.startswith('nec'):
            sys.exit('--cover only makes sense for the non-orientable families nec_N*')
        action = double_cover_action(action)
        label += ' (canonical double cover)'
    G = Geometry(g['N'])
    S = Surface(G, action, label=label)
    S.analyse(verbose=False)
    return G, S, label


def in_F(G, S, Ms, p, tol=mp.mpf(10) ** -12):
    """is the point p (upper half plane) in the fundamental domain? (Voronoi test in each polygon of F)"""
    for M in Ms:
        q = act(minv(M), p)                 # bring to P_0
        d0 = hdist(q, mp.mpc(0, 1))
        if all(d0 <= hdist(q, act(R, mp.mpc(0, 1))) + tol for R in G.Rota):
            return True
    return False


def packing_centres_in_F(G, S, Ms, alpha):
    """the centres alpha g(O) of the packing that lie in F (including copies on its boundary)"""
    found = []
    rc = G.F.H
    for M in Ms:
        x = act(minv(alpha), act(M, mp.mpc(0, 1)))      # centre of a polygon of F, seen from the packing frame
        # word v with x in v(T_0): reduce the point x into T_0
        q, word = x, []
        while True:
            vio = G.wall_violations(q)
            if not vio:
                break
            q = act(G.gens[vio[0]], q)
            word.append(vio[0])
        V = G.word_num(''.join(word))
        # centres of the tessellation near x: the polygon containing x, its neighbours and theirs
        near = [V] + [mmul(V, R) for R in G.Rota]
        near += [mmul(W, R) for W in near[1:] for R in G.Rota]
        for W in near:
            c = act(alpha, act(W, mp.mpc(0, 1)))
            if hdist(c, act(M, mp.mpc(0, 1))) <= rc + mp.mpf('1e-9') and in_F(G, S, Ms, c):
                if not any(hdist(c, f) < mp.mpf('1e-9') for f, _ in found):
                    found.append((c, mmul(alpha, W)))
    return found


def draw_domain(ax, G, S, polyM):
    """the k polygons of F: edges glued inside F thin grey, boundary edges thick black with the pairing labels"""
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=0.6, color='k'))
    ax.set_aspect('equal'); ax.axis('off')
    for M in polyM.values():
        draw_polygon(ax, G, M, color='0.6', lw=0.6)
    pairs, seen, num = {}, {}, 0
    for (P, kk, Q, m, e, w) in S.sp_raw:
        key = frozenset([(P, kk), (Q, m)])
        if key in seen:
            continue
        num += 1
        seen[key] = num
        for (Pp, j) in ((P, kk), (Q, m)):
            pairs[(Pp, j)] = (num, e)
    for (P, j), (lab, e) in pairs.items():
        z1, z2 = edge_segment(G, polyM[P], j)
        p = geodesic_pts(z1, z2)
        ax.plot(p.real, p.imag, color='k', lw=1.3)
        mid = disc(act(mmul(polyM[P], G.word_num('ca' * j)), G.vC))
        cen = disc(act(polyM[P], mp.mpc(0, 1)))
        t = mid + 0.2 * (cen - mid)
        fs = 7 if abs(mid) < 0.8 else (5.5 if abs(mid) < 0.9 else 4.5)
        ax.text(t.real, t.imag, '%d%s' % (lab, '~' if e else ''), ha='center', va='center', fontsize=fs,
                color='0.15')
    info = []
    done = set()
    for (P, kk, Q, m, e, w) in S.sp_raw:
        key = frozenset([(P, kk), (Q, m)])
        if key in done:
            continue
        done.add(key)
        info.append((seen[key], e, (P, kk), (Q, m), w))
    return sorted(info)


def label_polygons(ax, G, polyM):
    """small grey labels P1, ..., Pk next to the centres of the polygons of F"""
    for i, (P, M) in enumerate(polyM.items()):
        c = disc(act(M, mp.mpc(0, 1)))
        off = 0.045 * (1 - abs(c) ** 2)
        ax.text(c.real + off, c.imag + off, 'P%d' % (i + 1), fontsize=6 if abs(c) < 0.85 else 4.5, color='0.45',
                style='italic')


def fmt_point(w):
    """a point of the Poincare disc: coordinates, hyperbolic distance from O and argument"""
    w = complex(w)
    x, y = (0.0 if abs(w.real) < 5e-7 else w.real), (0.0 if abs(w.imag) < 5e-7 else w.imag)
    d = 2 * np.arctanh(min(abs(w), 1 - 1e-16))
    if d < 1e-9:
        return '%+.6f %+.6fi   (this is O)' % (x, y)
    return '%+.6f %+.6fi   (distance %.6f from O, angle %+8.3f deg)' % (x, y, d, np.degrees(np.angle(w)))


def write_points_file(path, header, sections):
    with open(path, 'w') as f:
        f.write(header + '\n')
        f.write('All the points are given approximately (6 decimals) in the Poincare disc model of the picture: O = 0, the\n'
                'first edge midpoint of the central polygon on the negative imaginary axis; "distance" is the hyperbolic\n'
                'distance from O and "angle" the argument in degrees.\n')
        for title, lines in sections:
            f.write('\n[%s]\n' % title)
            for l in lines:
                f.write(l + '\n')


ROT_MARK = {2: ('o', 4), 3: ('^', 6), 4: ('s', 5.5), 6: ('h', 6.5), 8: ('8', 6.5), 12: ('p', 7)}


def draw_isometries(ax, S, res, dom):
    from isometries import todisc, fromdisc
    handles = []
    # mirrors: the part inside F of the axis of every reflection of the cosets K gamma
    drawn_mirror = False
    for e in res['elements']:
        for (p, q) in e['mirrors']:
            pd, qd = todisc(p), todisc(q)
            u = (qd - pd) / (1 - np.conj(pd) * qd)
            u = u / abs(u)
            ts = np.tanh(np.linspace(-6, 6, 500))
            run = []
            for t in ts:
                w = (t * u + pd) / (1 + np.conj(pd) * t * u)
                z = fromdisc(w)
                if dom.contains(z, tol=1e-7):
                    run.append(w)
                else:
                    if len(run) > 1:
                        ax.plot([x.real for x in run], [x.imag for x in run], color='C0', lw=1.8, alpha=0.75)
                        drawn_mirror = True
                    run = []
            if len(run) > 1:
                ax.plot([x.real for x in run], [x.imag for x in run], color='C0', lw=1.8, alpha=0.75)
                drawn_mirror = True
    if drawn_mirror:
        handles.append(ax.plot([], [], color='C0', lw=1.8, label='mirror (fixed curve of a reflection)')[0])
    orders = sorted(set(m for _, _, m in res['points']))
    for z, allc, m in res['points']:
        mk, ms = ROT_MARK.get(m, ('*', 7))
        for w in allc:
            ax.plot([w.real], [w.imag], mk, color='C4', ms=ms, zorder=6)
    # genus 2, orientable: the hyperelliptic involution is the only one with 6 fixed points
    hyp = [e for e in res['elements'] if S.orientable and e['parity'] == 0 and e['order'] == 2
           and len(e['isolated']) == 6]
    if len(hyp) == 1:
        for key, (z, allc) in hyp[0]['isolated'].items():
            for w in allc:
                ax.plot([w.real], [w.imag], 'o', mfc='none', mec='C3', ms=9, mew=1.2, zorder=7)
        handles.append(ax.plot([], [], 'o', mfc='none', mec='C3', ms=9, mew=1.2, ls='none',
                               label='Weierstrass point (fixed by the hyperelliptic involution)')[0])
    for m in orders:
        mk, ms = ROT_MARK.get(m, ('*', 7))
        handles.append(ax.plot([], [], mk, color='C4', ms=ms, ls='none',
                               label='fixed point of a rotation, stabiliser of order %d' % m)[0])
    return handles


def summary_lines(S, res):
    n = res['order']
    ind = 1 + sum(1 for e in res['elements'] if e['source'] == 'induced')
    nofix = sum(1 for e in res['elements'] if not e['isolated'] and not e['mirrors'])
    withm = sum(1 for e in res['elements'] if e['mirrors'])
    out = []
    if S.orientable:
        conf = res['n_conformal']
        out.append('isometry group of order %d: %d conformal automorphisms and %d anticonformal ones; %d induced by the triangle group'
                   % (n, conf, n - conf, ind))
        out.append('%d anticonformal involutions with mirrors; %d nontrivial isometries without fixed points'
                   % (withm, nofix))
        out.append('check: chi(X/Aut+) = %s, chi(X/Isom) = %s (must be integers; the first one even)'
                   % (fmtnum(res['chi_quotient_conformal']), fmtnum(res['chi_quotient'])))
    else:
        out.append('isometry group of order %d; %d induced by the triangle group' % (n, ind))
        out.append('%d isometries with mirrors; %d nontrivial isometries without fixed points' % (withm, nofix))
        out.append('check: chi(X/Isom) = %s (must be an integer)' % fmtnum(res['chi_quotient']))
    return out


def fmtnum(x):
    return str(int(round(x))) if abs(x - round(x)) < 1e-9 else '%.4f' % x


def make_figure(family, n, cover=False, tess=False, iso=True, outdir=OUTDIR, png=True, quiet=False):
    G, S, label = load(family, n, cover)
    S.packings.sort(key=lambda P: G.reduce(P['alpha']) is None)   # the packing of the uniformization first
    N, k = G.N, S.k
    polyM = {P: G.word_num(S.gword[P]) for P in sorted(S.gword)}
    Ms = list(polyM.values())
    extent = max(abs(v) for M in Ms for v in polygon_disc_vertices(G, M))
    r = min(1.0, extent + 0.06)
    if iso:
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 6.4))
    else:
        fig, ax = plt.subplots(figsize=(6.2, 6.8))
    # left: packings
    pair_info = draw_domain(ax, G, S, polyM)
    label_polygons(ax, G, polyM)
    pnum = {P: i + 1 for i, P in enumerate(polyM)}
    sec_dom = []
    for P, M in polyM.items():
        sec_dom.append('polygon P%d: centre %s' % (pnum[P], fmt_point(disc(act(M, mp.mpc(0, 1))))))
        for j in range(N):
            sec_dom.append('    vertex %2d: %s' % (j, fmt_point(edge_segment(G, M, j)[0])))
    sec_sp = ['%d%s: edge %d of P%d  <->  edge %d of P%d   (%s; the side pairing is the word %s)'
              % (num, '~' if e else '', kk, pnum[P], m, pnum[Q], 'orientation reversing' if e else
                 'orientation preserving', compress(w)) for num, e, (P, kk), (Q, m), w in pair_info]
    sec_sp.insert(0, '(edge j and vertex j of a polygon g(P_0) are the images under g of the edge of P_0 with midpoint '
                     'at angle -90 + 360 j/N degrees and of its end point at angle -90 + 360 (j + 1/2)/N degrees)')
    sec_pk = []
    all_groups = []
    from isometries import Domain
    fdom = Domain(S)
    lines, handles = [], []
    for idx, Pk in enumerate(S.packings):
        col, mk = COLOURS[idx % len(COLOURS)], MARKERS[idx % len(MARKERS)]
        pts = packing_centres_in_F(G, S, Ms, Pk['alpha'])
        for c, W in pts:
            z = disc(c)
            ax.plot([z.real], [z.imag], mk, color=col, ms=6 if idx else 3.5, zorder=5 if idx else 7)
            if tess and idx > 0:
                draw_polygon(ax, G, W, color=col, lw=0.8, ls=(0, (4, 2)), alpha=0.8)
        # group the points drawn into discs of X (copies on the boundary of F are the same disc)
        groups = []
        for c, W in pts:
            zc = complex(c)
            for g in groups:
                if any(abs(zc - u) < 1e-7 * (1 + abs(u)) for u in g['copies']):
                    g['drawn'].append(disc(c))
                    break
            else:
                cps = fdom.copies(zc)
                groups.append({'copies': cps, 'drawn': [disc(c)]})
        all_groups.append(groups)
        sec_pk.append('packing %d%s: %d discs' % (idx + 1, ' (uniformization)' if idx == 0 else '', len(groups)))
        for gi, g in enumerate(groups):
            sec_pk.append('    disc %d: centre %s' % (gi + 1, fmt_point(g['drawn'][0])))
            for w2 in g['drawn'][1:]:
                sec_pk.append('            copy  %s' % fmt_point(w2))
        dz = hdist(mp.mpc(0, 1), Pk['z'])
        lines.append('packing %d: %d point%s drawn (one disc centred at distance %s from O)'
                     % (idx + 1, len(pts), '' if len(pts) == 1 else 's', mp.nstr(dz, 6)))
        handles.append(ax.plot([], [], mk, color=col, ms=6 if idx else 3.5, ls='none',
                               label='packing %d%s' % (idx + 1, ' (uniformization)' if idx == 0 else ''))[0])
    ax.legend(handles=handles, loc='upper left', fontsize=7, frameon=False, bbox_to_anchor=(-0.12, 1.0))
    ax.set_xlim(-r, r); ax.set_ylim(-r, r)
    ax.set_title('extremal %d-packings' % k, fontsize=9)
    kind = 'orientable' if S.orientable else 'non-orientable'
    title = '%s:  N = %d, k = %d, %s, %d extremal %d-packing%s' % (
        label, N, k, kind, len(S.packings), k, '' if len(S.packings) == 1 else 's')
    summary = {'surface': label, 'N': N, 'k': k, 'orientable': S.orientable, 'n_packings': len(S.packings)}
    shared = []
    for i1 in range(len(all_groups)):
        for i2 in range(i1 + 1, len(all_groups)):
            for gi, g in enumerate(all_groups[i1]):
                for hj, h in enumerate(all_groups[i2]):
                    if any(abs(u - v) < 1e-7 * (1 + abs(v)) for u in g['copies'] for v in h['copies']):
                        shared.append('disc %d of packing %d = disc %d of packing %d' % (gi + 1, i1 + 1, hj + 1, i2 + 1))
    sec_pk.append('discs shared by two different packings: %s' % ('; '.join(shared) if shared else 'none'))
    summary['shared_discs'] = len(shared)
    sections = [('fundamental domain of Step 1: polygons P1, ..., P%d' % k, sec_dom), ('side pairings', sec_sp),
                ('centres of the discs of the extremal packings (a centre on the boundary of the domain is listed with '
                 'its other copies)', sec_pk)]
    # right: isometries
    if iso:
        from isometries import analyse_isometries
        res, dom = analyse_isometries(S)
        draw_domain(ax2, G, S, polyM)
        h2 = draw_isometries(ax2, S, res, dom)
        sec_iso = summary_lines(S, res)[:2]
        for e in res['elements']:
            kind = ('conformal' if e['parity'] == 0 else 'anticonformal') if S.orientable else 'isometry'
            fx = len(e['isolated'])
            sec_iso.append('%s of order %d (%s): %d isolated fixed point%s%s' % (
                kind, e['order'], 'induced by the triangle group' if e['source'] == 'induced' else 'hidden', fx,
                '' if fx == 1 else 's', ', fixed curves (mirrors)' if e['mirrors'] else ''))
        sec_fix = []
        for z, allc, m in sorted(res['points'], key=lambda t: (-t[2], abs(t[0]))):
            sec_fix.append('stabiliser of order %d: %s%s' % (m, fmt_point(z), '' if len(allc) == 1 else
                                                              '   (%d copies on the boundary of the domain)' % len(allc)))
        hyp = [e for e in res['elements'] if S.orientable and e['parity'] == 0 and e['order'] == 2
               and len(e['isolated']) == 6]
        sec_w = [fmt_point(z) for z, allc in (hyp[0]['isolated'].values() if len(hyp) == 1 else [])]
        sections += [('isometries (identity excluded)', sec_iso),
                     ('points fixed by isometries acting as rotations', sec_fix)]
        if sec_w:
            sections.append(('Weierstrass points (fixed points of the hyperelliptic involution)', sec_w))
        ax2.legend(handles=h2, loc='upper left', fontsize=7, frameon=False, bbox_to_anchor=(-0.12, 1.0))
        ax2.set_xlim(-r, r); ax2.set_ylim(-r, r)
        ax2.set_title('isometries: fixed points and mirrors', fontsize=9)
        lines += summary_lines(S, res)
        summary.update({'isometries': res['order'], 'induced': 1 + sum(1 for e in res['elements'] if e['source'] == 'induced'),
                        'chi_quotient': res['chi_quotient'], 'conformal': res.get('n_conformal'),
                        'chi_quotient_conformal': res.get('chi_quotient_conformal')})
    fig.suptitle(title, fontsize=10)
    lines.append('edges with the same number are identified; n~ = orientation reversing pairing; '
                 'points on the boundary of the domain appear at each of their copies')
    fig.text(0.5, 0.07, '\n'.join(lines), ha='center', va='top', fontsize=7)
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, '%s_%d%s' % (family, n, '_cover' if cover else ''))
    fig.savefig(base + '.pdf', bbox_inches='tight')
    if png:
        fig.savefig(base + '.png', bbox_inches='tight', dpi=130)
    plt.close(fig)
    write_points_file(base + '.txt', title, sections)
    if not quiet:
        print(title)
        for l in lines[:-1]:
            print('   ' + l)
        print('written:', base + '.pdf' + (' and .png' if png else '') + ', coordinates in ' + base + '.txt')
    return summary


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        print(__doc__)
        sys.exit(1)
    if '--edcycles' in sys.argv:
        make_edcycles_figure(args[0], int(args[1]), cover='--cover' in sys.argv)
    elif '--dirichlet' in sys.argv or '--dirichlet-all' in sys.argv:
        make_dirichlet_figure(args[0], int(args[1]), which='all' if '--dirichlet-all' in sys.argv else 'O',
                              cover='--cover' in sys.argv)
    else:
        make_figure(args[0], int(args[1]), cover='--cover' in sys.argv, tess='--tess' in sys.argv,
                    iso='--no-isometries' not in sys.argv)



# ---------------------------------------------------------------- Dirichlet domains
def dirichlet_domains(G, S, which='O'):
    """list of (packing index, domain in the frame of S); which = 'O' (centre O only) or 'all' (one per packing)"""
    from dirichlet import dirichlet_at_O, transport
    from gapiface import coset_action_from_words, reduce_word
    chi = S.d // (2 * G.N) + S.d // 6 + S.d // 4 - 3 * S.d // 2 + S.d
    out = [(0, dict(dirichlet_at_O(G, S.perm, S.d, chi), centre=mp.mpc(0)))]
    if which == 'all':
        for idx, P in enumerate(S.packings[1:], 1):
            a2 = coset_action_from_words(G.N, [reduce_word(w) for w in P['words']])
            perm2 = {'a': a2['a'], 'b': a2['b'], 'c': a2['c']}
            dom = dirichlet_at_O(G, perm2, a2['d'], chi)
            out.append((idx, transport(dom, P['alpha'])))
    return out, chi


def draw_dirichlet_panel(ax, G, S, polyM, idx, dom):
    """the Dirichlet domain, recentred so that its centre is the centre of the picture"""
    from dirichlet import klein_to_poincare
    ctr = complex(klein_to_poincare(dom['centre']))
    T = lambda w: (w - ctr) / (1 - np.conj(ctr) * w)          # isometry of the disc taking the centre to 0
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=0.6, color='k'))
    ax.set_aspect('equal'); ax.axis('off')
    for M in polyM.values():                                     # the domain of Step 1, faint
        vs = [T(v) for v in polygon_disc_vertices(G, M)]
        c = T(disc(act(M, mp.mpc(0, 1))))
        vs.sort(key=lambda v: np.angle((v - c) / (1 - np.conj(c) * v)))
        for i in range(len(vs)):
            p = geodesic_pts(vs[i], vs[(i + 1) % len(vs)])
            ax.plot(p.real, p.imag, color='0.82', lw=0.5)
    V = dom['verts']
    m = len(V)
    lab, num = {}, 0
    for i in range(m):
        j = dom['partner'][i]
        if i in lab:
            continue
        num += 1
        lab[i] = num
        if j is not None:
            lab[j] = num
    ext = 0
    for i in range(m):
        k1, k2 = V[i], V[(i + 1) % m]
        pts = [T(complex(klein_to_poincare(k1 + t * (k2 - k1)))) for t in np.linspace(0, 1, 30)]
        ax.plot([p.real for p in pts], [p.imag for p in pts], color='C2', lw=1.4)
        ext = max(ext, max(abs(p) for p in pts))
        mid = pts[15]
        t = mid * (1 - 0.13 * (1 - abs(mid) ** 2))
        rev = mdet(dom['mats'][i]) < 0
        ax.text(t.real, t.imag, '%d%s' % (lab[i], '~' if rev else ''), ha='center', va='center', fontsize=6.5,
                color='C2')
    col, mk = COLOURS[idx % len(COLOURS)], MARKERS[idx % len(MARKERS)]
    ax.plot([0], [0], mk, color=col, ms=6 if idx else 4, zorder=6)
    c = dom['checks']
    ok = 'all checks passed' if c['ok'] else 'CHECK FAILED'
    ax.set_title('centred at a disc of packing %d: %d sides, vertex cycles %s\n%s' % (
        idx + 1, c['sides'], '+'.join(str(x) for x in c['cycles']), ok), fontsize=8)
    return ext


def make_dirichlet_figure(family, n, which='O', cover=False, outdir=OUTDIR, png=True, quiet=False):
    G, S, label = load(family, n, cover)
    S.packings.sort(key=lambda P: G.reduce(P['alpha']) is None)
    polyM = {P: G.word_num(S.gword[P]) for P in sorted(S.gword)}
    doms, chi = dirichlet_domains(G, S, which)
    nd = len(doms)
    cols = min(nd, 3)
    rows = (nd + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.2 * cols, 5.4 * rows + 0.6), squeeze=False)
    exts = []
    for (idx, dom), ax in zip(doms, axes.flat):
        exts.append(draw_dirichlet_panel(ax, G, S, polyM, idx, dom))
    for ax in list(axes.flat)[nd:]:
        ax.axis('off')
    r = min(1.0, max(exts) + 0.06)
    for ax in list(axes.flat)[:nd]:
        ax.set_xlim(-r, r); ax.set_ylim(-r, r)
    kind = 'orientable' if S.orientable else 'non-orientable'
    fig.suptitle('%s: N = %d, k = %d, %s, chi = %d -- Dirichlet domains (green; sides with the same number are '
                 'identified, n~ = orientation reversing)\neach picture is centred at the centre of its domain; the polygons of the fundamental domain of Step 1 in grey'
                 % (label, G.N, S.k, kind, chi), fontsize=9)
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, '%s_%d%s_dirichlet' % (family, n, '_cover' if cover else ''))
    fig.savefig(base + '.pdf', bbox_inches='tight')
    if png:
        fig.savefig(base + '.png', bbox_inches='tight', dpi=130)
    plt.close(fig)
    from dirichlet import klein_to_poincare
    sections = []
    for idx, d in doms:
        c = complex(klein_to_poincare(d['centre']))
        V = [complex(klein_to_poincare(v)) for v in d['verts']]
        m = len(V)
        lab, num = {}, 0
        for i in range(m):
            if i not in lab:
                num += 1
                lab[i] = num
                if d['partner'][i] is not None:
                    lab[d['partner'][i]] = num
        hd = lambda u, v: 2 * np.arctanh(abs((u - v) / (1 - np.conj(v) * u)))
        ch = d['checks']
        sec = ['centre: %s' % fmt_point(c),
               '%d sides, vertex cycles %s, largest distance from the centre to a vertex %.6f; checks: all sides paired %s, '
               'vertex cycles with angle sum 2 pi %s, area 2 pi |chi| %s' % (
                   m, '+'.join(map(str, ch['cycles'])), float(ch['rho_max']), ch['all_paired'], ch['cycle_sums_ok'],
                   ch['area_ok'])]
        for i in range(m):
            sec.append('vertex V%d: %s; distance %.6f from the centre' % (i + 1, fmt_point(V[i]), hd(V[i], c)))
        for i in range(m):
            j = d['partner'][i]
            rev = mdet(d['mats'][i]) < 0
            sec.append('side %d%s: from V%d to V%d, paired with the side from V%d to V%d by kappa^-1, kappa = %s' % (
                lab[i], '~' if rev else '', i + 1, (i + 1) % m + 1, (j if j is not None else -1) + 1,
                (j + 1) % m + 1 if j is not None else 0, compress(d['words'][i])))
        sections.append(('Dirichlet domain centred at a disc of packing %d' % (idx + 1), sec))
    write_points_file(base + '.txt', 'Dirichlet domains of %s (each picture is recentred at the centre of its domain; '
                      'the coordinates below are in the frame of the surface, as in %s_%d.txt)' % (label, family, n),
                      sections)
    summary = {'surface': label, 'dirichlet': [(idx, d['checks']) for idx, d in doms]}
    if not quiet:
        for idx, d in doms:
            c = d['checks']
            print('packing %d: Dirichlet domain with %d sides, vertex cycles %s, area = 2 pi |chi|: %s, '
                  'cycle sums 2 pi: %s, all sides paired: %s' % (idx + 1, c['sides'], c['cycles'], c['area_ok'],
                                                                c['cycle_sums_ok'], c['all_paired']))
        print('written:', base + '.pdf' + (' and .png' if png else ''))
    return summary


# ---------------------------------------------------------------- the ED-cycles of Step 3
def ed_cycle_curves(G, s, d, view):
    """the hypercycles of b_s(d) as lists of points of the disc (only the part with |w| <= view)"""
    curves = []
    for (A, B, C, D) in s.ed_cycle(d):
        A, B, C, D = (float(x) for x in (A, B, C, D))
        if abs(A) > 1e-14 * max(abs(B), abs(C), 1):           # circle x^2 + y^2 + (B/A) x + (C/A) y + D/A = 0
            cx, cy = -B / (2 * A), -C / (2 * A)
            r = np.sqrt(max(cx * cx + cy * cy - D / A, 0))
            t = np.linspace(0, 2 * np.pi, 4000)
            z = (cx + r * np.cos(t)) + 1j * (cy + r * np.sin(t))
        else:                                                 # line B x + C y + D = 0
            t = np.concatenate([-np.logspace(4, -4, 2000), np.logspace(-4, 4, 2000)])
            if abs(C) > abs(B):
                z = t + 1j * (-(D + B * t) / C)
            else:
                y = np.abs(t)
                z = (-(D + C * y) / B) + 1j * y
        z = z[z.imag > 0]
        w = (z - 1j) / (z + 1j) * 1j                         # the disc model of the paper (as to_disc)
        w = np.where(np.abs(w) <= view, w, np.nan)
        curves.append(w)
    return curves


def make_edcycles_figure(family, n, cover=False, outdir=OUTDIR, png=True, quiet=False):
    """Step 3 of the search for one surface: the ED-cycles of the two side pairings s_1, s_2 in the search disc,
    all their intersection points, the candidates and what happens to them. Also writes the report of
    check_surface.py --details to <...>_edcycles.txt."""
    from check_surface import details_lines
    G, S, label = load(family, n, cover)
    S.packings.sort(key=lambda P: G.reduce(P['alpha']) is None)   # packing 1 = the one of the uniformization
    se = S.search
    s1, s2 = se['s1'], se['s2']
    rs = float(mp.tanh(se['radius'] / 2))                          # Euclidean radius of the search disc
    rc = float(mp.tanh(G.F.H / 2))                                 # circumscribed disc of the central polygon
    view = rs * 1.12
    fig, ax = plt.subplots(figsize=(7.6, 8.2))
    ax.set_aspect('equal'); ax.axis('off')
    draw_polygon(ax, G, G.word_num(''), color='k', lw=1.0)
    ax.add_patch(Circle((0, 0), rc, fill=False, lw=0.8, ls='--', color='0.3'))
    ax.add_patch(Circle((0, 0), rs, fill=False, lw=0.6, ls=':', color='0.3'))
    for s, Dm, col in ((s1, se['D1'], 'C2'), (s2, se['D2'], 'C1')):
        for d in S.dists:
            if d > Dm:
                continue
            for w in ed_cycle_curves(G, s, d, view):
                ax.plot(w.real, w.imag, color=col, lw=0.55, alpha=0.85)
    pts = se['points']
    for num, p in enumerate(pts, 1):
        w = complex(to_disc(p['z']))
        if p['killed_by'] is not None:
            ax.plot([w.real], [w.imag], 'o', color='0.55', ms=2.5, zorder=4)
            continue
        if p['packings']:
            idx = sorted(S.packings.index(P) for P in p['packings'])
            i0 = idx[0]
            ax.plot([w.real], [w.imag], MARKERS[i0 % len(MARKERS)], color=COLOURS[i0 % len(COLOURS)], ms=7.5,
                    zorder=6, mec='k', mew=0.5)
            txt = '#%d: P%s' % (num, ',P'.join(str(i + 1) for i in idx))
        else:
            ax.plot([w.real], [w.imag], 'x', color='C3', ms=8, mew=1.6, zorder=6)
            txt = '#%d: not a centre' % num
        ax.annotate(txt, (w.real, w.imag), xytext=(5, 5), textcoords='offset points', fontsize=7, zorder=7,
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.8))
    ax.set_xlim(-view, view); ax.set_ylim(-view, view)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color='C2', lw=1.2, label='ED-cycles of s_1 = %s' % compress(s1.word)),
               Line2D([], [], color='C1', lw=1.2, label='ED-cycles of s_2 = %s' % compress(s2.word)),
               Line2D([], [], color='0.3', lw=0.8, ls='--', label='circumscribed disc of the central polygon'),
               Line2D([], [], color='0.3', lw=0.6, ls=':', label='search disc (radius + 0.1)'),
               Line2D([], [], marker='o', color='0.55', lw=0, ms=4, label='intersection discarded by the displacement test'),
               Line2D([], [], marker='o', color='w', mec='k', lw=0, ms=7,
                      label='candidate: centre of a disc of packing Pn (marker of packing n)'),
               Line2D([], [], marker='x', color='C3', lw=0, ms=7, mew=1.6, label='candidate that is not a centre')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.01), fontsize=7.5, frameon=False, ncol=2)
    ncand = sum(1 for p in pts if p['killed_by'] is None)
    ax.set_title('%s: N = %d, k = %d -- Step 3 of the search: %d intersection points, %d candidates, %d packings\n'
                 '(numbers #n as in the file %s_%d%s_edcycles.txt, which lists every point and every ED-cycle)'
                 % (label, G.N, S.k, len(pts), ncand, len(S.packings), family, n, '_cover' if cover else ''),
                 fontsize=8.5)
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, '%s_%d%s_edcycles' % (family, n, '_cover' if cover else ''))
    fig.savefig(base + '.pdf', bbox_inches='tight')
    if png:
        fig.savefig(base + '.png', bbox_inches='tight', dpi=130)
    plt.close(fig)
    with open(base + '.txt', 'w') as f:
        f.write('Step 3 of the search for %s (the picture %s_edcycles.pdf). Points in the Poincare disc model of the\n'
                'paper: O = 0, the first edge midpoint of the central polygon on the negative imaginary axis.\n\n'
                % (label, os.path.basename(base)[:-len('_edcycles')]))
        f.write('\n'.join(details_lines(G, S)) + '\n')
    if not quiet:
        print('%s: %d intersection points, %d candidates, %d packings' % (label, len(pts), ncand, len(S.packings)))
        print('written:', base + '.pdf' + (' and .png' if png else '') + ', ' + base + '.txt')
    decided = all(p['killed_by'] is not None or p['packings'] or p['rejection'] is not None for p in pts)
    return {'surface': label, 'points': len(pts), 'candidates': ncand, 'packings': len(S.packings),
            'all_decided': decided}


if __name__ == '__main__':
    main()
