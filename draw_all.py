"""Draw every surface of the paper and check the pictures against the data files.

Usage:  python3 draw_all.py [number of parallel processes, default 2] [--only-figures | --only-dirichlet |
                                                                       --only-edcycles]
Output: folder surface_figures/all/ with, for each of the 373 conjugacy classes, <family>_<n>.pdf (packings and
        isometries, as in draw_surface.py), <family>_<n>_dirichlet.pdf (Dirichlet domains centred at a disc of
        each packing, as in draw_surface.py --dirichlet-all) and <family>_<n>_edcycles.pdf with
        <family>_<n>_edcycles.txt (the ED-cycles, intersection points and candidates of Step 3, as in
        draw_surface.py --edcycles), the table surface_figures/all/index.csv (one line per surface: packings,
        isometries, Euler characteristic checks) and the report surface_figures/all/checks.txt.
Time: about an hour with two processes.

For each surface the program compares what it computed while drawing with the data of the paper:
  * the number of extremal packings and the order of the group induced by the triangle group, with
    groups/surfaces_summary.csv;
  * for orientable surfaces, the order of the full isometry group with the one implied by
    groups/surface_identifications.json (induced order times the number of packings related to the original one
    by an isometry);
  * the checks of every Dirichlet domain (all sides paired, vertex cycles with angle sum 2 pi, area 2 pi |chi|,
    at most 18 sides for orientable and 12 for non-orientable surfaces);
  * for the ED-cycle pictures, the number of candidates with groups/<family>_results.json (field n_candidates),
    and that every candidate is either the centre of a disc of a packing or shown not to be a centre;
  * the Euler characteristic of the quotient by the isometry group (and by the conformal automorphism group), which
    must be an integer (an even integer at most 2 for the conformal quotient), as a check that all the fixed points
    of all the isometries have been found.
The canonical double covers of the non-orientable surfaces are among the orientable classes (Appendix A of the
paper), so they are not drawn separately; `python3 draw_surface.py nec_N9 2 --cover` draws one of them on the
fundamental domain of the non-orientable surface.
"""
import csv
import json
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'surface_figures', 'all')
FAMILIES = ['nec_N12', 'nec_N9', 'nec_N8', 'nec_N7', 'fuchs_N18', 'fuchs_N12', 'fuchs_N10', 'fuchs_N9',
            'fuchs_N8', 'fuchs_N7']


def job(arg):
    fam, n, parts = arg
    import matplotlib
    matplotlib.use('Agg')
    from draw_surface import make_figure, make_dirichlet_figure, make_edcycles_figure
    t = time.time()
    s = {'surface': '%s#%d' % (fam, n), 'error': ''}
    try:
        if 'figures' in parts:
            s.update(make_figure(fam, n, outdir=OUT, png=False, quiet=True))
        if 'dirichlet' in parts:
            s['dirichlet'] = make_dirichlet_figure(fam, n, which='all', outdir=OUT, png=False, quiet=True)['dirichlet']
        if 'edcycles' in parts:
            s['edcycles'] = make_edcycles_figure(fam, n, outdir=OUT, png=False, quiet=True)
        s['seconds'] = round(time.time() - t, 1)
    except Exception as e:                                    # reported, not hidden
        s['error'] = repr(e)
    return s


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    nproc = int(args[0]) if args else 2
    parts = {'figures', 'dirichlet', 'edcycles'}
    for p in ('figures', 'dirichlet', 'edcycles'):
        if '--only-' + p in sys.argv:
            parts = {p}
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for fam in FAMILIES:
        n = sum(1 for l in open(os.path.join(HERE, 'groups', fam + '.jsonl')) if l.strip())
        jobs += [(fam, i, parts) for i in range(1, n + 1)]
    summary = {r['surface']: r for r in csv.DictReader(open(os.path.join(HERE, 'groups', 'surfaces_summary.csv')))}
    ident = {}
    for fam, d in json.load(open(os.path.join(HERE, 'groups', 'surface_identifications.json'))).items():
        ident.update(d)
    ncand = {}
    for fam in FAMILIES:
        for i, r in enumerate(json.load(open(os.path.join(HERE, 'groups', fam + '_results.json'))), 1):
            ncand['%s#%d' % (fam, i)] = (r['n_candidates'], r['n_packings'])
    results = []
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap(job, jobs), 1):
            results.append(r)
            print('%3d/%d  %-14s %s' % (i, len(jobs), r['surface'], r.get('error') or '%ss' % r.get('seconds')),
                  flush=True)
    problems = []
    rows = []
    for r in results:
        name = r['surface']
        if r.get('error'):
            problems.append('%s: error while drawing: %s' % (name, r['error']))
            continue
        for idx, c in r.get('dirichlet', []):
            maxsides = 18 if r['surface'].startswith('fuchs') else 12
            if not c['ok']:
                problems.append('%s: Dirichlet domain centred at packing %d fails a check: %s' % (name, idx + 1, c))
            if c['sides'] > maxsides:
                problems.append('%s: Dirichlet domain with %d sides' % (name, c['sides']))
        if 'edcycles' in r:
            e = r['edcycles']
            if (e['candidates'], e['packings']) != ncand[name]:
                problems.append('%s: ED-cycle picture with %d candidates and %d packings, %d and %d in the results file'
                                % (name, e['candidates'], e['packings'], *ncand[name]))
            if not e['all_decided']:
                problems.append('%s: some candidate is neither a centre nor shown not to be a centre' % name)
        if 'figures' not in parts:
            rows.append(r)
            continue
        ref = summary[name]
        if int(ref['n_packings']) != r['n_packings']:
            problems.append('%s: %d packings drawn, %s in surfaces_summary.csv' % (name, r['n_packings'], ref['n_packings']))
        if int(ref['induced_isometries']) != r['induced']:
            problems.append('%s: induced group of order %d, %s in surfaces_summary.csv'
                            % (name, r['induced'], ref['induced_isometries']))
        cq = r['chi_quotient']
        if abs(cq - round(cq)) > 1e-9:
            problems.append('%s: chi(X/Isom) = %.6f is not an integer' % (name, cq))
        if r['orientable']:
            cc = r['chi_quotient_conformal']
            if abs(cc - round(cc)) > 1e-9 or round(cc) % 2 or round(cc) > 2:
                problems.append('%s: chi(X/Aut+) = %.6f is not an even integer <= 2' % (name, cc))
            idd = ident.get(name)
            if idd is not None:
                expected = idd['induced_aut'] * (1 + idd['hidden_conf_links'] + idd['hidden_anti_links'])
                if expected != r['isometries']:
                    problems.append('%s: isometry group of order %d, %d expected from surface_identifications.json'
                                    % (name, r['isometries'], expected))
        rows.append(r)
    for r in rows:
        if 'dirichlet' in r:
            r['dirichlet_sides'] = ' '.join(str(c['sides']) for _, c in r['dirichlet'])
            r['dirichlet_vertex_cycles'] = ' '.join('+'.join(map(str, c['cycles'])) for _, c in r['dirichlet'])
    keys = ['surface', 'N', 'k', 'orientable', 'n_packings', 'induced', 'isometries', 'conformal',
            'chi_quotient', 'chi_quotient_conformal', 'shared_discs', 'dirichlet_sides', 'dirichlet_vertex_cycles',
            'ed_points', 'ed_candidates']
    for r in rows:
        if 'edcycles' in r:
            r['ed_points'], r['ed_candidates'] = r['edcycles']['points'], r['edcycles']['candidates']
    name = 'index.csv' if parts == {'figures', 'dirichlet', 'edcycles'} else 'index_%s.csv' % '_'.join(sorted(parts))
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    name = 'checks.txt' if parts == {'figures', 'dirichlet', 'edcycles'} else 'checks_%s.txt' % '_'.join(sorted(parts))
    with open(os.path.join(OUT, name), 'w') as f:
        f.write('%d surfaces drawn, %d problems\n' % (len(rows), len(problems)))
        for p in problems:
            f.write(p + '\n')
    print('%d surfaces drawn, %d problems (see %s)' % (len(rows), len(problems), os.path.join(OUT, name)))


if __name__ == '__main__':
    main()
