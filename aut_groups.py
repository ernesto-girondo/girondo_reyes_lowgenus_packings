"""The full isometry group of every surface of the paper, and its group of conformal automorphisms in the orientable
case, identified up to isomorphism by GAP.

Usage:  python3 aut_groups.py [number of parallel processes, default 2]
Output: groups/automorphism_groups.txt and groups/automorphism_groups.csv
Time:   about an hour with two processes.

For each of the 373 classes (the 50 canonical double covers are among the orientable classes, Appendix A of the
paper), the isometries of the surface X = D/K are obtained as in autfix.all_automorphisms: those induced by the
triangle group, and, for every hidden packing that is the image of the original one under an isometry, the
isometries that carry the original packing to it. Every isometry of X carries the original packing to an extremal
packing, so this list is complete; its length is compared with the order computed independently in
groups/surface_identifications.json.

The multiplication table is computed on X: an isometry of X is determined by the image of one generic point z0, so
the product of two elements is identified by reducing the image of z0 into the fundamental domain F. The table
defines the group as a group of permutations (its regular representation), and GAP gives its isomorphism type
(StructureDescription) and its number in the Small Groups Library (IdGroup). For orientable surfaces the subgroup
of conformal automorphisms (the elements that preserve the orientation) is identified in the same way.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool
import mpmath as mp
from geom import Geometry, act, mmul, minv
from surface import Surface
from autfix import all_automorphisms, reduce_mod_K
from isometries import reduce_mod_K_matrix

mp.mp.dps = 40
FAMILIES = ['nec_N12', 'nec_N9', 'nec_N8', 'nec_N7', 'fuchs_N18', 'fuchs_N12', 'fuchs_N10', 'fuchs_N9', 'fuchs_N8',
            'fuchs_N7']
Z0 = mp.mpc('0.1234567891', '1.3456789123')          # a generic point of the upper half plane


def point_key(S, z):
    w = reduce_mod_K(S, z)
    return (round(float(w.real), 7), round(float(w.imag), 7))


def cayley_table(S, mats):
    keys = [point_key(S, act(M, Z0)) for M in mats]
    if len(set(keys)) != len(keys):
        raise RuntimeError('two isometries move z0 to the same point')
    index = {k: i for i, k in enumerate(keys)}
    n = len(mats)
    table = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            # (M_i M_j)(z0) = M_i(M_j(z0)); keys[j] is M_j(z0) reduced mod K, and M_i is a lift of an isometry of X
            z = act(mats[i], act(mats[j], Z0))
            k = point_key(S, z)
            if k not in index:
                raise RuntimeError('product not found')
            table[i][j] = index[k]
    return table


def surface_group(fam, n, g):
    G = Geometry(g['N'])
    S = Surface(G, {'d': g['d'], 'a': [x - 1 for x in g['a']], 'b': [x - 1 for x in g['b']],
                    'c': [x - 1 for x in g['c']]}, label='%s#%d' % (fam, n))
    S.analyse(verbose=False)
    auts = all_automorphisms(S)
    mats = [reduce_mod_K_matrix(S, m) for m, _, _ in auts]
    table = cayley_table(S, mats)
    ident = [i for i in range(len(mats)) if all(table[i][j] == j for j in range(len(mats)))]
    if len(ident) != 1:
        raise RuntimeError('no unique identity')
    parity = [p for _, _, p in auts]
    return {'surface': S.label, 'N': g['N'], 'k': S.k if hasattr(S, 'k') else None, 'orientable': S.orientable,
            'n_packings': len(S.packings), 'table': table,
            'conformal': [i for i in range(len(mats)) if parity[i] == 0] if S.orientable else None}


def gap_identify(groups):
    """groups: list of lists of permutations (as lists of images, 0-based). Returns list of (description, id)."""
    lines = ['res := [];']
    for perms in groups:
        gens = ['PermList([%s])' % ','.join(str(x + 1) for x in p) for p in perms]
        lines.append('G := Group(%s);' % (', '.join(gens) if gens else '()'))
        lines.append('Add(res, [StructureDescription(G), IdGroup(G)]);')
    lines.append('for r in res do Print(r[1], "|", r[2], "\\n"); od;')
    lines.append('QUIT;')
    with tempfile.NamedTemporaryFile('w', suffix='.g', delete=False) as f:
        f.write('\n'.join(lines))
        name = f.name
    out = subprocess.run(['gap', '-q', '-o', '4g', name], capture_output=True, text=True,
                         stdin=subprocess.DEVNULL).stdout
    os.unlink(name)
    res = []
    for l in out.splitlines():
        if '|' in l:
            a, b = l.split('|')
            res.append((a.strip(), b.strip().replace(' ', '')))
    if len(res) != len(groups):
        raise RuntimeError('GAP returned %d answers for %d groups' % (len(res), len(groups)))
    return res


def regular_perms(table, elems):
    """the left regular representation of the subgroup formed by `elems` (indices into the table)"""
    pos = {e: i for i, e in enumerate(elems)}
    return [[pos[table[x][y]] for y in elems] for x in elems]


def job(arg):
    fam, n, g = arg
    return surface_group(fam, n, g)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    nproc = int(args[0]) if args else 2
    ident = {}
    for fam, d in json.load(open('groups/surface_identifications.json')).items():
        ident.update(d)
    summary = {r['surface']: r for r in csv.DictReader(open('groups/surfaces_summary.csv'))}
    jobs = []
    for fam in FAMILIES:
        groups = [json.loads(l) for l in open('groups/%s.jsonl' % fam) if l.strip()]
        jobs += [(fam, n, g) for n, g in enumerate(groups, 1)]
    rows = []
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap(job, jobs), 1):
            r['k'] = int(summary[r['surface']]['k'])
            rows.append(r)
            print('%3d/%d  %-14s |Isom| = %d' % (i, len(jobs), r['surface'], len(r['table'])), flush=True)
    todo = []
    for r in rows:
        allel = list(range(len(r['table'])))
        todo.append(regular_perms(r['table'], allel))
        if r['orientable']:
            todo.append(regular_perms(r['table'], r['conformal']))
    answers = iter(gap_identify(todo))
    problems = []
    for r in rows:
        r['isom'] = next(answers)
        r['aut'] = next(answers) if r['orientable'] else None
        idd = ident.get(r['surface'])
        if idd is not None:
            expected = idd['induced_aut'] * (1 + idd['hidden_conf_links'] + idd['hidden_anti_links'])
            if expected != len(r['table']):
                problems.append('%s: %d isometries, %d expected from surface_identifications.json'
                                % (r['surface'], len(r['table']), expected))
    out = ['Full isometry group Isom(X) of every surface of the paper and, for orientable surfaces, the group Aut+(X)',
           'of conformal automorphisms. Isomorphism type: GAP StructureDescription; [order, number] in the Small',
           'Groups Library (IdGroup). Cn = cyclic group of order n, D2n = dihedral group of order 2n (D_n in the paper), x = direct',
           'product, : = semidirect product.', '']
    for r in rows:
        line = '%-14s N=%-2d k=%-2d packings=%d  Isom: %-26s %-9s' % (r['surface'], r['N'], r['k'], r['n_packings'],
                                                                    r['isom'][0], r['isom'][1])
        if r['orientable']:
            line += '  Aut+: %-12s %s' % r['aut']
        out.append(line)
    out.append('')
    out.append('%d surfaces, %d problems' % (len(rows), len(problems)))
    out += problems
    open('groups/automorphism_groups.txt', 'w').write('\n'.join(out) + '\n')
    with open('groups/automorphism_groups.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['surface', 'N', 'k', 'orientable', 'n_packings', 'isometry_group', 'isometry_group_id',
                    'conformal_group', 'conformal_group_id'])
        for r in rows:
            w.writerow([r['surface'], r['N'], r['k'], r['orientable'], r['n_packings'], r['isom'][0], r['isom'][1],
                        r['aut'][0] if r['aut'] else '', r['aut'][1] if r['aut'] else ''])
    print('%d surfaces, %d problems' % (len(rows), len(problems)))


if __name__ == '__main__':
    main()
