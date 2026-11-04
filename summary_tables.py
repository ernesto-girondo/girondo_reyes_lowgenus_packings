"""Human-readable tables from groups/*_results.json:
  groups/surfaces_summary.csv   one row per surface (373 surfaces of the paper and 50 double covers)
  groups/surfaces_summary.xlsx  the same table as a spreadsheet (needs the Python package openpyxl)
  groups/double_covers.txt      the correspondence non-orientable surface -> canonical double cover (Appendix A)

In the tables, the side pairings are written with powers of repeated blocks: (ca)^3cb(ac)^2 is the word
cacacacbacac (read from left to right, as everywhere in the code). The column class_id is a short fingerprint of
the canonical form of the class (field `key` of the results files), which does not depend on the order in which
GAP lists the classes.
"""
import csv, hashlib, json
from wordfmt import compress

NEC = ['nec_N12', 'nec_N9', 'nec_N8', 'nec_N7']
FUCHS = ['fuchs_N18', 'fuchs_N12', 'fuchs_N10', 'fuchs_N9', 'fuchs_N8', 'fuchs_N7']


MAXP = 5          # largest number of packings of a surface


def centres(r):
    """one disc centre of each packing other than the one of the uniformization, padded to MAXP - 1 columns"""
    out = []
    for p in r['packings']:
        if float(p['dist_to_origin']) > 1e-9:
            x, y = (float(t) for t in p['centre_disc'])
            out.append('%.8f%+.8fi' % (x, y))
    assert len(out) <= MAXP - 1
    return out + [''] * (MAXP - 1 - len(out))


def class_id(key):
    return hashlib.sha1(json.dumps(key).encode()).hexdigest()[:10]


def sidepairings(words):
    return ' ; '.join(compress(w) for w in words)


rows = []
fuchs_by_key = {}
for fam in FUCHS + NEC:
    groups = [json.loads(l) for l in open('groups/%s.jsonl' % fam) if l.strip()]
    res = json.load(open('groups/%s_results.json' % fam))
    for g, r in zip(groups, res):
        if fam.startswith('fuchs'):
            fuchs_by_key.setdefault((r['N'], tuple(map(tuple, r['key']))), []).append(r['label'])
        rows.append([r['label'], r['N'], r['k'], r['orientable'], '', r['d'], r['n_packings'], r['n_aut'],
                     r['anticonformal_aut']] + centres(r) + [class_id(r['key']), ' ; '.join(g['gens']),
                                                              sidepairings(r['sidepairings'])])
        if r.get('double_cover'):
            c = r['double_cover']
            rows.append([c['label'], c['N'], c['k'], c['orientable'], r['label'], c['d'], c['n_packings'], c['n_aut'],
                         c['anticonformal_aut']] + centres(c) + [class_id(c['key']), '', sidepairings(c['sidepairings'])])
HEADER = (['surface', 'N', 'k', 'orientable', 'double_cover_of', 'index_in_Delta_pm', 'n_packings',
           'induced_isometries', 'induced_anticonformal'] +
          ['centre_packing_%d' % j for j in range(2, MAXP + 1)] + ['class_id', 'generators_of_K', 'side_pairings'])
with open('groups/surfaces_summary.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(HEADER)
    w.writerows(rows)
try:
    import datetime
    from openpyxl import Workbook
    from openpyxl.styles import Font
    wb = Workbook()
    ws = wb.active
    ws.title = 'surfaces'
    ws.append(HEADER)
    for row in rows:      # booleans as the words TRUE/FALSE, so that they do not change with the language of Excel
        ws.append([('TRUE' if x else 'FALSE') if isinstance(x, bool) else x for x in row])
    for cell in ws[1]:
        cell.font = Font(bold=True)
    ws.freeze_panes = 'B2'
    for col, width in zip('ABCDEFGHIJKLMNOP', [14, 4, 4, 10, 15, 8, 8, 8, 10, 26, 26, 26, 26, 12, 40, 60]):
        ws.column_dimensions[col].width = width
    wb.properties.creator = 'summary_tables.py'
    wb.properties.created = wb.properties.modified = datetime.datetime(2026, 1, 1)
    wb.save('groups/surfaces_summary.xlsx')
except ImportError:
    print('openpyxl is not installed: groups/surfaces_summary.xlsx not written')

lines = ['Canonical double covers (Appendix A). |ind| = order of the group of isometries induced by the triangle group.',
         '']
for fam in NEC:
    res = json.load(open('groups/%s_results.json' % fam))
    N, k = res[0]['N'], res[0]['k']
    lines.append('N=%d: non-orientable %d-extremal surfaces of genus 3 -> orientable %d-extremal surfaces of genus 2'
                 % (N, k, 2 * k))
    covers = set()
    for r in res:
        c = r['double_cover']
        tgt = fuchs_by_key[(N, tuple(map(tuple, c['key'])))]
        assert len(tgt) == 1
        covers.add(tgt[0])
        lines.append('   %-12s packings=%d |ind|=%-2d  ->  %-15s packings=%d |ind|=%d' % (
            r['label'], r['n_packings'], r['n_aut'], tgt[0], c['n_packings'], c['n_aut']))
    res_f = json.load(open('groups/fuchs_N%d_results.json' % N))
    others = [r['label'] for r in res_f if r['anticonformal_aut'] and r['label'] not in covers]
    lines.append('   Fuchsian classes with an induced anticonformal automorphism that are not double covers: %s'
                 % ', '.join(others))
    lines.append('')
open('groups/double_covers.txt', 'w').write('\n'.join(lines))
print('\n'.join(lines))
print('surfaces_summary.csv: %d rows' % len(rows))
