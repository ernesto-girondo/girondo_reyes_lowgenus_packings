"""Compare two copies of the data folder, e.g. the distributed one and a freshly computed one.

    python3 compare_outputs.py groups_distributed groups

Every file of the first folder is compared with the file of the same name in the second one. Running times are
ignored (field "time" of the results files, "time=" in the GAP log) and spreadsheets are compared cell by cell;
everything else must be identical, since the
computation is deterministic. Files present only in the second folder are listed but not treated as differences.
Exit status 0 if no differences were found.
"""
import json, os, re, sys


def strip_time(x):
    if isinstance(x, dict):
        return {k: strip_time(v) for k, v in x.items() if k != 'time'}
    if isinstance(x, list):
        return [strip_time(v) for v in x]
    return x


def load(path):
    if path.endswith('.json'):
        return strip_time(json.load(open(path)))
    if path.endswith('.xlsx'):                       # compare the cells, not the bytes (the file stores dates)
        try:
            import openpyxl
        except ImportError:
            return 'openpyxl not installed'
        return [[c.value for c in row] for row in openpyxl.load_workbook(path).active.iter_rows()]
    text = open(path, encoding='utf-8').read()
    return re.sub(r'time=\d+', 'time=*', text)


def main(a, b):
    files = sorted(f for f in os.listdir(a) if os.path.isfile(os.path.join(a, f)))
    bad, missing = [], []
    for f in files:
        pb = os.path.join(b, f)
        if not os.path.exists(pb):
            missing.append(f)
            continue
        if load(os.path.join(a, f)) != load(pb):
            bad.append(f)
    extra = sorted(f for f in os.listdir(b) if os.path.isfile(os.path.join(b, f)) and f not in files)
    for f in files:
        status = 'MISSING' if f in missing else ('DIFFERENT' if f in bad else 'identical')
        print('%-9s %s' % (status, f))
    if extra:
        print('only in %s: %s' % (b, ', '.join(extra)))
    print('%d files compared: %d identical, %d different, %d missing' % (
        len(files), len(files) - len(bad) - len(missing), len(bad), len(missing)))
    return 0 if not bad and not missing else 1


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
