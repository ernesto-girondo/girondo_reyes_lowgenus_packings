"""Copy everything printed by a script to a file in groups/, so that the output can be compared with the
distributed one (e.g. `diff`)."""
import sys


class _Tee:
    def __init__(self, path):
        self.f = open(path, 'w')
        self.out = sys.stdout

    def write(self, s):
        self.out.write(s)
        self.f.write(s)
        self.f.flush()

    def flush(self):
        self.out.flush()
        self.f.flush()


def tee(path):
    sys.stdout = _Tee(path)
