"""Short writing of words in a, b, c with powers of repeated blocks, for tables and reports:
compress('cacacacbacac') == '(ca)^3cb(ac)^2'. The words are read from left to right, as everywhere in the code."""
from functools import lru_cache


def compress(w, maxblock=6):
    """shortest writing of a word in a, b, c using powers of repeated blocks, e.g. cacacacb -> (ca)^3cb"""
    @lru_cache(None)
    def best(i):
        if i == len(w):
            return ''
        cands = [w[i] + best(i + 1)]
        for L in range(1, maxblock + 1):
            blk = w[i:i + L]
            if len(blk) < L:
                break
            k = 1
            while w[i + k * L:i + (k + 1) * L] == blk:
                k += 1
            for kk in range(2, k + 1):
                tok = ('(%s)' % compress(blk, maxblock) if L > 1 else blk) + '^%d' % kk
                cands.append(tok + best(i + kk * L))
        return min(cands, key=lambda t: (len(t), t))
    return best(0)
