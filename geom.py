"""Numeric hyperbolic geometry (upper half plane, mpmath) for the (2,3,N) tessellation."""
import mpmath as mp
from field import Field, pari, pmpf

mp.mp.dps = 40
EPS = mp.mpf(10) ** -18


# ---------- numeric 2x2 real matrices acting on H (orientation reversing: act on conj z) ----------
def mmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def mdet(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def minv(A):
    d = mdet(A)
    return [[A[1][1] / d, -A[0][1] / d], [-A[1][0] / d, A[0][0] / d]]


def mnorm(A):
    """normalise so that |det| = 1 (projective representative)"""
    d = mp.sqrt(abs(mdet(A)))
    return [[A[0][0] / d, A[0][1] / d], [A[1][0] / d, A[1][1] / d]]


def act(A, z):
    if mdet(A) < 0:
        z = mp.conj(z)
    return (A[0][0] * z + A[0][1]) / (A[1][0] * z + A[1][1])


def hdist(z, w):
    return mp.acosh(1 + abs(z - w) ** 2 / (2 * z.imag * w.imag))


def to_disc(z):
    """Cayley map i -> 0 with the polygon edge-midpoint i/u^2 going to the negative imaginary axis"""
    return 1j * (z - 1j) / (z + 1j)


def proj_equal(A, B, tol=mp.mpf(10) ** -12):
    """A == +-B projectively (both nonzero)"""
    A = mnorm(A)
    B = mnorm(B)
    e1 = max(abs(A[i][j] - B[i][j]) for i in range(2) for j in range(2))
    e2 = max(abs(A[i][j] + B[i][j]) for i in range(2) for j in range(2))
    return min(e1, e2) < tol


class Geometry:
    """Everything that depends only on N."""

    def __init__(self, N):
        self.N = N
        self.F = F = Field(N)
        self.gens = F.gens_num
        self.I = [[mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]]
        self.cot = mp.cot(mp.pi / N)
        self.csc = 1 / mp.sin(mp.pi / N)
        self.u = F.u_num
        # fundamental triangle T0: x>=0, |z|>=1/u, |z+cot|<=csc. vertices i, i/u, A
        # A: intersection of |z|=1/u with |z+cot|=csc, x>0
        # (x+cot)^2+y^2 = csc^2, x^2+y^2 = 1/u^2 -> 2x cot + cot^2 = csc^2 - 1/u^2
        x = (self.csc ** 2 - 1 / self.u ** 2 - self.cot ** 2) / (2 * self.cot)
        y = mp.sqrt(1 / self.u ** 2 - x ** 2)
        self.vA = mp.mpc(x, y)
        self.vB = mp.mpc(0, 1)
        self.vC = mp.mpc(0, 1 / self.u)
        # interior base point: centroid-ish, verified to be inside
        self.q0 = (self.vA + self.vB + self.vC) / 3
        assert self.wall_violations(self.q0) == []
        # check the angle at A is pi/3 (sanity of the conventions)
        # tangent directions at A of the two circles
        t1 = 1j * self.vA                          # tangent of |z|=1/u (rotation of radius vector)
        t2 = 1j * (self.vA + self.cot)             # tangent of circle centred -cot
        ang = abs(mp.acos(abs((t1.conjugate() * t2).real) / (abs(t1) * abs(t2))))
        assert abs(ang - mp.pi / 3) < 1e-20 or abs(ang - 2 * mp.pi / 3) < 1e-20, ang
        self.Rota = [self.word_num('ca' * k + 'ab') for k in range(N)]   # central polygon -> neighbour k

    def word_num(self, word):
        M = self.I
        for ch in word:
            M = mmul(M, self.gens[ch])
        return M

    # ---- reduction to the fundamental triangle ----
    def wall_violations(self, q, tol=mp.mpf(10) ** -25):
        v = []
        if q.real < -tol:
            v.append('a')
        if abs(q) < 1 / self.u - tol:
            v.append('b')
        if abs(q + self.cot) > self.csc + tol:
            v.append('c')
        return v

    def reduce(self, M, maxlen=2000):
        """Return the word w (string) with  M = word_matrix(w) projectively, or None if M is not
        (numerically) in the extended triangle group."""
        q = act(M, self.q0)
        word = []
        while True:
            v = self.wall_violations(q)
            if not v:
                break
            ch = v[0]
            q = act(self.gens[ch], q)
            word.append(ch)
            if len(word) > maxlen:
                return None
        w = ''.join(word)
        W = self.word_num(w)
        if proj_equal(W, M, tol=mp.mpf(10) ** -10):
            return w
        return None

    # ---- tessellation centres (BFS) ----
    def centres(self, D):
        """All centres g(i) (g in the extended triangle group) with d(i, g(i)) <= D, as a list of
        (dist, z, word).  Exact deduplication through the field.  Also returns the sorted list of
        distinct distances <= D with exact tanh^2(d/2)."""
        F = self.F
        H = F.H
        margin = mp.mpf('1e-6')
        from collections import deque
        start = ('', self.I, F.I)
        seen = {}
        queue = deque()

        def key_of(Mx):
            A, B, C, Dd = Mx[0, 0], Mx[0, 1], Mx[1, 0], Mx[1, 1]
            det = A * Dd - B * C
            den = C * C + Dd * Dd
            re = (A * C + B * Dd) / den
            im = det / den          # for det=-1 the point is M(conj i) = M(-i): imag part -det/den
            if F.num(det) < 0:
                im = -im
            return str(re) + '|' + str(im), re, im

        k, re, im = key_of(F.I)
        seen[k] = (mp.mpf(0), mp.mpc(0, 1), '', re, im)
        queue.append(('', self.I, F.I))
        results = []
        while queue:
            w, Mn, Mx = queue.popleft()
            z = act(Mn, mp.mpc(0, 1))
            d = hdist(mp.mpc(0, 1), z)
            for kk in range(self.N):
                w2 = w + 'ca' * kk + 'ab'
                Mn2 = mmul(Mn, self.Rota[kk])
                z2 = act(Mn2, mp.mpc(0, 1))
                d2 = hdist(mp.mpc(0, 1), z2)
                if d2 > D + H + margin:
                    continue
                Mx2 = Mx * F.word_matrix('ca' * kk + 'ab')
                key, re, im = key_of(Mx2)
                if key in seen:
                    continue
                seen[key] = (d2, z2, w2, re, im)
                queue.append((w2, Mn2, Mx2))
        # distances
        dists = {}
        for key, (d, z, w, re, im) in seen.items():
            if d <= D + margin:
                results.append((d, z, w))
                # tanh^2(d/2) = |z-i|^2/|z+i|^2 exactly
                t2 = (re * re + (im - 1) ** 2) / (re * re + (im + 1) ** 2)
                ks = str(t2)
                if ks not in dists:
                    dists[ks] = (d, t2)
        dl = sorted(dists.values(), key=lambda t: t[0])
        # drop the zero distance
        dl = [t for t in dl if t[0] > margin]
        results.sort(key=lambda t: t[0])
        return results, dl
