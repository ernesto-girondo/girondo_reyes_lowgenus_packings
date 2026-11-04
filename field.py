"""Exact arithmetic for the extended triangle group (2,3,N) acting on the upper half plane.

Conventions (same as the Sage notebooks of Girondo-Reyes):
  u = e^R,  cosh R = 1/(2 sin(pi/N)),  R = distance from the centre of the regular
  N-gon (angles 2pi/3) to the midpoint of an edge.
  Matrices (acting on the upper half plane, orientation reversing ones act on conj(z)):
     a = [[1,0],[0,-1]]          reflection in the imaginary axis
     b = [[0,1/u],[u,0]]         reflection in the circle |z| = 1/u   (an edge of the polygon)
     c = [[-cs,sn],[sn,cs]]      reflection in a line through i at angle pi/N (through a vertex)
  with sn = sin(pi/N) = u/(u^2+1), cs = cos(pi/N).
  Relations: a^2=b^2=c^2=(ab)^2=(bc)^3=(ca)^N=1.
The field L = Q(u) or Q(u, cos(pi/N)) is handled through PARI; absolute primitive element 'y'.
"""
import cypari2
import mpmath as mp

pari = cypari2.Pari()
pari.allocatemem(2**30)
pari.set_real_precision(80)
mp.mp.dps = 40


def pmpf(x):
    return mp.mpf(str(x).replace(' ', ''))


class Field:
    def __init__(self, N):
        self.N = N
        # minimal polynomial of c0 = cos(pi/N), in variable x
        mc = pari('minpoly(Mod((x^2+1)/(2*x), polcyclo(%d, x)), x)' % (2 * N))
        # u satisfies (1-c0^2)(u^2+1)^2 - u^2 = 0 ; eliminate x, keep y (=u)
        R = pari('polresultant((1-x^2)*(y^2+1)^2 - y^2, %s, x)' % mc)
        s_num = mp.sin(mp.pi / N)
        c_num = mp.cos(mp.pi / N)
        u_num = (1 + mp.sqrt(1 - 4 * s_num ** 2)) / (2 * s_num)
        fac = pari.factor(R)
        best = None
        for i in range(int(pari.matsize(fac)[0])):
            f = fac[i, 0]
            if pari.poldegree(f) < 1:
                continue
            val = abs(pmpf(pari.subst(f, 'y', pari(str(u_num)))))
            if best is None or val < best[0]:
                best = (val, f)
        assert best[0] < mp.mpf(10) ** (-30), "could not identify the minimal polynomial of u"
        red = pari.polredbest(best[1], 1)
        Pu = red[0]
        u_in_Pu = pari('Mod(%s, %s)' % (pari.lift(red[1]), Pu))
        # numeric root of Pu corresponding to u: evaluate candidates
        roots = [pmpf(r.real()) for r in pari.polroots(Pu, precision=400) if abs(pmpf(r.imag())) < mp.mpf(10) ** -20]
        tnum = None
        for r in roots:
            if abs(pmpf(pari.subst(pari.lift(u_in_Pu), 'y', pari(str(r)))) - u_num) < mp.mpf(10) ** -30:
                tnum = r
        assert tnum is not None
        nf = pari.nfinit(Pu)
        cos2 = 1 - u_in_Pu ** 2 / (u_in_Pu ** 2 + 1) ** 2
        rel = pari('x^2') - cos2
        roots2 = pari.nfroots(nf, rel)
        if len(roots2) > 0:
            self.P = Pu
            self.ext = False
            self.t_num = tnum
            self.u = u_in_Pu
            cand = [pari('Mod(%s, %s)' % (pari.lift(r), Pu)) for r in roots2]
            vals = [self.num(cc) for cc in cand]
            self.cs = cand[0] if abs(vals[0] - c_num) < abs(vals[1] - c_num) else cand[1]
        else:
            self.ext = True
            # primitive element theta = u + k*cos(pi/N); minimal polynomial by resultant
            # Pu0 = minimal polynomial of u itself (before polredbest)
            Pu0 = best[1] / pari.pollead(best[1])
            found = None
            for k in (1, 2, 3):
                Rt = pari.polresultant(pari.subst(Pu0, 'y', pari('y') - k * pari('x')), mc, 'x')
                th_num = u_num + k * c_num
                fac2 = pari.factor(Rt)
                for i in range(int(pari.matsize(fac2)[0])):
                    f = fac2[i, 0]
                    if pari.poldegree(f) < 1:
                        continue
                    if abs(pmpf(pari.subst(f, 'y', pari(str(th_num))))) < mp.mpf(10) ** -30 and int(pari.poldegree(f)) == 2 * int(pari.poldegree(Pu0)):
                        found = (k, f)
                        break
                if found:
                    break
            assert found is not None, "no primitive element found"
            k, f = found
            red2 = pari.polredbest(f, 1)
            Pabs = red2[0]
            theta = pari('Mod(%s, %s)' % (pari.lift(red2[1]), Pabs))      # theta inside Q[y]/Pabs
            self.P = Pabs
            # numeric root of Pabs matching theta_num
            self.t_num = None
            for r in pari.polroots(Pabs, precision=400):
                if abs(pmpf(r.imag())) > mp.mpf(10) ** -20:
                    continue
                self.t_num = pmpf(r.real())
                if abs(self.num(theta) - th_num) < mp.mpf(10) ** -25:
                    break
            assert abs(self.num(theta) - th_num) < mp.mpf(10) ** -25
            nf2 = pari.nfinit(Pabs)
            uroots = pari.nfroots(nf2, pari.subst(Pu0, 'y', 'x'))
            self.u = None
            for r in uroots:
                cand = pari('Mod(%s, %s)' % (pari.lift(r), Pabs))
                if abs(self.num(cand) - u_num) < mp.mpf(10) ** -25:
                    self.u = cand
            assert self.u is not None
            self.cs = (theta - self.u) / k
            assert self.cs ** 2 == 1 - self.u ** 2 / (self.u ** 2 + 1) ** 2
        self.sn = self.u / (self.u ** 2 + 1)
        assert abs(self.num(self.sn) - s_num) < mp.mpf(10) ** -25
        assert abs(self.num(self.cs) - c_num) < mp.mpf(10) ** -25
        self.deg = int(pari.poldegree(self.P))
        self.zero = pari('Mod(0, %s)' % self.P)
        self.one = pari('Mod(1, %s)' % self.P)
        one, zero, u, cs, sn = self.one, self.zero, self.u, self.cs, self.sn
        self.A = pari.matrix(2, 2, [one, zero, zero, -one])
        self.B = pari.matrix(2, 2, [zero, 1 / u, u, zero])
        self.C = pari.matrix(2, 2, [-cs, sn, sn, cs])
        self.gens = {'a': self.A, 'b': self.B, 'c': self.C}
        self.I = pari.matrix(2, 2, [one, zero, zero, one])
        self.u_num, self.cs_num, self.sn_num = u_num, c_num, s_num
        self.gens_num = {g: self.mat_num(m) for g, m in self.gens.items()}
        self.R = mp.acosh(1 / (2 * s_num))                                   # inradius
        self.H = mp.acosh(1 / (mp.tan(mp.pi / N) * mp.tan(mp.pi / 3)))      # circumradius r_c
        ab, bc, ca = self.A * self.B, self.B * self.C, self.C * self.A
        assert self.is_scalar(ab * ab) and self.is_scalar(bc ** 3) and self.is_scalar(ca ** N)

    def num(self, e):
        return pmpf(pari.subst(pari.lift(e), 'y', pari(str(self.t_num))))

    def mat_num(self, M):
        return [[self.num(M[0, 0]), self.num(M[0, 1])], [self.num(M[1, 0]), self.num(M[1, 1])]]

    def is_scalar(self, M):
        return M[0, 1] == 0 and M[1, 0] == 0 and M[0, 0] == M[1, 1]

    def word_matrix(self, word):
        M = self.I
        for ch in word:
            M = M * self.gens[ch]
        return M

    def det(self, M):
        return M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
