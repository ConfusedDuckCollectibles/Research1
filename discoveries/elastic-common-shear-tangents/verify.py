#!/usr/bin/env python3
"""Exact fraction checks; no packages, input data, fitting, FEM, or network."""
from fractions import Fraction as F


class Series:
    """Constant, linear and quadratic coefficients in Q[t]/(t^3)."""
    def __init__(self, *c):
        self.c = tuple(map(F, c)) + (F(0),) * (3 - len(c))

    @classmethod
    def lift(cls, x):
        return x if isinstance(x, cls) else cls(x)

    def __add__(self, other):
        return Series(*(a + b for a, b in zip(self.c, self.lift(other).c)))

    __radd__ = __add__

    def __neg__(self):
        return Series(*(-x for x in self.c))

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        a, b = self.c, self.lift(other).c
        return Series(*(sum(a[j] * b[i-j] for j in range(i + 1)) for i in range(3)))

    __rmul__ = __mul__

    def inv(self):
        a, b, c = self.c
        return Series(1/a, -b/a**2, b*b/a**3-c/a**2)

    def __truediv__(self, other):
        return self * self.lift(other).inv()

    def __rtruediv__(self, other):
        return self.lift(other) / self


def theorem(d, f, mu, k1, k2, j, r2):
    """Moment formulas, r2=||R||^2."""
    m = F(d*(d+1), 2) - 1
    v = 1 - F(1, d)
    a = 2*mu*v
    w1, w2 = 1/(k1+a), 1/(k2+a)
    h = f*w1+(1-f)*w2
    delta = w1-w2
    beta = delta/h
    tot = v*f*(1-f)
    q0 = (2*(d-1)*f*(1-f)/mu+4*w2*tot)/m
    bulk0 = 1/h-a
    b = 2*beta**2*j
    c2 = 4*beta**2*r2
    q = q0+(4*delta*j-4*delta**2*r2/h)/m
    if delta:
        assert m*q == m*q0+(2*h*h/delta)*b-h*c2
        assert c2/(2*f) <= b <= 2*beta**2*tot-c2/(2*(1-f))
    return bulk0, b, c2, q


# Independent exact direct 2D layer minimization, from local elastic
# energy in the three strain components, not the moment perturbation.
def laminate(f, mu, k1, k2):
    p, n = Series(mu, 1), Series(mu)
    ll = (k1+p, k2+n)
    mm = (k1-p, k2-n)
    avg = lambda u, w: f*u+(1-f)*w
    a = avg(1/ll[0], 1/ll[1])
    b = avg(mm[0]/ll[0], mm[1]/ll[1])
    c = avg(ll[0]-mm[0]*mm[0]/ll[0], ll[1]-mm[1]*mm[1]/ll[1])
    inv_mu = avg(1/p, 1/n)

    def energy(xx, yy, xy):
        u = xx+yy*b
        return u*u/a+yy*yy*c+4*xy*xy/inv_mu

    hydro = energy(F(1,2), F(1,2), 0)
    shear = (energy(1,-1,0)+energy(0,0,1))/8
    coupling_on_diag = (energy(F(3,2), F(-1,2), 0)
                        -energy(F(-1,2), F(3,2), 0))/4
    assert coupling_on_diag.c[0] == 0
    assert shear.c[:2] == (mu, f)
    return hydro.c[0], hydro.c[1], coupling_on_diag.c[1]**2/2, -2*shear.c[2]


for f, mu, k1, k2 in (
    (F(1,2), F(1), F(1,10), F(10)),
    (F(3,8), F(7,6), F(5,3), F(3,2)),       # reversed ordering
    (F(9,10), F(1,3), F(43,7), F(1,7)),     # skew share / contrast
    (F(4,7), F(7,3), F(19,3), F(19,3)),     # coincident bulks
):
    tot = f*(1-f)/2
    j, r2 = (1-f)*tot, f*(1-f)*tot
    result = laminate(f, mu, k1, k2)
    assert theorem(2, f, mu, k1, k2, j, r2) == result
    print('direct layer benchmark (B0,b,||c||^2,Q):', *result)


# Finite benchmark including both hard reading errors and finite interfaces.
f, mu, k1, k2, t = F(1,2), F(1), F(1,10), F(10), F(1,2)
w1, w2, h, delta, beta = F(10,11), F(1,11), F(1,2), F(9,11), F(18,11)
tot = F(1,8)
q0 = F(3,11)
q_end = theorem(2,f,mu,k1,k2,tot,0)[3]
assert q_end == F(21,44) and q_end/q0 == F(7,4)


def cap(bulk_observed_lower):
    j_min = (bulk_observed_lower-1)/(t*2*beta**2)
    r2_max = (1-f)*(tot-j_min)
    q_min = theorem(2, f, mu, k1, k2, j_min, r2_max)[3]
    shear_upper = mu+f*t-t*t*q_min/(2*(1+t/mu))
    return q_min, shear_upper


bulk_goal, shear_goal, eps = F(13,10), F(607,500), F(1,4000)
assert cap(bulk_goal) == (F(637,1452), F(21143,17424))
assert cap(bulk_goal-eps)[1] < shear_goal-eps
print('bulk goal 1.30 exact shear cap:', cap(bulk_goal)[1])
print('hard simultaneous reading-error excess margin:', shear_goal-eps-cap(bulk_goal-eps)[1])

assert (F(5,8))**20 < F(1,10000)
shell = mu+t
h_shell = f/(k1+shell)+(1-f)/(k2+shell)
bulk_shell = 1/h_shell-shell
q_shell,q_core = (1/((ki+shell)*h_shell) for ki in (k1,k2))
assert (bulk_shell,q_shell,q_core) == (F(343,262),F(230,131),F(32,131))
stress_shell_sq = 2*(k1*q_shell)**2+2*shell**2*(q_shell-q_core)**2
stress_core_sq = 2*(k2*q_core)**2
assert max(stress_shell_sq,stress_core_sq) < 12
assert max(abs(1/(2*k1)-1/(2*k2)),abs(1/(2*mu)-1/(2*shell))) < 5
eta = F(1,10000)
core_bound = mu+f*t-t*t*(q0+4*delta*eta)/2
shell_bound = bulk_shell-60*eta
assert core_bound == F(107,88)-F(9,22)*eta
assert core_bound-eps > shear_goal and shell_bound-eps > bulk_goal
print('separate finite-layout guaranteed S and B:', core_bound, shell_bound)


# An exact positive, anisotropic ellipse: 2D depolarization factor
# Hessian(U_a)_11 = b/(a+b). Verify physical confocality and core matrix.
ao, bo, ac, bc = F(4), F(2), F(7,2), F(1,2)
assert ao*ao-bo*bo == ac*ac-bc*bc == 12
g = ac*bc/(ao*bo)
m1 = (bc/(ac+bc)-g*bo/(ao+bo))/(1-g)
assert g == F(7,32) and m1 == F(1,15)
print('exact confocal ellipse core shape (g,M eigenvalues):', g, m1, 1-m1)


# 3D positive anisotropic compact confocal ellipsoid with inner squares
# (1,1,4), outer squares (6,6,9), g=1/9. Substitution v=sqrt(4+u)
# reduces M_3 to an integral of a positive decreasing rational function.
# Rational upper/lower sums enclose it, independently of any quadrature.
def kernel(v):
    return F(9,4)/(v*v*(v*v-3))


n = 128
nodes = [2+F(i,n) for i in range(n+1)]
low = sum(map(kernel,nodes[1:]))/n
high = sum(map(kernel,nodes[:-1]))/n
assert F(158,1000) < low < high < F(163,1000)
print('3D confocal rigorous enclosure for small M eigenvalue: (.158,.163)')
print('All checks passed.')
