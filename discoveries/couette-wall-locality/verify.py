"""Exact rational controls for the wall-dominance theorem; Python standard library."""
from fractions import Fraction as Q
from math import factorial


def exp_lower(x, n):
    assert x >= 0
    return sum((x**j / factorial(j) for j in range(n+1)), Q(0))


def exp_upper(x, n):
    assert 0 <= x < n+2
    first_omitted = x**(n+1) / factorial(n+1)
    return exp_lower(x,n) + first_omitted / (1-x/(n+2))


def gamma_integral(j, n=54):
    """Enclose integral_0^infty u^j exp(-u^3)du for j=0,1."""
    assert j in (0,1) and n >= 8
    # Integral over [0,2] of Taylor expansion, Lagrange residual alternating
    # in its sign for positive argument even before the terms decrease.
    terms = [Q((-1)**k * 2**(3*k+j+1), (3*k+j+1)*factorial(k))
             for k in range(n+1)]
    a, b = sum(terms), sum(terms[:-1])
    # u^j <= u^2 / 2^(2-j) on [2,infty); exp(8)>2000 via terms.
    assert exp_lower(Q(8),20) > 2000
    return (min(a,b), max(a,b)+Q(1,2000*3*2**(2-j)))


def airy_c():
    # Ai(x)/Ai(0) = F0(x) - c F1(x).
    # c = 3^(1/3)*Gamma(2/3)/Gamma(1/3)
    a0,b0=gamma_integral(0)
    a1,b1=gamma_integral(1)
    cube_lo,cube_hi=Q(14422,10000),Q(14423,10000)
    assert cube_lo**3<3<cube_hi**3 and a0>0 and a1>0
    return cube_lo*a1/b0,cube_hi*b1/a0


def basis(x, r, derivative=False, n=18):
    """F_r nominal Airy basis coefficient has degrees r+3i, r=0 or 1."""
    assert r in (0,1) and abs(x) <= Q(5,2)
    a, total = Q(1),Q(0)
    for i in range(n+2):
        deg=r+3*i
        term= (a*deg*x**(deg-1) if deg else Q(0)) if derivative else a*x**deg
        if i<=n:
            total+=term
        else:
            assert abs(x)**3 < Q(1,2)*(deg if derivative else deg+3)*(deg+2)
            return (total-2*abs(term),total+2*abs(term))
        a/= (deg+2)*(deg+3)


def times(a,b):
    z=[x*y for x in a for y in b]
    return min(z),max(z)


def airy(x, derivative=False):
    lo,hi=basis(x,0,derivative)
    u,v=times(airy_c(),basis(x,1,derivative))
    return lo-v,hi-u


def main():
    a,b=airy_c()
    assert Q(728,1000)<a<b<Q(730,1000)
    # All signs exact; any root inside is a sufficient Airy witness.
    assert airy(-Q(102,100),True)[0] > 0
    assert airy(-Q(100,100),True)[1] < 0
    assert airy(-Q(24,10))[1] < 0
    assert airy(-Q(23,10))[0] > 0
    print("exact roots by opposite signs: Neumann beta in (1.00,1.02), Dirichlet (2.30,2.40)")
    print("scalar Airy initial slope coefficient is in (0.728,0.730)")
    t,L,d=Q(3),Q(128),Q(12)
    beta_upper=Q(102,100)
    wall_lb = Q(216,1000)
    bulk_ub = Q(106,1000)
    assert exp_upper(beta_upper*t/2,22) < 1/wall_lb
    assert exp_lower(t**3/12,22) > 1/bulk_ub
    # Universal explicit error: pi < 22/7 and exp(v) > conservative finite sum.
    epsilon_ub = Q(22,7)**2*t/(L-2*d)**2 + 16/exp_lower(d**2/(4*t),25)
    assert epsilon_ub < Q(3,1000)
    assert wall_lb-epsilon_ub > Q(213,1000)
    echo_ratio_lb = (Q(213,1000)/bulk_ub)**2
    assert echo_ratio_lb > Q(403,100)
    print("t=3, L=128, d=12: wall > 0.216, bulk < 0.106, finite-width error < 0.003")
    print("finite Neumann channel > 0.213; echo norm ratio > %s > 4.03" % echo_ratio_lb)
    print("PASS")


if __name__=="__main__":
    main()
