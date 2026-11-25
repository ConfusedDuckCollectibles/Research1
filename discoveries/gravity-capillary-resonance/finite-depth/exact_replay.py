#!/usr/bin/env python3
"""Independent finite-depth controls; sign assertions use only integer/rational arithmetic."""
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from itertools import permutations, product
from math import isqrt


def mul(P, Q_):
    out = Counter()
    for x, a in P.items():
        for y, b in Q_.items():
            out[tuple(i + j for i, j in zip(x, y))] += a*b
    return {x: a for x, a in out.items() if a}


def rational_enclosure_cubic():
    def p(x):
        return x**13 + 89*x**5 - 89*x**8 - 1

    def T(x, r):
        return (x**r - 1)/(x**r + 1)

    lo, hi = Q('2.41829503'), Q('2.41829504')
    assert p(lo) < 0 < p(hi)
    assert T(lo,5)/T(lo,8) < Q(44,45) < T(hi,5)/T(hi,8)
    for _ in range(70):
        mid = (lo + hi)/2
        if p(mid) < 0:
            lo = mid
        else:
            hi = mid
    # Independent enumeration of the η choice and both orders of gradient/normal ψ.
    # Polynomial variables: λ5,λ8,ρ8/ρ5. The two ψ factors were divided by ρ5².
    waves = ((4,3),(4,-3),(-8,0))
    assert tuple(map(sum,zip(*waves))) == (0,0)
    assert tuple(sum(x*x for x in v) for v in waves) == (25,25,64)
    psi_phases_with_i_removed = (-1,-1,1)
    symbolic = Counter()
    for eta, i, j in permutations(range(3)):
        phase = -psi_phases_with_i_removed[i]*psi_phases_with_i_removed[j]
        dot = sum(a*b for a,b in zip(waves[i],waves[j]))
        n8 = int(i == 2 or j == 2)
        symbolic[(0,0,n8)] += -Q(phase*dot,2)
        symbolic[(2-n8,n8,n8)] += -Q(phase,2)
    assert dict(symbolic) == {(0,0,0):7,(2,0,0):1,(0,0,1):64,(1,1,1):-2}
    ratio = Q(5,8)*Q(44,45)
    assert ratio == Q(11,18)
    # Lower bound uses λ5 upper, not a rounded critical depth.
    l5_max = 5*T(hi,5)
    bound = 7 - 3*l5_max**2 + 128*ratio
    assert bound > Q(92,9)
    with localcontext() as ctx:
        ctx.prec = 30
        dec = (Decimal(lo.numerator)/Decimal(lo.denominator)).ln()/2
    print('triad: exact polynomial signs, frequencies and momentum; critical h ≈', dec)
    print('cubic: six ordered placements give terms', dict(symbolic), '; positive, > 92/9')


def merging_example():
    # Half-product in square roots for nonzero coefficients with sum=0.
    rad = [1,2,3]
    c = [1,1,-2]
    F = {(0,0,0):1}
    for eps in product((-1,1),repeat=2):
        F = mul(F, {(1,0,0):c[0],(0,1,0):eps[0]*c[1],(0,0,1):eps[1]*c[2]})
    assert all(all(v % 2 == 0 for v in a) for a in F)
    F = {tuple(v//2 for v in a):b for a,b in F.items()}
    assert all(sum(a) == 2 for a in F)
    # Substitute (z_i−1)/(z_i+1), clear (z_i+1)^2 and then z_i=z^r_i.
    E = Counter()
    for a,b in F.items():
        Q_ = {(0,0,0):b}
        for i, j in enumerate(a):
            ei = tuple(int(k==i) for k in range(3))
            for k in range(2):
                Q_ = mul(Q_, {(0,0,0): -1 if k < j else 1, ei: 1})
        E.update(Q_)
    M = Counter()
    for a,b in E.items():
        M[sum(x*y for x,y in zip(a,rad))] += b
    M = {a:b for a,b in M.items() if b}
    assert len(M) < len([v for v in E.values() if v])
    assert max(M) == 11 and M[11] == 16
    assert sum(M.values()) == sum(j*b for j,b in M.items()) == 0
    assert all(M.get(j,0) == M.get(12-j,0) for j in range(13))
    print('commensurate control: merged',len([v for v in E.values() if v]),
          'monomials into',len(M),'nonzero terms; leading coefficient 16 at degree 11')


def sqrt_enclose(x, digits=260):
    scale = 10**digits
    l = isqrt(x.numerator*scale**2//x.denominator)
    lo,hi = Q(l,scale),Q(l+1,scale)
    assert lo*lo <= x < hi*hi
    return lo,hi


def tsqrt_bounds(zlo,zhi):
    return (sqrt_enclose((1-zhi)/(1+zhi))[0],
            sqrt_enclose((1-zlo)/(1+zlo))[1])


def e2_bounds(n=160):
    term = out = Q(1)
    for k in range(1,n+1):
        term *= Q(2,k)
        out += term
    return out, out + term*Q(2,n+1)/(1-Q(2,n+2))


def trident():
    kk = ((49,0),(-9,0),(20,15),(20,-15))
    r = [isqrt(sum(x*x for x in v)) for v in kk]
    assert r == [49,9,25,25]
    assert [v*v for v in r] == [sum(x*x for x in v) for v in kk]
    assert tuple(map(sum,zip(*kk[:2]))) == tuple(map(sum,zip(*kk[2:])))
    assert sum(map(isqrt,r[:2])) == sum(map(isqrt,r[2:]))
    elo,ehi = e2_bounds()
    # h=1 is algebraic. h=log(2)/2 is a separate exactly rational z control.
    for label,Elo,Ehi in [('h=1',elo,ehi),('h=log(2)/2',Q(2),Q(2))]:
        for N in [1,2,4,8]:
            def zin(m):
                return Ehi**(-(m*N)),Elo**(-(m*N))
            a,b,c = (tsqrt_bounds(*zin(m)) for m in [49,9,25])
            d0,d1 = 7*a[0] + 3*b[0] - 10*c[1], 7*a[1] + 3*b[1] - 10*c[0]
            assert d0 < d1 < 0
            z0,z1 = zin(9)
            norm0,norm1 = d0/z0,d1/z1
            # Since sqrt(N) cancels, these certify Δ_N/(sqrt(N) exp(-18hN)).
            tol = Q(4,2**(9*N)) if label == 'h=log(2)/2' else Q(1,10**(6*N))
            assert -3-tol < norm0 <= norm1 < -3+tol
        print('trident:',label,', exact negative detuning and leading −3 certified at N=1,2,4,8')


if __name__ == '__main__':
    rational_enclosure_cubic()
    merging_example()
    trident()
