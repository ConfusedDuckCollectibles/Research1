#!/usr/bin/env python3
"""Exact arithmetic replay of polynomial and finite arithmetic checks."""
from itertools import product
from math import gcd, isqrt

if not __debug__:
    raise RuntimeError("Run without python -O")

class Poly(dict):
    """Sparse exact polynomials with integral coefficients in two variables."""
    @staticmethod
    def aspoly(x):
        return x if isinstance(x, Poly) else Poly({(0, 0): x} if x else {})

    def __add__(self, other):
        out = Poly(self)
        for e, a in Poly.aspoly(other).items():
            out[e] = out.get(e, 0) + a
            if not out[e]:
                del out[e]
        return out

    __radd__ = __add__

    def __neg__(self):
        return Poly({e: -a for e, a in self.items()})

    def __sub__(self, other):
        return self + -Poly.aspoly(other)

    def __rsub__(self, other):
        return Poly.aspoly(other) + -self

    def __mul__(self, other):
        out = Poly()
        for (i, j), a in self.items():
            for (k, l), b in Poly.aspoly(other).items():
                e = i+k, j+l
                out[e] = out.get(e, 0) + a*b
                if not out[e]:
                    del out[e]
        return out

    __rmul__ = __mul__

    def __pow__(self, n):
        assert type(n) is int and n >= 0
        out = Poly({(0, 0): 1})
        for _ in range(n):
            out = out * self
        return out


def common(p, q):
    return p*p+q*q, p*p-q*q, p*q


def formulas(p, q):
    """Use the expanded factors and integer lengths in arXiv:2401.06784 §7.4."""
    p2, q2 = p*p, q*q
    p4, q4 = p2*p2, q2*q2
    P, D, U = common(p,q)
    cplus = p4+2*p**3*q+2*p2*q2-2*p*q**3+q4
    cminus = p4-2*p**3*q+2*p2*q2+2*p*q**3+q4
    fprod = cplus*cminus
    long = p**8+8*p**6*q2-2*p4*q4+8*p2*q**6+q**8
    a_minus = p4-2*p**3*q-2*p2*q2-2*p*q**3+q4
    a_plus  = p4+2*p**3*q-2*p2*q2+2*p*q**3+q4
    b_minus = p4-2*p**3*q+6*p2*q2-2*p*q**3+q4
    b_plus  = p4+2*p**3*q+6*p2*q2+2*p*q**3+q4
    wide_minus = p4-4*p2*q2-q4
    wide_plus = p4+4*p2*q2-q4
    r1 = -4*p2*q2*wide_minus*wide_plus*(3*p4+q4)*(p4+3*q4)
    r2 = D*D*a_minus*b_minus*a_plus*b_plus
    I  = (D*fprod, 4*U*(p4-q4)**2, r1, D*long, P*a_minus*a_plus, 2*U*fprod, P*fprod)
    II = (2*U*fprod, 8*U*U*D*P*P, r2, 2*U*long, P*wide_minus*wide_plus, D*fprod, P*fprod)
    return I, II, fprod, (a_minus*a_plus, b_minus*b_plus)


def reduced(P, D, U):
    F = P**4-4*U*U*D*D
    Fm, Fp = F-8*U*U*P*P, F+8*U*U*P*P
    return F**2-4*D**4*P**4, Fm*Fp


def verify_poly():
    p, q = Poly({(1,0):1}), Poly({(0,1):1})
    P,D,U = common(p,q)
    I,II,F0, factors0 = formulas(p,q)
    F=P**4-4*U*U*D*D
    assert F == F0
    A, B = F-8*U*U*P*P, F+8*U*U*P*P
    assert factors0 == (A,B)
    assert B == P**4+P*P*(2*U)**2+(2*U)**4
    R1, R2 = reduced(P,D,U)
    assert I[2] == (2*U)**2*R1 and II[2] == D*D*R2
    for a,b,c2,d,e,f,g in (I, II):
        assert a*a+b*b == d*d
        assert a*a+c2 == e*e
        assert b*b+c2 == f*f
        assert a*a+b*b+c2 == g*g
    # Odd–odd normalization, expressed as p=r+s,q=r-s.
    P0,D0,U0=common(p+q,p-q)
    a,b=reduced(P0,D0,U0)
    assert a==256*R2 and b==256*R1
    # Duplication numerator in y^2=x^3+x^2+x, after denominator 4y^2.
    x=p
    assert (3*x*x+2*x+1)**2-4*(x**3+x*x+x)*(1+2*x)==(x*x-1)**2
    print('Exact polynomial identities: both expanded seven-length tuples, '
          'family-I residue form, family-II factors, parity swap, and doubling.')


def verify_finite():
    squares={x*x%16 for x in range(16)}
    vals=set()
    for p,q in product(range(16),repeat=2):
        if (p-q)%2==1:
            vals.add(reduced(*common(p,q))[0]%16)
    assert vals=={13} and 13 not in squares
    print('All 128 opposite-parity residue pairs modulo 16: normalized I = 13.')
    for p in (5,7):
        pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x**3-x*x-x)%p==0]
        assert len(pts)+1==8 and (-48)%p!=0
        print('Good reduction:',p,'affine points',pts,'plus O, total 8.')
    for p,q,j,expected in ((2,1,0,(1443,1800,2405)),(4,1,1,(552968,554880,1175057))):
        a,b,c2,d,e,f,g=formulas(p,q)[j]
        assert (a,b,g)==expected and c2>0 and isqrt(c2)**2!=c2
        assert (a*a+b*b,a*a+c2,b*b+c2,a*a+b*b+c2)==(d*d,e*e,f*f,g*g)
        P,D,U=common(p,q)
        F=P**4-4*U*U*D*D
        assert gcd(F-8*U*U*P*P,F+8*U*U*P*P)==1
        print(f'Source example {p,q}: edges {a},{b},sqrt({c2}); space diagonal {g}.')
    a,b,c=44,117,240
    assert gcd(gcd(a,b),c)==1
    assert (a*a+b*b,a*a+c*c,b*b+c*c)==(125**2,244**2,267**2)
    assert sorted([a*a%8,b*b%8,c*c%8])==[0,0,1]
    print('Ordinary Euler-brick control: squared edges are 0,0,1 mod 8.')


if __name__ == '__main__':
    verify_poly()
    verify_finite()
