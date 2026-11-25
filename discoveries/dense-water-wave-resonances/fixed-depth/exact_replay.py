#!/usr/bin/env python3
"""Exact dyadic phase enclosures and integer cyclotomic controls."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as R
from math import gcd, isqrt

Q=1<<220

@dataclass(frozen=True)
class D:
    lo: int
    hi: int
    def __add__(a,b): return D(a.lo+b.lo,a.hi+b.hi)
    def __sub__(a,b): return D(a.lo-b.hi,a.hi-b.lo)
    def __mul__(a,b):
        z=[x*y for x in (a.lo,a.hi) for y in (b.lo,b.hi)]
        return D(min(z)//Q,-(-max(z)//Q))
    def __truediv__(a,b):
        assert b.lo>0
        return a*D(Q*Q//b.hi,-(-Q*Q//b.lo))
    def sqrt(a):
        assert a.lo>=0
        u=isqrt(a.hi*Q)
        return D(isqrt(a.lo*Q),u+(u*u!=a.hi*Q))
    def nearest(a):
        lo=(a.lo+Q//2)//Q
        assert (a.hi+Q//2)//Q==lo and a.lo>lo*Q-Q//2
        return lo

def c(x):
    x=R(x)*Q
    return D(x.numerator//x.denominator,-(-x.numerator//x.denominator))

def exp_negative(x):
    assert 0<x<=2
    s=t=R(1)
    for n in range(1,101): t*=-x/n; s+=t
    return D(c(s+t*(-x/101)).lo,c(s).hi)

def tanh_low(h):
    X=exp_negative(2*h)
    return (c(1)-X)/(c(1)+X)

def omega_high(x,b,alpha):
    # For all controls h>.3 and b>=1000, hence rh>128.
    assert b>=1000
    N=x*x+b*b
    return (c(N).sqrt()*c(1+alpha*N)*D(Q-1,Q)).sqrt()

def scaled_phase(w,alpha,m,b):
    return (omega_high(m+1,b,alpha)-omega_high(m,b,alpha)-w)*c(b)/w

def rational_control():
    alpha=R(1,2)
    # h=log(2)/2>.3 since e^(-.6) > 1-.6+.18-.036 > .5.
    assert exp_negative(R(3,5)).lo>Q//2
    w=c(R(1,2)).sqrt()
    A=R(9,4)
    assert A==R(9,4)*alpha/((1+alpha)*R(1,3))
    tally=Counter()
    flips=0
    for m in list(range(100,104))+list(range(1000,1004))+list(range(10000,10004)):
        M=R(2*m+1,2)
        x=A*M*M
        g=x-1/(2*A)
        b=c(g).nearest()
        predicted=(g-b)/2
        assert predicted==(R(49,288) if (2*m+1)**2%16==1 else -R(23,288))
        definite=scaled_phase(w,alpha,m,b)
        assert c(predicted-R(2,b)).hi<definite.lo<=definite.hi<c(predicted+R(2,b)).lo
        assert definite.lo>0 if predicted>0 else definite.hi<0
        tally[predicted]+=1
        old_b=c(x).nearest()
        if b!=old_b:
            assert definite.lo>0>scaled_phase(w,alpha,m,old_b).hi
            flips+=1
    assert tally=={R(49,288):6,-R(23,288):6} and flips==6
    for P in range(1,25):
        for q in range(1,25):
            if gcd(P,q)!=1: continue
            for m in range(q):
                u=P*P*(2*m+1)**2-2*q*q
                assert u%(2*P*q)!=0  # excludes integers and half-integers
    print("rational A=9/4: six exact sign brackets around each limit -23/288 and 49/288; six unshifted-rounding sign reversals")

def irrational_control():
    alpha=R(1,2)
    w=(c(1+alpha)*tanh_low(R(3,5))).sqrt()
    A=c(R(9,4)*alpha)/(w*w)
    shift=c(R(1,2))/A
    bins=[0]*10
    for m in range(100,20100):
        x=A*c(R((2*m+1)**2,4))
        g=x-shift
        b=g.nearest()
        z=scaled_phase(w,alpha,m,b)
        deviation=(g-c(b))*c(R(1,2))
        error=z-deviation
        assert -c(R(2,b)).lo<error.lo<=error.hi<c(R(2,b)).lo
        assert z.hi<0 or z.lo>0
        which=next(i for i in range(10) if z.hi<c(R(-4+i,20)).lo)
        assert z.lo>c(R(-5+which,20)).hi
        bins[which]+=1
    assert bins==[1984,1959,2012,1992,2052,1978,2029,2004,2038,1952]
    print("algebraic h=3/5 alpha=1/2: 20,000 exact phase bins over [-1/4,1/4]:",*bins)

def divide_poly(f,g):
    assert g[-1]==1
    f=f[:]
    degree=len(f)-len(g)
    qu=[0]*max(degree+1,0)
    for i in range(degree,-1,-1):
        qu[i]=v=f[i+len(g)-1]
        for j,w in enumerate(g): f[i+j]-=v*w
    return qu,f

def gauss_control():
    phi=[None]
    for q in range(1,25):
        f=[-1]+[0]*(q-1)+[1]
        for d in range(1,q):
            if q%d==0:
                f,rem=divide_poly(f,phi[d])
                assert not any(rem)
        phi.append(f)
    num=0
    had_tight=False
    for q in (7,8,9,12,16,24):
        for a in range(q):
            residue=[(a*x*x)%q for x in range(q)]
            ker=[d for d in range(q) if 2*a*d%q==0]
            assert len(ker)==gcd(2*a,q)<=2*gcd(a,q)
            if a==1 and q==8: had_tight=(len(ker)==2)
            for j in range(q):
                es=[(r+j*x)%q for x,r in enumerate(residue)]
                diff=[0]*q
                for x in es:
                    for y in es: diff[(x-y)%q]+=1
                for d in ker: diff[es[d]]-=q
                # Evaluate as integer polynomial at primitive q-th root of unity.
                _,rem=divide_poly(diff,phi[q])
                assert not any(rem)
                num+=1
    assert had_tight
    print(f"quadratic Gauss squared identities checked exactly in cyclotomic quotients: {num} odd/even cases")

if __name__=="__main__":
    rational_control(); irrational_control(); gauss_control()
