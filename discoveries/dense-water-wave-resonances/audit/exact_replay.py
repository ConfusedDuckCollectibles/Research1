#!/usr/bin/env python3
"""Integer dyadic enclosures + direct DNO Fourier convolution; imports no author code."""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt

Q = 1 << 220

def upper(n, d):
    return -((-n)//d)

@dataclass(frozen=True)
class D:
    lo: int
    hi: int
    def __add__(a,b):
        return D(a.lo+b.lo,a.hi+b.hi)
    def __neg__(a):
        return D(-a.hi,-a.lo)
    def __sub__(a,b):
        return a+-b
    def __mul__(a,b):
        ends=[x*y for x in (a.lo,a.hi) for y in (b.lo,b.hi)]
        return D(min(ends)//Q, upper(max(ends),Q))
    def __truediv__(a,b):
        assert b.lo>0
        return a*D(Q*Q//b.hi,upper(Q*Q,b.lo))
    def sqrt(a):
        assert a.lo>=0
        v=isqrt(a.hi*Q)
        return D(isqrt(a.lo*Q),v+(v*v!=a.hi*Q))

def number(x):
    r=Fraction(x)
    return D(r.numerator*Q//r.denominator, upper(r.numerator*Q,r.denominator))

def box(x,y):
    return D(number(x).lo, number(y).hi)

@lru_cache(None)
def e_negative(x):  # exp(-x), Taylor alternating tail starts far after decreasing.
    assert 0<=x<=2
    t=s=Fraction(1)
    for n in range(1,101):
        t *= -x/n
        s += t
    return number(s+t*(-x/101)).lo,number(s).hi

def tanh(a):
    if a.lo>=128*Q: # 1-tanh(x)<2*e^-256<2^-255<2^-220
        return D(Q-1,Q)
    assert 0<a.lo<=a.hi<=Q
    e=D(e_negative(Fraction(2*a.hi,Q))[0],e_negative(Fraction(2*a.lo,Q))[1])
    return (number(1)-e)/(number(1)+e)

def add(k,l): return tuple(x+y for x,y in zip(k,l))
def neg(k): return tuple(-x for x in k)
def dot(k,l): return sum(x*y for x,y in zip(k,l))

def mode(k,a,h,only_omega=False):
    r2=dot(k,k)
    r=number(r2).sqrt()
    assert r.lo>0
    t=tanh(r*h) if h is not None else D(Q-1,Q)
    L=r*t
    w=(L*number(1+a*r2)).sqrt()
    if only_omega: return w
    E=(L/(number(2)*w)).sqrt()
    return L,w,E,w/L

def dno_cubic(ks,a,h):
    """Evaluate 1/2[psi G1(eta) psi]_0, G1= -G0 eta G0 - div(eta grad)."""
    directed=ks[:2]+(neg(ks[2]),)
    values=[mode(k,a,h) for k in ks]
    L={k:z[0] for d,z in zip(directed,values) for k in (d,neg(d))}
    # Keep only variables a1,a2,bar(a3), as formal coefficient extraction for degree 3.
    # psi entries store coefficient with a common i removed.
    eta=[(k,1<<n,z[2]) for n,(k,z) in enumerate(zip(directed,values))]
    psi=[(k,1<<n,number(-1 if n<2 else 1)*z[2]*z[3]) for n,(k,z) in enumerate(zip(directed,values))]
    G1={}
    for p,i,x in eta:
        for q,j,y in psi:
            if i&j: continue
            k=add(p,q)
            assert k in L
            z=x*y*(number(dot(k,q))-L[k]*L[q])
            G1[k,i|j]=G1.get((k,i|j),number(0))+z
    result=number(0)
    for k,j,y in psi:
        result+=number('-0.5')*y*G1[neg(k),7^j] # i^2=-1
    return result

def phase(a,m,b,h):
    return (mode((m+1,b),a,h,True)-mode((m,b),a,h,True)
            -mode((1,0),a,h,True))

def certify(alpha,m,b,lo,hi,c_min,c_max):
    a=Fraction(alpha)
    ks=((1,0),(m,b),(m+1,b))
    assert add(*ks[:2])==ks[2] and b>0
    assert phase(a,m,b,number(lo)).lo>0>phase(a,m,b,number(hi)).hi
    h=box(lo,hi)
    coefficient=dno_cubic(ks,a,h)
    L,w,E,rho=zip(*(mode(k,a,h) for k in ks))
    normalized=coefficient/(E[0]*E[1]*E[2]*number(b))
    assert number(c_min).hi<normalized.lo<=normalized.hi<number(c_max).lo
    print(f'alpha={alpha} pair=({m},{b}) certified root ({lo},{hi}); {c_min}<C/b<{c_max}; G1 Fourier coefficient positive')

def coupling_at_any_root(a,m,b,h):
    # Here omega_3-omega_2=omega_1 exactly: only call after locating a root.
    k2,k3=(m,b),(m+1,b)
    l1,w1,_,p1=mode((1,0),a,h)
    l2,_,_,p2=mode(k2,a,None)
    l3,_,_,p3=mode(k3,a,None)
    rr=(number(dot(k2,k2))*number(dot(k3,k3))).sqrt()
    q=number(dot(k2,k3))
    high_pair=-(number(b*b)/(rr+q))+(rr-l2*l3)
    return p1*(p2*number(m)+p3*number(m+1)-l1*w1)+p2*p3*high_pair

def count(b):
    # Outside m<4sqrt(b) no partners on I, using the analytic bound in AUDIT.md.
    a=Fraction(1,2)
    lo,hi=number('0.3'),number('0.8')
    wl,wh=(mode((1,0),a,h,True) for h in (lo,hi))
    last=None
    good=[]
    for m in range(4*isqrt(b)):
        highdiff=mode((m+1,b),a,None,True)-mode((m,b),a,None,True)
        if last is not None: assert (highdiff-last).lo>0
        last=highdiff
        test1,test2=highdiff-wl,highdiff-wh
        assert (test1.hi<0 or test1.lo>0) and (test2.hi<0 or test2.lo>0)
        if test1.lo>0 and test2.hi<0: good.append(m)
    assert good==list(range(min(good),max(good)+1))
    assert all(coupling_at_any_root(a,m,b,box("0.3","0.8")).lo>0 for m in good)
    print(f'alpha=1/2 I=[0.3,0.8] b={b}: {len(good)} pairwise-distinct active roots, m={min(good)}..{max(good)}')

def main():
    certify("0.5",100,14105,"0.5999911374537900","0.5999911374537902","1.749999","1.750001")
    certify("2",100,101759,"0.1499997535178385","0.1499997535178388","2.999997","2.999999")
    certify("0.35",103,16446,"0.4000034182185347","0.4000034182185350","1.6249","1.6251")
    for b in (1024,4096,16384,65536): count(b)
    assert phase(Fraction(0),100,14105,number('0.6')).hi<number('-0.5').lo
    print('alpha=0: exact negative phase control')

if __name__=='__main__': main()
