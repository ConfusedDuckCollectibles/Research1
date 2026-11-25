#!/usr/bin/env python3
"""Exact integer/rational interval checks. No floating point in any assertion."""
from fractions import Fraction as F
from math import isqrt
from itertools import permutations
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

def I(x):return F(x),F(x)
def add(X,Y):return X[0]+Y[0],X[1]+Y[1]
def neg(X):return -X[1],-X[0]
def sub(X,Y):return add(X,neg(Y))
def mul(X,Y):
    A=[x*y for x in X for y in Y];return min(A),max(A)
def div(X,Y):
    assert Y[0]>0
    return mul(X,(1/Y[1],1/Y[0]))
def sqrt(X):
    assert X[0]>=0
    D=10**58
    def fl(x):return isqrt(x.numerator*D*D//x.denominator)
    return F(fl(X[0]),D),F(fl(X[1])+1,D)
def square(X):return mul(X,X)

def exp_point(x,N=110):
    assert 0<=x<=2
    S=t=F(1)
    for k in range(1,N+1):
        t*=x/k
        S+=t
    tail=t*x/(N+1)/(1-x/(N+2))
    return S,S+tail
def tanh_positive(X):
    if X[0]>=200:
        # exp(2X)>=2^400 using e>2; 1-tanh(X)<2 exp(-2X).
        return 1-F(1,2**399),F(1)
    assert 0<=X[0]<=X[1]<=1
    e_lo=exp_point(2*X[0])[0];e_hi=exp_point(2*X[1])[1]
    return (e_lo-1)/(e_lo+1),(e_hi-1)/(e_hi+1)
def dot(p,q):return p[0]*q[0]+p[1]*q[1]
def cross(p,q):return p[0]*q[1]-p[1]*q[0]
def mode(k,alpha,h):
    n=dot(k,k);r=sqrt(I(n));T=tanh_positive(mul(r,h))
    lam=mul(r,T);omega=sqrt(mul(lam,I(1+alpha*n)))
    return r,lam,omega,div(omega,lam),sqrt(div(lam,mul(I(2),omega)))
def coeff_by_momenta(ks,phases,lam,rhos):
    # Two phase factors are each (-i) for a or (+i) for a*.
    # Represent +/-i as its sign; product of the signs times -1 is real.
    C=I(0)
    for k,i,j in permutations(range(3)):
        T=sub(I(-dot(ks[i],ks[j])),mul(lam[i],lam[j]))
        C=add(C,mul(I(F(-phases[i]*phases[j],2)),mul(T,mul(rhos[i],rhos[j]))))
    return C
def values(alpha,m,b,h):
    ks=((1,0),(m,b),(m+1,b))
    assert tuple(ks[0][k]+ks[1][k] for k in range(2))==ks[2]
    assert cross(ks[0],ks[1])==b and cross(ks[1],ks[2])==-b
    r,l,w,p,E=zip(*(mode(k,alpha,h) for k in ks))
    phase=sub(sub(w[2],w[1]),w[0]);q=dot(ks[1],ks[2])
    # Stable high-pair term: q-r2r3=-b²/(q+r2r3).
    rr=mul(r[1],r[2])
    small=add(neg(div(I(b*b),add(I(q),rr))),sub(rr,mul(l[1],l[2])))
    lowbr=sub(add(mul(I(m),p[1]),mul(I(m+1),p[2])),mul(l[0],add(w[0],phase)))
    C=add(mul(p[0],lowbr),mul(mul(p[1],p[2]),small))
    raw=coeff_by_momenta((ks[0],ks[1],tuple(-n for n in ks[2])),(-1,-1,1),l,p)
    assert max(raw[0],C[0])<=min(raw[1],C[1])
    full=mul(C,mul(E[0],mul(E[1],E[2])))
    return phase,div(C,I(b)),div(full,sqrt(I(b))),w[0]
def disp(x,prec=16,rounding=ROUND_FLOOR):
    with localcontext() as c:
        c.prec=prec; c.rounding=rounding
        return str(Decimal(x.numerator)/Decimal(x.denominator))
def interval(X,prec=16):return "["+", ".join(disp(x,prec,r) for x,r in zip(X,(ROUND_FLOOR,ROUND_CEILING)))+"]"

def main():
    controls=[
        (F(1,2),100,14105,F(3,5), F("0.5999911374537900"), F("0.5999911374537902"),F("1.7499999"),F("1.7500002")),
        (F(2),100,101759,F(3,20), F("0.1499997535178385"), F("0.1499997535178388"),F("2.999997"),F("2.999999")),
    ]
    for alpha,m,b,h0,lo,hi,clo,chi in controls:
        wm=mode((1,0),alpha,I(h0))[2]
        rb=div(I(F(9,4)*alpha*(m+F(1,2))**2),square(wm))
        assert b-F(1,2)<rb[0]<=rb[1]<b+F(1,2)
        F0=values(alpha,m,b,I(lo))[0]; F1=values(alpha,m,b,I(hi))[0]
        assert F0[0]>0>F1[1]
        _,C,canonical,_=values(alpha,m,b,(lo,hi))
        assert 0<clo<C[0]<C[1]<chi and canonical[0]>0
        print("alpha, m, b, target:",alpha,m,b,h0)
        print(" rounding argument:",interval(rb))
        print(" root bracket:",interval((lo,hi)))
        print(" phase signs + / -:",interval(F0),interval(F1))
        print(" normalized C/b:",interval(C))
        print(" canonical cubic / sqrt(b):",interval(canonical))
    # A zero-tension parameter check, plus the analytic gravity bound in proof.md.
    Fg=values(F(0),100,14105,I(F(3,5)))[0]
    assert Fg[1]<-F(1,2)
    print("pure gravity finite-depth phase control:",interval(Fg))
    print("PASS: all reported signs, rounding choices and lower coefficients have exact rational enclosures.")
if __name__=="__main__":main()
