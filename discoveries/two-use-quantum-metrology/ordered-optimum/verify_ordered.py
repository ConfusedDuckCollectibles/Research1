#!/usr/bin/env python3
"""Exact rational QFI witness and directed enclosures; optional numeric author-SDP comparison."""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import argparse

@dataclass(frozen=True)
class I:
    l: Q
    h: Q
    @staticmethod
    def of(x): return x if isinstance(x,I) else I(Q(x),Q(x))
    def __add__(self,other):
        other=I.of(other)
        return I(self.l+other.l,self.h+other.h)
    __radd__=__add__
    def __neg__(self): return I(-self.h,-self.l)
    def __sub__(self,other): return self+(-I.of(other))
    def __rsub__(self,other): return I.of(other)+(-self)
    def __mul__(self,other):
        other=I.of(other)
        v=[a*b for a in (self.l,self.h) for b in (other.l,other.h)]
        return I(min(v),max(v))
    __rmul__=__mul__
    def reciprocal(self):
        assert self.l*self.h>0
        return I(1/self.h,1/self.l)
    def __truediv__(self,other): return self*I.of(other).reciprocal()
    def __rtruediv__(self,other): return I.of(other)*self.reciprocal()
    def square(self):
        return I(0 if self.l<=0<=self.h else min(self.l**2,self.h**2),max(self.l**2,self.h**2))
    def sqrt(self):
        assert self.l>=0
        B=2**115
        l=isqrt(self.l.numerator*B**2//self.l.denominator)
        h=isqrt(self.h.numerator*B**2//self.h.denominator)
        return I(Q(l,B),Q(h if Q(h*h,B*B)==self.h else h+1,B))

def j1(s):
    return 4*s/(1+I.of(s).sqrt()).square()

def quartic_root_point(s,d,bits=80):
    lo,hi=Q(0),Q(1)
    ell=1-s
    for _ in range(bits):
        t=(lo+hi)/2
        f=s-t**4-d*ell*t*(1+t*t)
        if f==0:return I.of(t)
        elif f>0:lo=t
        else:hi=t
    return I(lo,hi)

def detail(a,s,p):
    P=1-(1-a)*p
    d=I.of(a*p*(1-p)).sqrt()/P
    # t decreases with d; evaluate opposite endpoints to control all radicals.
    t=I(quartic_root_point(s,d.h).l,quartic_root_point(s,d.l).h)
    e=t/(1+t.square())
    W=(s+t.square())/(1+t.square())
    V=4*s*(d+e).square()/W
    f=P*V+(1-P)*j1(s)
    dp=(a*(1-2*p))/(2*d*P*P)+(1-a)*d/P
    deriv=(1-a)*(j1(s)-V)+P*(8*s*(d+e)/W)*dp
    return f,deriv,t

def max_bracket(a,s):
    if a==0 or s==0:return j1(s if a==0 else a),I(Q(0),Q(1))
    if a==s==1:return I.of(4),I.of(Q(1,2))
    if a==1:return detail(a,s,Q(1,2))[0],I.of(Q(1,2))
    left,right=Q(0),Q(1)
    for _ in range(32):
        p=(left+right)/2
        _,g,_=detail(a,s,p)
        if g.l>0:left=p
        elif g.h<0:right=p
        else: break
    p=(left+right)/2
    f,g,_=detail(a,s,p)
    # Concavity independently proved: every maximum lies in [left,right];
    # its value lies below the tangent line determined by the exact derivative.
    above=max(Q(0), (g*I(left-p,right-p)).h)
    return I(f.l,f.h+above),I(left,right)

def fmt_box(b,n=12):
    B=10**n
    f=lambda x:f"{x//B}.{x%B:0{n}d}"
    l=(b.l*B).__floor__();h=(b.h*B).__ceil__()
    return "["+f(l)+", "+f(h)+"]"

def rational_control():
    a=Q(16,25);s=Q(81,100);p=Q(1,2);t=Q(9,10)
    P=Q(41,50);d=Q(20,41);e=t/(1+t*t)
    j=4*s/(1+Q(9,10))**2
    val=P*4*s*(d+e)**2*(1+t*t)/(s+t*t)+(1-P)*j
    assert P==1-(1-a)*p and d*d==a*p*(1-p)/(P*P)
    assert val==Q(203724139,66974525)
    assert val/Q(5184,1849)-1>Q(84937,1000000)
    print(f"Exact rational witness J={val}; gain over exact parallel >8.4937%.")

def proofs_checks():
    import sympy as S
    x,z,d,s,t=S.symbols("x z d s t",real=True)
    W=1-(1-s)*x
    total=d*d+s*x+2*s*d*z-(s*x-(1-s)*d*z)**2/W-(1-s)*d*d*z*z/x
    want=d*d+(s*x*(1-x)+2*s*d*z-(1-s)*d*d*z*z/x)/W
    assert S.cancel(total-want)==0
    assert S.cancel(want.subs(z,s*x/((1-s)*d))-(d*d+s*x/(1-s)))==0
    # With z^2=x(1-x), visible branch formula is s*(d+z)^2/W.
    reduced=S.cancel((want-s*(d+z)**2/W)*W*x)
    assert S.rem(reduced,z*z-x*(1-x),z)==0
    R=d*(1+t*t)+t
    D=(1+t*t)*(s+t*t)
    assert S.cancel(S.diff(R,t)*D-R*S.diff(D,t)/2-(s-t**4-d*(1-s)*t*(1+t*t)))==0
    print("Symbolic branch QFI identity, unsaturated maximum, saturation and quartic: PASS.")

def numeric_sdp_reference():
    import numpy as np
    from scipy.optimize import brentq,minimize_scalar
    data=Path(__file__).resolve().parent/"source/QFI_seq_rangep_0_1_N_2.txt"
    q=np.loadtxt(data,delimiter=",")
    assert q.shape==(201,)
    def calc(a,s):
        if a==s==0:return 0.
        if a==s==1:return 4.
        J=4*s/(1+s**.5)**2
        def f(p):
            w=1-(1-a)*p
            d=(a*p*(1-p))**.5/w
            t=brentq(lambda t:s-t**4-d*(1-s)*t*(1+t*t),0,1)
            return w*4*s*(d*(1+t*t)+t)**2/((1+t*t)*(s+t*t))+(1-w)*J
        return -minimize_scalar(lambda p:-f(p),bounds=(0,1),method="bounded",
                                options={"xatol":1e-13}).fun
    analytic=np.array([calc(1-j/200,1-j/200) for j in range(201)])
    diff=np.abs(q-analytic)
    for j in (0,40,80,100,150,180,199,200):
        print(f"SDP equal decay {j/200:.3f}: source {q[j]:.12f}, scalar {analytic[j]:.12f}, absdiff {diff[j]:.3e}")
    j=int(np.argmax(diff[1:-1]))+1
    print(f"Source comparison maximum interior difference {diff[j]:.6e} at decay {j/200:.3f} ({len(q)-2} rows).")
    assert diff[1:-1].max()<1e-8

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--sdp",action="store_true")
    opt=arg.parse_args()
    proofs_checks();rational_control()
    cases=[(Q(16,25),Q(81,100)),(Q(81,100),Q(16,25)),(Q(3,5),Q(3,5)),(Q(1,2),Q(1,2)),(Q(1,4),Q(1)),(Q(1),Q(1,4))]
    boxes=[]
    for a,s in cases:
        box,p=max_bracket(a,s);boxes.append(box)
        assert box.h-box.l<Q(1,10**15)
        print(f"Certified {a} then {s}: J in {fmt_box(box,14)}, p in {fmt_box(p,10)}.")
    gap=boxes[0]-boxes[1]
    assert gap.l>Q(8974639,100000000) and gap.h<Q(897464,10000000)
    print(f"Certified unequal-order gap: {fmt_box(gap,14)}.")
    if opt.sdp:numeric_sdp_reference()
if __name__=="__main__":main()
