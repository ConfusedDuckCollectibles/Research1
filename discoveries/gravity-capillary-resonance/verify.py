#!/usr/bin/env python3
"""Finite exact-arithmetic controls; no search over resonance graphs."""
from collections import Counter
from fractions import Fraction as F

def sqsplit(n):
    """Return independent square-free part and positive square coefficient."""
    assert n > 0
    residue=n; sf=v=1;p=2
    while p*p <= residue:
        ex=0
        while residue % p == 0:
            residue //= p
            ex += 1
        sf *= p**(ex % 2)
        v *= p**(ex // 2)
        p += 1
    sf *= residue
    assert sf*v*v == n
    return sf,v

def mode(alpha, aspect2, k):
    A,B=alpha.numerator,alpha.denominator
    U,D=aspect2.numerator,aspect2.denominator
    i,j=k
    N=D*i*i+U*j*j
    q,l=sqsplit(N)
    R,m=sqsplit(l*(B*D+A*q*l*l))
    # Four-power-free radical and its un-reduced fourth power.
    S=q*R*R
    assert S*m**4 == N*(B*D+A*N)**2
    # The quartic factorization means omega = m * fourth_root(S) / fourth_root(D) / sqrt(BD).
    c=A*q
    X,Y=c*R*l,c*R*R*m
    C=c*B*D*R*R
    assert Y*Y==X**3+C*X
    assert -64*C**3!=0 and X%(c*R)==0 and Y%(c*R*R)==0
    assert (X//(c*R),Y//(c*R*R))==(l,m)
    return (N,(q,R),l,m,(X,Y,C))

def exact(alpha, aspect2, signed, expected):
    slots=Counter()
    wave=[0,0]
    for eps,k in signed:
        N,label,l,m,point=mode(alpha,aspect2,k)
        slots[label]+=eps*m
        for i in range(2):
            wave[i]+=eps*k[i]
    assert wave==[0,0]
    balanced=all(v==0 for v in slots.values())
    assert balanced is expected
    return slots

def tri(p,q,r):
    return [(1,p),(1,q),(-1,r)]

def main():
    # Noncollinear and crossing the capillary-gravity energy-term threshold.
    a=F(7,1259); d=F(1); triad=tri((12,5),(12,-5),(24,0))
    assert a*169<1<a*576
    assert (12*(-5)-5*12)!=0
    assert exact(a,d,triad,True)=={(1,31746):0}
    T=[mode(a,d,x) for _,x in triad]
    assert [x[3] for x in T]==[1,1,2]
    print('gravity-to-capillary noncollinear triad:',a,T)

    # Class excludes a momentum-conserving harmonic as tension is changed.
    triad=tri((1,0),(1,0),(2,0))
    assert exact(F(1,2),d,triad,True)=={(1,3):0}
    assert exact(F(1),d,triad,False)=={(1,2):2,(1,10):-1}
    print('harmonic:',exact(F(1,2),d,triad,True),'; tuned-off:',exact(F(1),d,triad,False))

    # Rational squared aspect with irrational actual aspect; q on physical radius must include D.
    a=F(7,1978);d=F(1,2);triad=tri((17,12),(17,-12),(34,0))
    opts=exact(a,d,triad,True)
    assert len(opts)==1 and next(iter(opts))[0]==2
    T=[mode(a,d,x) for _,x in triad]
    assert [x[0] for x in T]==[722,722,2312]
    assert [x[2] for x in T]==[19,19,34]
    print('rectangle with squared aspect 1/2:',a,T)

    quat=[(1,(1,1)),(1,(2,-1)),(-1,(1,-1)),(-1,(2,1))]
    assert exact(F(1),F(1),quat,True)=={(2,3):0,(5,6):0}
    print('two allowed sectors in one angular 2+2 resonance:',exact(F(1),F(1),quat,True))

    # The quadratic radical identity shows no parameter-uniform positive gap.
    triad=tri((4,3),(4,-3),(8,0))
    assert exact(F(1),F(1),triad,True)=={(1,130):0}
    for M in (5,50,500):
        a=1+F(1,M)
        assert 4*5*(1+25*a)-8*(1+64*a)==-F(12,M)
        assert exact(a,F(1),triad,False)
    print('near-resonance rational detuning: 12(1-alpha) exact difference of squared totals')
    print('PASS: all validations use exact integer or rational arithmetic.')

if __name__=='__main__':
    main()
