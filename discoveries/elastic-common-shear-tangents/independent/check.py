#!/usr/bin/env python3
"""Exact rational checks: direct laminate minimization, radial matching, and finite/prior gaps."""
from fractions import Fraction as F
from itertools import product

def norm(E):
    return sum(z*z for row in E for z in row)

def tr(E):
    return sum(E[i][i] for i in range(len(E)))

def inner(K,mu,E,G):
    return K*tr(E)*tr(G)+2*mu*(sum(x*y for a,b in zip(E,G) for x,y in zip(a,b))-tr(E)*tr(G)/len(E))

def add(E,G):
    return tuple(tuple(x+y for x,y in zip(a,b)) for a,b in zip(E,G))

def diag(*vals):
    return tuple(tuple(F(vals[i]) if i==j else F(0) for j in range(len(vals))) for i in range(len(vals)))

def shear(d,a,b):
    return tuple(tuple((F(i==a and j==b)+F(i==b and j==a))/2 for j in range(d)) for i in range(d))

def direct(d,f,K0,K1,mu,E):
    """Return exact finite energy and first 3 Taylor coefficients via direct layered minimization.
    Phase 1 (fraction f) changes shear by t; layers normal to e_(d-1).
    """
    g=1-f
    free0=f*inner(K0,mu,E,E)+g*inner(K1,mu,E,E)
    free1=f*inner(0,1,E,E)
    args=[]
    for i in range(d):
        L=shear(d,i,d-1)
        for j in range(i):
            G=shear(d,j,d-1)
            assert inner(K0,mu,L,G)==inner(K1,mu,L,G)==inner(0,1,L,G)==0
        A=f*inner(K0,mu,L,L)+f*f/g*inner(K1,mu,L,L)
        B=f*inner(0,1,L,L)
        C=f*(inner(K0,mu,L,E)-inner(K1,mu,L,E))
        D=f*inner(0,1,L,E)
        args.append((A,B,C,D))
    def energy(t):
        assert all(A+B*t>0 for A,B,C,D in args)
        return free0+free1*t-sum((C+D*t)**2/(A+B*t) for A,B,C,D in args)
    change0=sum(C*C/A for A,B,C,D in args)
    change1=sum(2*C*D/A-B*C*C/A**2 for A,B,C,D in args)
    change2=sum((D-C*B/A)**2/A for A,B,C,D in args)
    return energy,(free0-change0,free1-change1,-change2)

def tracefree(d):
    D=[diag(1,-1)] if d==2 else [diag(1,-1,0),diag(1,1,-2)]
    return D+[shear(d,i,j) for i in range(d) for j in range(i+1,d)]

def laminate_test(d,f,K0,K1,mu):
    v=F(d-1,d); g=1-f; h=diag(*([F(1,d)]*d))
    aa=2*mu*v; w0=1/(K0+aa); w1=1/(K1+aa)
    w=f*w0+g*w1; D=w0-w1; beta=D/w
    dim=(d-1)*(d+2)//2
    T=v*f*g; J=g*T; R2=v*(f*g)**2
    H=diag(*([F(-1,d)]*(d-1)+[v]))
    B,Bjet=direct(d,f,K0,K1,mu,h)
    assert Bjet[:2]==(1/w-aa,2*beta**2*J)
    Sb=[]
    for E in tracefree(d):
        SE,jet=direct(d,f,K0,K1,mu,E)
        _,mixed=direct(d,f,K0,K1,mu,add(h,E))
        dot=sum(x*y for a,b in zip(H,E) for x,y in zip(a,b))
        assert (mixed[1]-jet[1]-Bjet[1])/2==2*beta*f*g*dot
        Sb.append((SE,tuple(x/(2*norm(E)*dim) for x in jet),2*norm(E)*dim))
    S0,S1,S2=map(sum,zip(*(q[1] for q in Sb)))
    expect=(2*(d-1)*f*g/mu+4*w1*T+4*D*J-4*D*D*R2/w)/dim
    assert S0==mu and S1==f and -2*S2==expect
    if D:
        uncorrected=expect+4*D*D*R2/(w*dim)
        assert -2*S2<uncorrected
    for t in (-mu/2,mu/7,mu/2,2*mu):
        A=(B(t)-Bjet[0])/t
        S=sum(fun(t)/n for fun,_,n in Sb)
        loss=(mu+f*t-S)/(t*t)
        factors=sorted([F(1),1/(1+t/mu)])
        assert factors[0]*Bjet[1]<=A<=factors[1]*Bjet[1]
        assert factors[0]*(-S2)<=loss<=factors[1]*(-S2)
    return Bjet,-2*S2

def radial(d,f,coreK,coremu,shellK,shellmu):
    """Unit ball hydrostatic displacement is Ax inside, (Bx+C|x|^-d)x outside.
    Solve its exact traction equations in volume coordinates.
    """
    a=2*shellmu*F(d-1,d)
    q=f*(shellK-coreK)/(coreK+a)
    B=F(1,d)/(1+q)
    C=B*q
    A=B+C/f
    assert B+C==F(1,d)
    assert coreK*A==shellK*B-a*C/f
    assert coremu>0 # constant core displacement gradient is hydrostatic, no shear
    bulk=d*(shellK*B-a*C)
    return bulk

def H(m,n,s,f):
    return 1/(f/(m+s)+(1-f)/(n+s))-s

if __name__=='__main__':
    for d,f,(k,l),mu in product((2,3),(F(1,4),F(1,2),F(8,13)),
                                  ((F(1,10),F(10)),(F(8),F(1,3)),(F(7,2),F(7,2))),
                                  (F(1),F(7,5))):
        laminate_test(d,f,k,l,mu)
    b,sc=laminate_test(2,F(1,2),F(1,10),F(10),F(1))
    print('36 exact two-phase laminate controls; includes both signs of contrast and omitted-correction negative controls.')
    print('Anisotropic layer fixture:',b,sc)

    for d,f,u,v,m,n in product((2,3),(F(1,2),F(5,7)),(F(1,10),F(3)),
                                (F(7,3),F(4)),(F(1),F(3,2)),(F(2),)):
        assert radial(d,f,u,m,v,n)==H(u,v,2*n*F(d-1,d),f)
    low,high,f,mu,t=F(1,10),F(10),F(1,2),F(1),F(1,2)
    top=radial(2,f,high,mu,low,mu+t)
    bottom=radial(2,f,low,mu+t,high,mu)
    assert (bottom,top)==(F(1),F(343,262))
    w0,w1=1/(low+mu),1/(high+mu); D=w0-w1; w=(w0+w1)/2
    beta=D/w
    T=f*(1-f)/2
    Jfloor=(top-bottom)/(2*t*beta**2)
    assert T/2<Jfloor<=T
    background=f*(1-f)/(2*mu)+w1*T
    curvature_floor=background+D*Jfloor-D*D*(1-f)*(T-Jfloor)/w
    sup=mu+f*t-t*t*curvature_floor/(1+t/mu)
    inf=mu+f*t-t*t*background
    assert (Jfloor,curvature_floor,sup,inf,inf-sup)==(F(121,1048),F(28461,126808),F(307533,253616),F(107,88),F(841,253616))
    print('32 exact radial controls; finite illustrative gap:',(bottom,top,Jfloor,sup,inf,inf-sup))

    shift=lambda k,m:k*m/(k+2*m)
    M_low=H(mu,mu+t,shift(high,mu),f)
    M_high=H(mu,mu+t,shift(low,mu+t),f)
    assert (M_low,M_high,M_low-M_high)==(F(61,50),F(387,322),F(73,4025))
    assert M_high<sup<inf<M_low
    print('Bogdan edge moment/Jensen comparison is strictly sharper:',(M_low,M_high,M_low-M_high))

    # Confocal core/outer axes with exactly equal squared offsets, independent
    # planar depolarization formula and exact shape inversion on both labels.
    core=(F(5),F(1)); outer=(F(7),F(5))
    assert outer[0]**2-core[0]**2==outer[1]**2-core[1]**2==24
    f=(core[0]*core[1])/(outer[0]*outer[1]); g=1-f
    L=lambda ab:(ab[1]/sum(ab),ab[0]/sum(ab))
    M=tuple((a-f*b)/g for a,b in zip(L(core),L(outer)))
    assert (f,M)==(F(1,7),(F(1,8),F(7,8)))
    for fraction in (f,g):
        c=(1-fraction)*(M[1]-M[0])
        y=(1-fraction**2-c*c)/(2*fraction*c)
        assert y>1
        assert (c+fraction*y)**2==1+fraction**2*(y*y-1)
    R2=(f*g)**2 * sum((i-F(1,2))**2 for i in M)
    T=f*g/2
    lower,upper=R2/f,T-R2/g
    assert 0<lower<upper<T
    assert (lower+upper)/2>=R2/f
    print('Exact anisotropic confocal control (inner 5,1; outer 7,5):', (f,M,R2,lower,upper))

    # Implicit shell coordinate derivative cancellation, expressed through
    # squared coordinates so no irrational arithmetic is used.
    for axes2,lam,shares in (
            ((F(25),F(1)),F(2),(F(1,4),F(3,4))),
            ((F(25),F(25),F(1)),F(13),(F(1,7),F(2,7),F(4,7))),
            ((F(2),F(3),F(7),F(11)),F(4),(F(1,5),F(2,5),F(1,5),F(1,5)))):
        x2=[(a+lam)*b for a,b in zip(axes2,shares)]
        assert sum(x/(a+lam) for x,a in zip(x2,axes2))==1
        terms=[x/(a+lam)**2 for x,a in zip(x2,axes2)]
        assert sum(2*v/sum(terms) for v in terms)==2
    print('Confocal shell endpoint-derivative cancellation checked in dimensions 2, 3, 4.')
