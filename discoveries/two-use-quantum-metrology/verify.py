#!/usr/bin/env python3
"""Exact algebra, Fraction certificates, and direct accessible-qubit density/SLD replay."""
from fractions import Fraction as F
from itertools import product
import sympy as S

def scores(s,t,m,h,k):
    return (m*m,
            s*(1-m)**2+(1-s)*(1-h)**2,
            t*(1-m)**2+(1-t)*(1-k)**2,
            s*t*(2-m)**2+(1-s)*t*(2-h)**2+s*(1-t)*(2-k)**2)

def f1(s,y):
    return 4*s*y*(1-y)/(1-y+s*y)

def verify_symbolic():
    r,t,m,y,z=S.symbols("r t m y z")
    s1,s2=r*r,t*t
    u,v=r/(1+r),t/(1+t)
    h1,h2=s1/(1-s1),s2/(1-s2)
    Z=S.factor((v*v-u*u)/(h2-h1))
    Z_expected=(1-r)*(1-t)*(2*r*t+r+t)/((1+r)*(1+t)*(r+t))
    assert S.cancel(Z-Z_expected)==0
    residual=r*r*t*t*(1-r)*(1-t)/((1+r)*(1+t)*(r+t)**2)
    assert S.cancel(((1-Z)/2)**2-v*v-h1*Z-residual)==0
    g=scores(s1,s2,m,1+v,1+u)
    A=s1*s2*(2-m)**2-(s1+s2)*(1-m)**2+m*m
    assert S.cancel(g[3]-g[1]-g[2]+g[0]-A)==0
    assert S.expand(A.subs(t,r)-(((1-r*r)*m+2*r*r)**2-2*r*r))==0
    assert S.cancel((g[0]-g[1]).subs(t,r)-
                    (m*((1-r*r)*m+2*r*r)-2*r*r/(1+r)))==0
    assert S.cancel(4*s1/(1+r)**2-f1(s1,y)-
                    4*s1*(y-1/(1+r))**2/(1-y+s1*y))==0
    assert S.cancel((scores(s1,s2,2*z/(1+z),2,2)[1]-4*z*z/(1+z)**2)*(1+z)**2-
                    (1+2*z-3*z*z-4*z*s1))==0
    edge=scores(S.Integer(1),s2,m,S.Integer(0),S.Rational(3,2))
    assert S.expand(edge[3]-edge[2]-s2*(3-2*m))==0
    assert edge[0].subs(m,S.Rational(3,2))==S.Rational(9,4)
    assert S.expand(edge[2].subs(m,S.Rational(3,2)))==S.Rational(1,4)
    print("Symbolic identities: heterogeneous residual, diagonal and edge contradictions, GHZ and one-call certificates: PASS.")

def meanfrac(s,t,p):
    a,b,c,d=p
    return ((s*b+t*c+2*s*t*d)/(a+s*b+t*c+s*t*d),
            (b+2*t*d)/(b+t*d) if b+t*d else 0,
            (c+2*s*d)/(c+s*d) if c+s*d else 0)

def frac_par(s,t,p):
    return 4*sum(pi*gi for pi,gi in zip(p,scores(s,t,*meanfrac(s,t,p))))

def frac_ad(s,t,p,first,q):
    a,b,c,d=p
    wt,z,ps=((1-s)*(b+d),b+d,t) if first==0 else ((1-t)*(c+d),c+d,s)
    return frac_par(s,t,p)+(wt*(f1(ps,q)-f1(ps,d/z)) if wt else 0)

def idx(bits):
    res=0
    for b in bits:
        res=(res<<1)+b
    return res

STR=list(product(range(2),repeat=4))

def embed(i,A):
    return S.kronecker_product(*[A if j==i else S.eye(2) for j in range(4)])

def direct_density(s,t,p,first,q=None):
    """Physical two sensor and two tag qubits. At theta=0 need rho and d rho."""
    ss=(S.Rational(s),S.Rational(t))
    first=int(first); second=1-first
    ps=(S.Rational(q) if q is not None else None)
    v=S.zeros(16,1)
    for (a,b),pr in zip(product(range(2),repeat=2),(p[0],p[2],p[1],p[3])):
        v[idx((a,b,a,b))]=S.sqrt(pr)
    P=S.diag(0,1)
    G=[embed(i,P) for i in range(2)]
    B=[(embed(i,S.diag(1,S.sqrt(e))),embed(i,S.Matrix([[0,S.sqrt(1-e)],[0,0]])))
       for i,e in enumerate(ss)]
    P_loss=S.diag(*[int(bs[first]!=bs[first+2]) for bs in STR])
    rho=S.zeros(16); drho=S.zeros(16)
    allzero=lambda A: all(S.simplify(f)==0 for f in A)
    for e in range(2):
        v1=B[first][e]*v
        v1d=S.I*G[first]*v1
        assert allzero(P_loss*v1-e*v1)  # Kraus jump readable on sensor-tag parity
        if e==1 and q is not None:
            assert allzero(v1d)  # no first-call phase available in this sector
            w=(v1.conjugate().T*v1)[0].simplify()
            v1=S.zeros(16,1)
            for y in range(2):
                b=[0]*4
                b[first+2]=1  # preserve actual first-loss syndrome
                b[second]=b[second+2]=y
                v1[idx(b)]=S.sqrt(w*(ps if y else 1-ps))
            v1d=S.zeros(16,1)
        for f in range(2):
            vv=B[second][f]*v1
            dd=B[second][f]*v1d+S.I*G[second]*vv
            rho+=vv*vv.conjugate().T
            drho+=dd*vv.conjugate().T+vv*dd.conjugate().T
    rho=rho.applyfunc(S.simplify); drho=drho.applyfunc(S.simplify)
    assert S.trace(rho)==1 and S.trace(drho)==0
    # Recover sectors by the accessible sensor/tag parity, not an environment flag.
    LB=S.zeros(16)
    keys=[(b[0]^b[2],b[1]^b[3]) for b in STR]
    assert all((rho[i,j]==drho[i,j]==0) for i in range(16) for j in range(16) if keys[i]!=keys[j])
    for k in set(keys):
        ii=[i for i in range(16) if keys[i]==k]
        RR=rho.extract(ii,ii); DD=drho.extract(ii,ii)
        w=S.trace(RR)
        if w==0:
            assert allzero(DD)
            continue
        assert w>0 and allzero(RR*RR-w*RR)
        for i,j in product(range(len(ii)),repeat=2):
            LB[ii[i],ii[j]]=2*DD[i,j]/w
    assert allzero((rho*LB+LB*rho)/2-drho)
    return S.simplify(S.trace(rho*LB*LB))

def check(s,t,p,fu,fa,f1st,q):
    assert sum(p)==1
    result=frac_ad(s,t,p,f1st,q)
    assert result==fa and fa>fu
    assert direct_density(s,t,p,0)==S.Rational(frac_par(s,t,p))
    assert direct_density(s,t,p,1)==S.Rational(frac_par(s,t,p))
    assert direct_density(s,t,p,f1st,q)==S.Rational(fa)
    print(f'Exact physical density/SLD witness at ({s}, {t}): parallel upper {fu}, adaptive {fa}, absolute gap {fa-fu}.')
    return result

def witnesses():
    s=t=F(3,5); z=s
    p=(z/(1+z),0,0,1/(1+z))
    u=4*max(scores(s,t,2*z/(1+z),2,2))
    assert u==frac_par(s,t,p)==F(9,4)
    a=check(s,t,p,u,F(2421,992),0,F(9,16))
    assert a/u-1==F(21,248)

    s,t,z=F(16,25),F(81,100),F(18,25)
    p=(z/(1+z),0,0,1/(1+z))
    assert z*z==s*t
    u=4*max(scores(s,t,2*z/(1+z),2,2))
    assert u==frac_par(s,t,p)==F(5184,1849)
    a=check(s,t,p,u,F(1996812,667489),0,F(10,19))
    assert a/u-1==F(387,5776)
    b=frac_ad(s,t,p,1,F(5,9))
    assert direct_density(s,t,p,1,F(5,9))==S.Rational(b)
    assert b-u==F(304,3483) and a-u==F(2916,15523)

    s=t=F(1,4)
    p=(F(2,11),F(3,22),F(3,22),F(6,11))
    u=4*max(scores(s,t,F(491,1000),F(743,500),F(743,500)))
    assert u==F(967669,1000000)
    a=check(s,t,p,u,F(1087,1100),0,F(2,3))
    assert a-u>F(2051,100000)
    assert a/u-1>F(21,1000)
    print("All checks passed; no floating-point optimizer or channel-environment readout.")

if __name__=="__main__":
    verify_symbolic()
    witnesses()
