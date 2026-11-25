#!/usr/bin/env python3
"""Exact independent Legendre-basis continuum calculation; no source imports."""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from math import factorial

N,R=Q(8),Q(7,100)
c1,c2=Q(3,5),Q(6,5)

def step(v,c,N,R):
    """A applied in basis (-i)^k P_k, not a velocity grid."""
    out=[Q(0)]*(len(v)+1)
    for k,w in enumerate(v):
        if k:
            out[k]-=(R + (N if k>=2 else 0))*w
            out[k-1]-=c*k*w/(2*k+1)
        out[k+1]+=c*(k+1)*w/(2*k+1)
    return out

def exact_theta(c,N,R,n=52):
    v=[Q(1)]
    val=Q(1)
    for j in range(1,n+1):
        v=[w/j for w in step(v,c,N,R)]
        val+=v[0]
    err=(N+R+c)**(n+1)/factorial(n+1)
    assert err<Q(1,10**17)
    return val-err,val+err

def altern(c,power_offset,terms=14):
    # sin(c)/c if offset=1, cos(c) if zero. Terms decrease from j=0.
    assert c*c < (power_offset+2)*(power_offset+1)
    s=sum((-1)**i*c**(2*i)/factorial(2*i+power_offset)
          for i in range(terms))
    s2=s+(-1)**terms*c**(2*terms)/factorial(2*terms+power_offset)
    return min(s,s2),max(s,s2)

def sub(x,y):
    return x[0]-y[1],x[1]-y[0]

def dec(x):
    return str(Decimal(x.numerator)/Decimal(x.denominator))

if __name__=="__main__":
    T1,T2=[exact_theta(c,N,R) for c in (c1,c2)]
    B1,B2=[altern(c,1) for c in (c1,c2)]
    D1,D2=sub(T1,B1),sub(T2,B2)
    a,b=altern(c1,0),altern(c2,0)
    K=(4*b[0]/a[1],4*b[1]/a[0])
    margin=(K[0]*D1[0]-D2[1],K[1]*D1[1]-D2[0])
    third=(Q(1361354,10**9),Q(1361355,10**9))
    assert third[0] < margin[0] <= margin[1] < third[1]
    assert Q(1047567,10**9)<D1[0]<=D1[1]<Q(1047568,10**9)
    assert Q(478353,10**9)<D2[0]<=D2[1]<Q(478354,10**9)
    eps=Q(1,2500)
    assert min(D1[0],D2[0])>eps
    assert margin[0]-(1+K[1])*eps>Q(25,10**5)

    # Positive displacement-only control: Y=cos(c1 Z), Y∈[cos c1,1].
    # Exact interval proof of strict feasibility of first/second moments.
    y1lo,y1hi=T1
    y2lo,y2hi=map(lambda x:(1+x)/2,T2)
    assert a[1]<y1lo<y1hi<1
    assert y1hi**2<y2lo
    assert y2hi<(1+a[0])*y1lo-a[1]
    # Choose U=1-Y∈{0,r/e}, mass e²/r at r/e, where e=EU,r=EU².
    # Bound support strictly inside [-1,1] and inversion-symmetric.
    elo,ehi=1-y1hi,1-y1lo
    rlo,rhi=1-2*y1hi+y2lo,1-2*y1lo+y2hi
    assert 0<elo<=ehi
    assert ehi**2<rlo
    assert rhi/elo < 1-a[1]
    # Useful conservative bound |Z|<0.9, certified by monotone cos.
    assert 1-rhi/elo>altern(c1*Q(9,10),0)[1]

    # Source-independent R-only null control, and check its two inequalities tight.
    U1,U2=[sub(exact_theta(c,Q(0),R),altern(c,1)) for c in (c1,c2)]
    assert U1[0]>0 and U2[0]>K[1]*U1[1] and U2[1]<4*U1[0]

    # Low-order exact identities via the Legendre recurrence.
    for c in [Q(0),Q(1,7),Q(4,5),Q(1)]:
        v=[Q(1)]
        for k in range(1,6):
            v=step(v,c,N,R)
            if k==3: assert v[0]==R*c*c/3
        assert v[0]==R**3*c*c/3-(4*N+22*R)*c**4/45
    with localcontext() as ctx:
        ctx.prec=18
        for name, vals in [('Theta1',T1),('Theta2',T2),('D1',D1),('D2',D2),
                           ('K',K),('margin',margin)]:
            print(name,*[dec(x) for x in vals])
    print("PASS exact finite separator and nonzero resistive-only control")
    print("PASS nonstationary-positive displacement law compatible with full finite witness")
    print("PASS 0.0004 each residual error, even keeping both residuals positive")
