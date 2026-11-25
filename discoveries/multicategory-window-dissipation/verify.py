#!/usr/bin/env python3
"""Rational certificates for a fixed finite-window CTMC example."""
from fractions import Fraction as F
from math import isqrt
from decimal import Decimal, localcontext

def I(x): return F(x),F(x)
def add(a,b): return a[0]+b[0],a[1]+b[1]
def sub(a,b): return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    x=[i*j for i in a for j in b]
    return min(x),max(x)
def div(a,b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))
def root(a):
    assert a[0]>=0
    D=10**28
    def floor(x): return isqrt(x.numerator*D*D//x.denominator)
    return F(floor(a[0]),D),F(floor(a[1])+1,D)
def artlow(t,n=36):
    assert 0<=t<=F(1,3)
    power=t; s=F(0)
    for i in range(n):
        s+=power/(2*i+1)
        power*=t*t
    return s,s+power/((2*n+1)*(1-t*t))
def rlog(t):
    assert t>0
    k=0
    while t>=2:t/=2;k+=1
    while t<1:t*=2;k-=1
    return mul(I(2),add(artlow((t-1)/(t+1)),mul(I(k),artlow(F(1,3)))))
def log(a):
    D=10**28
    lo=F(int(a[0]*D),D); hi=F(int(a[1]*D)+1,D)
    return rlog(lo)[0],rlog(hi)[1]
def art(a):
    assert 0<=a[0]<1 and a[1]<1
    return mul(I(F(1,2)),log(div(add(I(1),a),sub(I(1),a))))
def phi(a):
    y=root(sub(I(1),mul(a,a)))
    return mul(y,art(y))
def fmt(a):
    with localcontext() as ctx:
        ctx.prec=13; x=sum(a)/2
        return str(Decimal(x.numerator)/Decimal(x.denominator))
def pivot_psd(A):
    """Check by exact elimination: all pivots positive or zero with a zero row."""
    A=[r[:] for r in A]; p=[]
    for k in range(len(A)):
        x=A[k][k]; assert x>=0; p.append(x)
        if x==0:
            assert all(b==0 for b in A[k][k:])
        else:
            for i in range(k+1,len(A)):
                for j in range(k+1,len(A)):
                    A[i][j]-=A[i][k]*A[k][j]/x
    return p

def simulate_algebra(e,z,a):
    # R_e = average of two rank-one matrices with offdiagonal m±sqrt(m²-de df).
    # u, v are their row, column marginals. Verify enough identities to
    # ensure stationary reverse generator and exact classification table.
    scale=(e,F(1),e,F(1),e); theta=a*e/2
    H=[(i,i+1,scale[i]**2,scale[i+1]**2,e) for i in range(4)]
    H[0]=(0,1,e**2*(1-theta),F(1),e)
    H[3]=(3,4,F(1),e**2*(1-theta),e)
    if a:H.append((4,0,e*e*theta,e*e*theta,a*e*e))
    P=[[F(0) for _ in range(5)] for _ in range(5)]
    for i,j,di,dj,m in H:
        assert m*m>=di*dj>0
        disc=m*m-di*dj
        # Both algebraic offdiagonals m±sqrt(disc) are strictly positive.
        assert F(0)<=disc<m*m
        assert m*m-disc==di*dj     # determinant identity in Q[sqrt(disc)]
        P[i][i]+=di; P[j][j]+=dj; P[i][j]+=m; P[j][i]+=m
        # Quadratic forms of the two product marginals average to R_e.
        T=di+dj+2*m
        assert (di+m)*(dj+m)+disc==T*m
        assert (di+m)**2-disc==T*di
        assert (dj+m)**2-disc==T*dj
    Z=sum(map(sum,P)); assert Z==4+8*e+(4+2*a)*e*e
    p=1-z
    T=[[z*z]+[z*p*sum(row)/Z for row in P]]
    T.extend([[z*p*sum(row)/Z]+[p*p*x/Z for x in row] for row in P])
    assert sum(map(sum,T))==1
    pivot_psd(T)
    for mask in range(64):
        S=[i for i in range(6) if (mask>>i)&1]
        r=sum(T[i][j] for i in S for j in S)
        s=sum(sum(T[i]) for i in S)
        assert r>=s*s
    # CP minorant with censor row: no loss wherever not forced by its
    # dependency on the deleted rare edge.
    P[0][4]=P[4][0]=F(0)
    G=[[z*z]+[z*p*sum(row)/Z for row in P]]
    G.extend([[z*p*sum(row)/Z]+[p*p*x/Z for x in row] for row in P])
    assert all(0<=G[i][j]<=T[i][j] for i in range(6) for j in range(6))
    deficit=1-sum(map(sum,G))
    alpha=2*a*e*e/Z
    assert deficit==(1-z*z)*alpha
    return alpha,Z

def calculate():
    e,z,a=F(1,32),F(1,16),F(1,4)
    alpha,Z=simulate_algebra(e,z,a)
    assert simulate_algebra(e,z,F(0))[0]==0
    theta=a*e/2; p=1-z
    t=mul(I(4),log(I(2)))
    bound=mul(I(p*p*alpha),add(I(1),log(I(2/e))))
    d=(1-z*z)*alpha
    r=root(I(d)); old=mul(I(2),mul(r,art(r)))
    ratio=div(bound,old)
    E0=add(mul(I(4*e/Z),phi(root(I(1-theta)))),
           mul(I(alpha),phi(I(e/2))))
    # p_H=1/2, lambda=1, hence actual sigma=E0.
    rate=div(mul(I(F(1,2)),bound),t)
    oldrate=div(mul(I(F(1,2)),old),t)
    recovered=div(rate,E0)
    excursion=mul(I(alpha),add(I(1),log(I(2/e))))
    exc_recovered=div(excursion,E0)

    mu=a*e/100
    eta=(1-z)**2/Z*e*e*mu
    noisy=mul(I(2*(1-z)**2/Z*e*e*(a-mu)),
              sub(add(I(1),log(I(2/e))),I(23*mu/(e*(a-mu)))))
    noisy_ratio=div(noisy,old)
    assert ratio[0]>F(9,4)
    assert noisy_ratio[0]>F(21,10)
    assert F(13,100)<recovered[0]<recovered[1]<F(15,100)
    assert rate[0]>oldrate[1]
    for label,val in [('half-window tau',t),('conditional uncensored rare pair mass',I(alpha)),
                      ('new EPR rate lower bound',rate),('Petkov Eq11 rate ceiling',oldrate),
                      ('guaranteed improvement factor',ratio),('actual EPR rate',E0),
                      ('true rate fraction certified',recovered),
                      ('richer uncensored excursion EPR lower bound',excursion),
                      ('true rate fraction with richer excursion rate data',exc_recovered),
                      ('matched noise L_inf allowance',I(eta)),
                      ('noisy guaranteed improvement factor',noisy_ratio)]:
        assert val[1]-val[0]<F(1,10**16)
        print(label+':',fmt(val))
    print('PASS: exact finite table, 64 blind binary partitions, equilibrium control, and interval gates.')

if __name__=='__main__':calculate()
