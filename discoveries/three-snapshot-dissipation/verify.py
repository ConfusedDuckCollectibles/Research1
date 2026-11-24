#!/usr/bin/env python3
"""Standard-library rational enclosures, not simulation.

Probability enclosures use uniformization & Poisson tail. Log enclosures
use atanh series after reducing by powers of 2. Decimals print midpoints;
all PASS conditions are asserted against rational interval endpoints.
"""
from fractions import Fraction as F
from math import isqrt
from decimal import Decimal, localcontext

def I(x):
    return (F(x), F(x))

def plus(a, b):
    return (a[0] + b[0], a[1] + b[1])

def neg(a):
    return (-a[1], -a[0])

def sub(a, b):
    return plus(a, neg(b))

def mul(a,b):
    c = [x*y for x in a for y in b]
    return (min(c), max(c))

def div(a,b):
    assert not b[0] <= 0 <= b[1]
    return mul(a,(1/b[1],1/b[0]))

def comp(a):
    return sub(I(1),a)

def sqrt_f(x):
    assert x>=0
    D=10**34
    y=isqrt((x.numerator*D*D)//x.denominator)
    return F(y,D), F(y+1,D)

def sqrt_i(x):
    return sqrt_f(x[0])[0],sqrt_f(x[1])[1]

def expneg(z,n=28):
    assert 0<=z<=1
    term=tot=F(1)
    for i in range(1,n+1):
        term *= -z/i
        tot += term
    following = tot+term*(-z/(n+1))
    return min(tot,following),max(tot,following)

def atanh_pos(z,n=40):
    assert 0<=z<=F(1,3)
    power=z; acc=F(0)
    for i in range(n):
        acc += power/(2*i+1)
        power *= z*z
    return acc,acc+power/((2*n+1)*(1-z*z))

def ln_pos(x):
    D=10**34
    lower=F(int(x*D),D)
    upper=lower+F(1,D)
    assert 0<lower<=x<=upper
    return _ln(lower)[0],_ln(upper)[1]

def _ln(x):
    j=0
    while x>=2:
        x/=2
        j+=1
    while x<1:
        x*=2
        j-=1
    return mul(I(2),plus(atanh_pos((x-1)/(x+1)),
                         mul(I(j),atanh_pos(F(1,3)))))

def ln(x):
    return ln_pos(x[0])[0],ln_pos(x[1])[1]

def atanh(x):
    assert 0<=x[0] and x[1]<1
    return mul(I(F(1,2)),ln(div(plus(I(1),x),comp(x))))

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def mat_column(M,x):
    return [dot(row,x) for row in M]

def kernel(H,k,j):
    rates=[H,k,H,k]
    Q=[[F(0) for _ in rates] for _ in rates]
    for i,t in enumerate(rates):
        p=(i+1)%4
        Q[i][p]=t+j
        Q[p][i]=t-j
    for i in range(4):
        Q[i][i]=-(H+k)
    return Q

def check_structure(Q):
    z=[0,0,1,1]
    reflect=lambda i:(1-i)%4
    rot=lambda i:(i+2)%4
    assert all(sum(row)==0 and all(v>0 for n,v in enumerate(row)
                                  if n!=i and v!=0)
               for i,row in enumerate(Q))
    assert all(sum(Q[i][j] for i in range(4))==0 for j in range(4))
    for i in range(4):
        assert z[reflect(i)]==z[i] and z[rot(i)]==1-z[i]
        for j in range(4):
            assert Q[reflect(i)][reflect(j)]==Q[j][i]
            assert Q[rot(i)][rot(j)]==Q[i][j]
            assert (Q[i][j]!=0)==(Q[j][i]!=0)

def p_forward(H,k,j,t,N=28):
    Q=kernel(H,k,j); check_structure(Q)
    lam=H+k
    z=lam*t
    assert 0<z<=1
    M=[[F(int(i==j))+q/lam for j,q in enumerate(row)] for i,row in enumerate(Q)]
    assert all(sum(row)==1 and min(row)>=0 for row in M)
    f=list(map(F,[0,0,1,1]))
    c=F(1); cumulative=[c*x for x in f]
    for n in range(1,N+1):
        f=mat_column(M,f)
        c*=z/n
        cumulative=[a+c*b for a,b in zip(cumulative,f)]
    # Every M^n f coordinate lies in [0,1].
    rem=c*z/(N+1)/(1-z/(N+2))
    E=expneg(z)
    return [(E[0]*p,E[1]*(p+rem)) for p in cumulative[:2]]

def scale(a,b):
    return mul(a,I(b))

def ep(H,k,j):
    return scale(plus(ln(I((k+j)/(k-j))),ln(I((H+j)/(H-j)))),j)

def witness_from_uv(u,v,t):
    hdiff=sub(u,v)
    s=scale(plus(u,v),F(1,2)); r=mul(u,v)
    d2=scale(hdiff,F(1,2))
    g=mul(d2,d2)
    if hdiff[0]>0:
        w=sub(plus(u,v),scale(r,2))
        ours=scale(mul(hdiff,atanh(div(hdiff,w))),1/t)
        cpdelta=sub(w,scale(sqrt_i(mul(r,mul(comp(u),comp(v)))),2))
        x=sqrt_i(cpdelta)
        cp=scale(mul(x,atanh(x)),2/t)
    else:
        assert hdiff[0]<=0<=hdiff[1],hdiff
        ours=cp=None
    linear=scale(g,8/t)
    return s,r,ours,cp,linear

def printable(x,n=12):
    if x is None:return "0 (exact symmetry)"
    a,b=x
    with localcontext() as ctx:
        ctx.prec=n
        q=(a+b)/2
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def confirm_product_controls():
    # s=1/10, r=21/2500 < s^2, roots 7/50 and 3/50
    u,v=F(7,50),F(3,50)
    s=(u+v)/2; r=u*v
    # For equality controls use exact arithmetic on roots, not floating math.
    assert s==F(1,10) and r==F(21,2500) and s*s-r==F(1,625)
    assert u*(1-v)-v*(1-u)==u-v
    assert (u*(1-v))*(v*(1-u))==r*(1-u)*(1-v)
    _,_,ours,cp,linear=witness_from_uv(I(u),I(v),F(1))
    jhalf=scale(ln(I((u*(1-v))/(v*(1-u)))),(u-v)/2)
    assert max(ours[0],jhalf[0])<min(ours[1],jhalf[1])
    assert ours[0]>cp[1] and ours[0]>linear[1]
    print('exact two-class product, tight certificate:',printable(ours))
    # Zero cost at positive covariance via identical conditionals.
    s,r=F(1,2),F(17,50)
    ps=r/s; weighting=s*s/r
    assert s*s<r<=s and weighting*ps==s and weighting*ps*ps==r

def main():
    confirm_product_controls()
    cases=[(F(2),F(1),F(0),F(1,8)),
           (F(10),F(1),F(1,2),F(1,1000)),
           (F(10),F(1),F(19,20),F(1,10000))]
    for H,k,j,t in cases:
        v,u=p_forward(H,k,j,t)  # state 0 has backward visible edge (k-j)
        s,r,ours,cp,linear=witness_from_uv(u,v,t)
        e=ep(H,k,j)
        if j:
            assert ours[1]<e[0] and ours[0]>cp[1] and ours[0]>linear[1]
            assert ours[1]-ours[0]<F(1,10**15)
            if j==F(19,20):
                assert ours[0]>F(5,4)*cp[1] and ours[0]>F(100)*linear[1]
        else:
            assert e[0]<=0<=e[1] and ours is None
        print("H,k,j,t:",H,k,j,t)
        print(" s,r:",printable(s),printable(r))
        print(" sigma:",printable(e) if j else "0 (exact symmetry)")
        print(" three-reading floor:",printable(ours),
              "binary CP floor:",printable(cp),
              "same-triple Pinsker floor:",printable(linear) if j else "0 (exact symmetry)")
        if j: print(" recovered fraction:",printable(div(ours,e)))
    print('PASS: all asserted inequalities are certified by rational enclosures.')

if __name__=='__main__':
    main()
