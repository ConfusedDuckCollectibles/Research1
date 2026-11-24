#!/usr/bin/env python3
"""Exact interval replay. Direct signed Taylor of exp(tQ), not uniformization."""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from math import factorial, isqrt

@dataclass(frozen=True)
class Box:
    lo: Q
    hi: Q
    def __post_init__(self):
        assert isinstance(self.lo,Q) and isinstance(self.hi,Q) and self.lo <= self.hi
    def __add__(self, other):
        other = interval(other)
        return Box(self.lo+other.lo,self.hi+other.hi)
    __radd__ = __add__
    def __neg__(self):
        return Box(-self.hi,-self.lo)
    def __sub__(self, other):
        return self+(-interval(other))
    def __rsub__(self, other):
        return interval(other)+(-self)
    def __mul__(self, other):
        other=interval(other)
        x=[j*k for j in (self.lo,self.hi) for k in (other.lo,other.hi)]
        return Box(min(x),max(x))
    __rmul__ = __mul__
    def __truediv__(self,other):
        other=interval(other)
        assert not other.lo <= 0 <= other.hi
        return self*Box(1/other.hi,1/other.lo)
    def quantize(self,scale=10**70):
        return Box(Q((self.lo*scale)//1,scale),Q(-((-self.hi*scale)//1),scale))
    def sqrt(self):
        assert self.lo>=0
        scale=10**65
        a=isqrt((self.lo*scale*scale)//1)
        b=isqrt((self.hi*scale*scale)//1)+1
        assert Q(a,scale)**2<=self.lo<=self.hi<=Q(b,scale)**2
        return Box(Q(a,scale),Q(b,scale))
    def log(self):
        x=self.quantize()
        assert 0<x.lo
        return Box(log_bound(x.lo).lo,log_bound(x.hi).hi).quantize()
    def atanh(self):
        assert 0<=self.lo<=self.hi<1
        return ((1+self)/(1-self)).log()/2
    def approx(self):
        x=(self.lo+self.hi)/2
        with localcontext() as ctx:
            ctx.prec=15
            return f"{Decimal(x.numerator)/Decimal(x.denominator):.12g}"

def interval(x):
    return x if isinstance(x,Box) else Box(Q(x),Q(x))

def log_bound(x):
    # Integrate the finite geometric expansion for 2/(1-z^2) from zero to z.
    # The remaining denominator is at least 1-z^2 and z<=1/3.
    power2=0
    while x<1:
        x *= 2
        power2-=1
    while x>2:
        x /= 2
        power2+=1
    def integrate(y,n=78):
        z=(y-1)/(y+1)
        s=Q(0)
        u=z
        for j in range(n):
            s+=2*u/(2*j+1)
            u*=z*z
        return Box(s,s+2*u/((2*n+1)*(1-z*z)))
    return integrate(x)+power2*integrate(Q(2))

def make_ring(forward,backward):
    n=len(forward)
    assert n>=3 and n==len(backward)
    q=[[Q(0)]*n for _ in range(n)]
    for i in range(n):
        j=(i+1)%n
        q[i][j]=Q(forward[i])
        q[j][i]=Q(backward[i])
    for i in range(n):
        q[i][i]=-sum(q[i])
    return q

def transpose(q):
    return [list(row) for row in zip(*q)]

def check_generator(q,pi):
    assert sum(pi)==1 and min(pi)>0
    n=len(q)
    for i in range(n):
        assert sum(q[i])==0
        for j in range(n):
            if i!=j:
                assert q[i][j]>=0 and (q[i][j]==0)==(q[j][i]==0)
    assert all(sum(pi[i]*q[i][j] for i in range(n))==0 for j in range(n))

def check_reflection(q,labels,perm):
    assert sorted(perm)==list(range(len(q)))
    assert all(labels[perm[i]]==labels[i] for i in range(len(q)))
    assert all(q[perm[i]][perm[j]]==q[j][i]
               for i in range(len(q)) for j in range(len(q)))

def Taylor_vec(q,t,f,N=32):
    # Direct signed series; norm majorizes each vector term and the geometric tail.
    z=t*max(sum(abs(qij) for qij in row) for row in q)
    assert z < N+2 and max(abs(x) for x in f)<=1
    term=list(map(Q,f))
    accum=term[:]
    for k in range(1,N+1):
        term=[t*sum(a*b for a,b in zip(row,term))/k for row in q]
        accum=[a+b for a,b in zip(accum,term)]
    tail=z**(N+1)/factorial(N+1)/(1-z/(N+2))
    return [Box(max(Q(0),c-tail),min(Q(1),c+tail)).quantize() for c in accum]

def ep(q,pi):
    # Unordered stationary-current convention only, independent of ring simplification.
    total=interval(0)
    n=len(q)
    for i in range(n):
        for j in range(i+1,n):
            u=pi[i]*q[i][j]
            v=pi[j]*q[j][i]
            if u and u!=v:
                total+=(u-v)*interval(u/v).log()
    return total

def scores(s,r,t,force_zero=False):
    if force_zero:
        assert (s*s-r).lo<=0<=(s*s-r).hi
        return interval(0),interval(0),interval(0)
    change=s-r
    gap=(s*s-r).sqrt()
    ours=2*gap*(gap/change).atanh()/t
    # Below the diagonal CP minorant uses the geometric mean of the two diagonals.
    deficit=2*(change-((1-2*s+r)*r).sqrt())
    x=deficit.sqrt()
    cp=2*x*x.atanh()/t
    pins=8*(s*s-r)/t
    return ours,cp,pins

def probe(q,labels,pi,t,certificate_zero=()):
    check_generator(q,pi)
    reverse=[[pi[j]*q[j][i]/pi[i] for j in range(len(q))] for i in range(len(q))]
    p=[]; parts=[]; S=[]
    for label in (0,1):
        f=[int(z!=label) for z in labels]
        u=Taylor_vec(q,t,f)
        v=Taylor_vec(reverse,t,f)
        mass=sum(pi[i] for i,z in enumerate(labels) if z==label)
        avg_u=sum(pi[i]*u[i] for i,z in enumerate(labels) if z==label)/mass
        avg_v=sum(pi[i]*v[i] for i,z in enumerate(labels) if z==label)/mass
        assert max(avg_u.lo,avg_v.lo)<=min(avg_u.hi,avg_v.hi)
        s=Box(max(avg_u.lo,avg_v.lo),min(avg_u.hi,avg_v.hi))
        r=sum(pi[i]*u[i]*v[i] for i,z in enumerate(labels) if z==label)/mass
        p.append(mass);S.append(s)
        parts.append(scores(s,r,t,label in certificate_zero))
    # Verify stationary label-change flow on the two visible blocks.
    common_flow=p[0]*S[0]-p[1]*S[1]
    assert common_flow.lo<=0<=common_flow.hi
    return (ep(q,pi),)+tuple(sum(p[i]*parts[i][j] for i in range(2)) for j in range(3))

def record(key,record,window=Q(1,10**17)):
    for score in record:
        assert score.hi-score.lo<window,(key,score.approx())
    print(f"{key}: (sigma, three-point, CP-by-block, Pinsker) = "
          f"{', '.join(x.approx() for x in record)}")

def exact_inequality(S,B,C,P):
    assert 0 <= P.lo <= P.hi < B.lo < B.hi < S.lo
    assert C.hi < B.lo

def family(H,k,j):
    w=[H,k,H,k]
    q=make_ring([v+j for v in w],[v-j for v in w])
    check_reflection(q,[0,0,1,1],[1,0,3,2])
    return q

def main():
    p4=[Q(1,4)]*4
    j=Q(19,20)
    q=family(10,1,j)
    fixed=probe(q,[0,0,1,1],p4,Q(1,10000))
    exact_inequality(*fixed)
    S,B,C,P=fixed
    assert Q("3.45679138911")<B.lo<B.hi<Q("3.45679138913")
    assert Q("2.73865902408")<C.lo<C.hi<Q("2.73865902410")
    assert B.lo>5*C.hi/4
    record("rare finite-lag four-state ring",fixed)
    moderate=probe(family(10,1,Q(1,2)),[0,0,1,1],p4,Q(1,1000))
    exact_inequality(*moderate)
    record("moderate finite-lag four-state ring",moderate)
    db=probe(family(2,1,0),[0,0,1,1],p4,Q(1,8),(0,1))
    assert all(x==interval(0) for x in db)
    record("detailed balance",db)
    # Unequal midpoint probabilities: 0 block has two states, 1 is a singleton.
    q=make_ring([2]*3,[1]*3)
    check_reflection(q,[0,0,1],[1,0,2])
    unequal=probe(q,[0,0,1],[Q(1,3)]*3,Q(1,100),(1,))
    exact_inequality(*unequal)
    record("unequal midpoint mass, singleton block exactly zero",unequal)
    # Independent hidden driven orbit times a reversible two-level process.
    Q8=[[Q(0)]*8 for _ in range(8)]
    ring=family(2,1,Q(1,2))
    for a in range(2):
        for i in range(4):
            k=4*a+i
            for j in range(4):
                Q8[k][4*a+j]=ring[i][j]-int(i==j)
            Q8[k][4*(1-a)+i]=1
    blind=probe(Q8,[0]*4+[1]*4,[Q(1,8)]*8,Q(1,10),(0,1))
    assert blind[0].lo>0 and all(x==interval(0) for x in blind[1:])
    record("positive hidden EP, reversible Markov readout",blind)
    ratios=[]
    for H in (10,100,1000):
        x=probe(family(H,1,Q(1,2)),[0,0,1,1],p4,Q(1,H**3))
        exact_inequality(*x)
        ratios.append(x[1]/x[0])
    assert ratios[0].hi<ratios[1].lo<ratios[1].hi<ratios[2].lo
    assert ratios[-1].lo>Q(997,1000)
    print("heterogeneous sequence rigid finite enclosures for bound/sigma:",
          ", ".join(x.approx() for x in ratios))
    # A positive covariance relaxed pair is exactly modeled by equal latent rates.
    s,r=Q(1,2),Q(17,50)
    w,a=s*s/r,r/s
    assert 0<w<1 and 0<a<1 and w*a==s and w*a*a==r
    # Positive anticovariance input whose equal-weight swapped atoms are rational.
    u,v=Q(7,50),Q(3,50)
    s,r=(u+v)/2,u*v
    B,C,P=scores(interval(s),interval(r),Q(1))
    separate=(u-v)*interval(u*(1-v)/(v*(1-u))).log()/2
    assert max(B.lo,separate.lo)<=min(B.hi,separate.hi) and B.lo>max(C.hi,P.hi)
    print("exact swapped-atom product moments verified, nonlinear cost:",B.approx())
    print("PASS: probability and transcendental output brackets use exact rational arithmetic.")

if __name__=='__main__':
    main()
