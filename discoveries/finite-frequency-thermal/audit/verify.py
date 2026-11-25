"""Independent exact finite controls; Python standard library only."""
from fractions import Fraction as Q
from math import prod,isqrt

def derivatives(freq,nodes):
    """R'(0)'s interpolation weights by differentiating Lagrange polynomials."""
    n=2*len(freq)
    assert len(nodes)==n+1 and all(0<x for x in nodes)
    assert sorted(set(nodes))==nodes
    den=lambda x: prod(x*x+w*w for w in freq)
    return [den(x)/den(0)*sum(prod(-z for z in nodes if z not in (x,y))
                               for y in nodes if y!=x)/
            prod(x-y for y in nodes if y!=x) for x in nodes]

def splits(freq,nodes):
    W=derivatives(freq,nodes)
    assert sum(W)==0 and all((-1)**(i+1)*w>0 for i,w in enumerate(W))
    # Independent sign tests on the physical measurement basis.
    for w in freq:
        assert sum(q/(x*x+w*w) for x,q in zip(nodes,W))==0
        assert sum(q*w*w*x/(x*x+w*w) for x,q in zip(nodes,W))==1
    d=2/sum(map(abs,W))
    plus=[(x,q*d) for x,q in zip(nodes,W) if q>0]
    minus=[(x,-q*d) for x,q in zip(nodes,W) if q<0]
    assert sum(q for x,q in plus)==sum(q for x,q in minus)==1
    # Direct physical one-hidden-box Schur fractions, summed as star branches.
    def Y(w,M,lam):
        r,i=lam,2*w
        for x,c in M:
            k=c*x
            a,b=k,c*w
            r+=k-k*k*a/(a*a+b*b)
            i+=k*k*b/(a*a+b*b)
        return r,i
    assert all(Y(w,plus,Q(1))==Y(w,minus,1+d) for w in freq)
    return d,plus,minus

f1=[Q(1)]
for N in (10,100,1000):
    d,*_=splits(f1,[Q(1,N),Q(1),Q(N)])
    assert d==Q(1,2)-Q(N,1+N*N)
    assert 0<Q(1,2)-d<Q(1,N)
print("m=1 whole-axis: exact matching; width -> 1/2")

freq=[Q(1),Q(4)]
def ABC(x):
    A=x*x-9*x+4
    B=3*x-6
    P=(x*x+1)*(x*x+16)
    assert A*A+2*x*B*B==P
    # Independent partial fraction coefficients for affine estimator.
    r=(A*A-2*x*B*B)/P
    assert r==1-Q(48,5)/(x*x+1)+Q(768,5)/(x*x+16)-Q(36,5)*x/(x*x+1)-Q(9,5)*16*x/(x*x+16)
    return r
for N in (100,10000,1000000):
    sq=Q(isqrt(65*N*N),N)
    nodes=[Q(1,N),(9-sq)/2,Q(2),(9+sq)/2,Q(N)]
    d,*_=splits(freq,nodes)
    assert Q(2,9)-Q(2,N)<d<Q(2,9)
    assert all(-1<=ABC(x)<=1 for x in nodes)
print("m=2 frequencies 1,4: exact matching; width -> 2/9; expansion verifies slope -9")

# omega=1 and b=3/4: phase's middle knot is exactly 1/3; norm D=12.
for N in (100,10000):
    d,*_=splits(f1,[Q(1,N),Q(1,3),Q(3,4)])
    assert Q(1,6)-Q(2,N)<d<Q(1,6)
print("m=1 ceiling 3/4: exact matching; width -> 1/6")

# omega=1, [a,b]=[3/4,4/3]: middle knot exactly 1 and both endpoints physical.
d,*_=splits(f1,[Q(3,4),Q(1),Q(4,3)])
assert d==Q(1,50),d
print("m=1 window [3/4,4/3]: exact matching; width 1/50 attained")

# Generic complex rational control with a positive nonstar conductance network.
class G:
    def __init__(self,a=0,b=0): self.a,self.b=Q(a),Q(b)
    def __add__(x,y):
        y=y if isinstance(y,G) else G(y)
        return G(x.a+y.a,x.b+y.b)
    __radd__=__add__
    def __neg__(x): return G(-x.a,-x.b)
    def __sub__(x,y): return x+-y
    def __mul__(x,y):
        y=y if isinstance(y,G) else G(y)
        return G(x.a*y.a-x.b*y.b,x.a*y.b+x.b*y.a)
    __rmul__=__mul__
    def __truediv__(x,y):
        y=y if isinstance(y,G) else G(y)
        return x*G(y.a,-y.b)*G(1/(y.a*y.a+y.b*y.b))
    def __eq__(x,y): return (x.a,x.b)==(y.a,y.b)
def invform(M,u):
    a,b,c=M
    return (c*(u[0]**2)-b*(2*u[0]*u[1])+a*(u[1]**2))/(a*c-b*b)

# surface links 7,1 and hidden link 4; capacities 2,3,5; lambda=2.
for w in (Q(1),Q(4),Q(19,7)):
    s=G(0,w)
    K=(G(11)+3*s,G(-4),G(5)+5*s)
    Y=G(10)+2*s-invform(K,[7,1])
    assert Y==G(2)+10*s-s*s*invform(K,[3,5])
print("Exact nonstar Schur control: total hidden capacity recovered before diagonalization")
