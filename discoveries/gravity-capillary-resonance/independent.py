#!/usr/bin/env python3
"""Independent integer controls using the single fourth-power numerator."""
from collections import defaultdict
from fractions import Fraction
from math import isqrt, lcm

def factor(n):
    assert n>0
    out={}
    while n>1:
        p=next((d for d in range(2,isqrt(n)+1) if n%d==0),n)
        out[p]=out.get(p,0)+1
        n//=p
    return out

def kernel_parts(rsq, alpha, D):
    A,B=alpha.numerator,alpha.denominator
    N=rsq*D
    assert N.denominator==1 and N>0
    N=int(N)
    b=B*D
    L=b+A*N
    # Compute valuations on the two factors to avoid factoring the entire product.
    f,g=factor(N),factor(L)
    radical=weight=q=rho=1
    for p in f.keys()|g.keys():
        e=f.get(p,0)+2*g.get(p,0)
        radical*=p**(e%4)
        weight*=p**(e//4)
        if e%4 in (1,3):
            q*=p
        if e%4 in (2,3):
            rho*=p
    assert radical*weight**4==N*L*L
    # Recover rather than compute the two squarefree stages.
    assert q*rho*rho==radical and N%q==0
    ell=isqrt(N//q)
    assert q*ell*ell==N
    assert rho*weight**2==ell*(b+A*q*ell*ell)
    assert (rsq*(1+alpha*rsq)**2)==Fraction(radical*weight**4,B**2*D**3)
    c=A*q
    C=c*b*rho*rho
    X=c*rho*ell
    Y=c*rho*rho*weight
    assert X>0 and Y>0 and C>0 and Y*Y==X**3+C*X
    xx,xmod=divmod(X,c*rho)
    yy,ymod=divmod(Y,c*rho*rho)
    assert (xx,yy,xmod,ymod)==(ell,weight,0,0)
    return radical,weight,q,rho,ell

class Lattice:
    def __init__(self, g11, g12, g22):
        a,b,c=(Fraction(z) for z in (g11,g12,g22))
        assert a>0 and a*c>b*b
        self.g=(a,b,c)
        self.D=lcm(*(p.denominator for p in (a,2*b,c)))
    def radius2(self,k):
        a,b,c=self.g
        x,y=k
        return a*x*x+2*b*x*y+c*y*y

def tri(a,b,c):
    return [(1,a),(1,b),(-1,c)]

def inspect(name,lattice,alpha,list_,match,expected=None):
    assert all(sum(sign*k[i] for sign,k in list_)==0 for i in range(2))
    charge=defaultdict(int)
    found=[]
    for sign,k in list_:
        cls,j,q,rho,l=kernel_parts(lattice.radius2(k),alpha,lattice.D)
        charge[cls]+=sign*j
        found.append((q,rho,l,j))
    assert (not any(charge.values()))==match
    if expected:
        assert found==expected
    print("{}: (q,rho,ell,j) {}; sector charges {}".format(name,found,
                                                           sorted(charge.items())))
    return found

def gaussian_counts(maxe=10):
    # Enumerate only the chosen exact shells, not the resonance graph or arbitrary radii.
    for e in range(maxe+1):
        g=5**e
        S=[]
        for x in range(-isqrt(g),isqrt(g)+1):
            v=g-x*x
            y=isqrt(v)
            if y*y==v:
                S.append((x,y))
                if y:
                    S.append((x,-y))
        assert len(S)==len(set(S))==4*(e+1)
    print("square angular shell 5**{}: {} distinct modes in one sector".format(maxe,len(S)))

def main():
    F=Fraction
    square=Lattice(1,0,1)
    x=tri((12,5),(12,-5),(24,0)); alpha=F(7,1259)
    assert alpha*square.radius2(x[0][1])<1<alpha*square.radius2(x[2][1])
    assert x[0][1][0]*x[1][1][1]!=x[0][1][1]*x[1][1][0]
    assert [F(k)*(1+alpha*k*k) for k in (13,24)]==[F(31746,1259),F(126984,1259)]
    inspect("crossover square",square,alpha,x,True,
            [(1,31746,13,1),(1,31746,13,1),(1,31746,24,2)])
    x=tri((1,0),(1,0),(2,0))
    inspect("collinear repeated input",square,F(1,2),x,True,
            [(1,3,1,1),(1,3,1,1),(1,3,2,2)])
    inspect("one-squarefree-label false positive",square,F(1),x,False,
            [(1,2,1,1),(1,2,1,1),(1,10,2,1)])
    rect=Lattice(1,0,F(1,2))
    x=tri((17,12),(17,-12),(34,0))
    inspect("irrational aspect / rational squared aspect",rect,F(7,1978),x,True,
            [(2,171190,19,1),(2,171190,19,1),(2,171190,34,2)])
    # Non-diagonal rational Gram with norm^2=(m+n/2)^2+n^2.
    skew=Lattice(1,F(1,2),F(5,4))
    x=tri((19,10),(29,-10),(48,0))
    assert [skew.radius2(k) for _,k in x]==[26**2,26**2,48**2]
    inspect("non-diagonal Gram control",skew,alpha/4,x,True)
    quartet=[(1,(1,1)),(1,(2,-1)),(-1,(1,-1)),(-1,(2,1))]
    x=inspect("angular two-sector quartic",square,F(1),quartet,True,
              [(2,3,1,1),(5,6,1,1),(2,3,1,1),(5,6,1,1)])
    assert len(set(v[:2] for v in x))==2
    assert all(sum(v[:2] == s for v in x)==2 for s in {v[:2] for v in x})
    t=tri((4,3),(4,-3),(8,0))
    inspect("tuned square",square,F(1),t,True)
    t+= [(1,(1,1)),(-1,(1,1))]
    x=inspect("formal quintic with triad and a distinct conjugate pair",square,F(1),t,True)
    assert len({v[:2] for v in x})==2
    # Rational sign and explicit parameter-dependent enclosure for detuning, no floating point.
    for m in (7,101,100001):
        beta=1+F(1,m)
        t1=5*(1+25*beta); t2=8*(1+64*beta)
        assert 4*t1-t2==-F(12,m)
        # Positive denominator is >4 sqrt(130)>44; square comparison certified.
        assert t1>130 and t2>520 and 130>11**2
        # Therefore detuning lies strictly between -3/(11m) and 0.
        energy=inspect(f"detuned same geometry M={m}",square,beta,t[:3],False)
        assert len({v[:2] for v in energy}) == 2
    gaussian_counts()
    try:
        kernel_parts(F(0),F(1),1)
    except AssertionError:
        pass
    else:
        raise AssertionError("zero mode not allowed")
    print("PASS: independent integer/fraction arithmetic on all fixed controls.")

if __name__=='__main__':
    main()
