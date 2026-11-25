"""Exact finite-network and rational interpolation controls. Python standard library."""
from fractions import Fraction as F
from math import prod


def f(r,w):
    return r*r/(r*r+w*w)


def h(r,w):
    return r*w*w/(r*r+w*w)


def D(r,waves):
    return prod(r*r+w*w for w in waves)


def weights(nodes,waves):
    """Derivative at 0 of rational cardinal interpolants, exact positive nodes."""
    assert len(nodes)==2*len(waves)+1 and all(r>0 for r in nodes)
    assert len(set(nodes))==len(nodes)
    out=[]
    for r in nodes:
        other=[x for x in nodes if x!=r]
        out.append(D(r,waves)/D(F(0),waves)*prod(-x/(r-x) for x in other)
                   *(-sum(1/x for x in other)))
    for k in range(len(nodes)):
        target=1/D(F(0),waves) if k==1 else 0
        assert sum(w*r**k/D(r,waves) for r,w in zip(nodes,out))==target
    return out


def Y(star, lam, c0, w):
    """Direct phasor from actual positive conductance G=c*r for each branch."""
    R,I=lam,w*c0
    for r,c in star:
        G=c*r
        assert G>0 and c>0
        # Given surface phasor 1, reservoir phasor G/(G+i*w*c).
        re_u=G*G/(G*G+(w*c)**2)
        im_u=-G*w*c/(G*G+(w*c)**2)
        R+=G*(1-re_u)
        I-=G*im_u
    return R,I


def two_wave():
    w1,w2=waves=F(1),F(4)
    nodes=tuple(map(F,("1/1000","469/1000","2","8531/1000","1000")))
    ww=weights(nodes,waves)
    assert [w>0 for w in ww]==[False,True,False,True,False]
    C=F(1)
    c0=F(1,2)
    V=sum(abs(w) for w in ww)
    lam_gap=2*C/V
    assert F(199,900)<lam_gap<F(2,9)
    A=[(r,lam_gap*w) for r,w in zip(nodes,ww) if w>0]
    B=[(r,-lam_gap*w) for r,w in zip(nodes,ww) if w<0]
    assert sum(c for r,c in A)==sum(c for r,c in B)==C
    assert sum(c*h(r,w1) for r,c in A)-sum(c*h(r,w1) for r,c in B)==lam_gap
    for w in waves:
        assert Y(A,F(2),c0,w)==Y(B,F(2)+lam_gap,c0,w)
    print("finite connected (1,4) star pair: both complex admittances exactly equal")
    print("restoring gap >199/900, <2/9 (display only: %.12f)" % float(lam_gap))
    for r in map(F,(0,1,2,3,5)):
        RE=(r*r-9*r+4)**2
        IM=18*r*(r-2)**2
        assert RE+IM==D(r,waves)
        g=h(r,w1)+(RE-IM)/(9*D(r,waves))
        assert g==F(1,9)+F(16,15)*(f(r,w1)-f(r,w2))-F(1,5)*(h(r,w2)-h(r,w1))
    # Numerator degree at most 4, so the five rational checks certify identity.
    assert sum(map(abs,(F(4,5),F(1,5),-F(16,15),F(4,15))))==F(7,3)
    print("explicit minimax estimator expansion exact; complex component error factor 7/3")


def one_wave():
    w=F(1)
    c0=F(1,2)
    # Strictly positive finite rate interval [3/4, 4/3].
    lo,hi=F(3,4),F(4,3)
    A=[(F(1),F(1))]
    B=[(lo,F(1,2)),(hi,F(1,2))]
    assert Y(A,F(2),c0,w)==Y(B,F(101,50),c0,w)
    assert 2/sum(abs(x) for x in weights((lo,F(1),hi),(w,)))==F(1,50)
    # Known ceiling 3/4. Fixed one-box exact record, no disconnected physical branch.
    # Formula: kappa=4/b*(1+sqrt(1+b^2))=12; limiting W=1/6.
    r0=F(1,3)
    A=[(r0,F(1))]
    target=f(r0,w)
    assert (4/lo)*(1+F(5,4))==12
    assert target/f(lo,w)==F(5,18)
    limiting_gap=h(r0,w)-(target/f(lo,w))*h(lo,w)
    assert limiting_gap==F(1,6)
    for r in (F(1,10),F(1,100),F(1,1000)):
        high=(target-f(r,w))/(f(lo,w)-f(r,w))
        B=[(r,1-high),(lo,high)]
        gap=h(r0,w)-sum(c*h(rate,w) for rate,c in B)
        assert 0<gap<limiting_gap and Y(A,F(2),c0,w)==Y(B,F(2)+gap,c0,w)
    print("one-tone finite band width exactly 1/50; ceiling width tends to 1/6 on one fixed record")


def cp(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def cq(a,b):
    u,v=cp(a,(b[0],-b[1]))
    d=b[0]**2+b[1]**2
    return u/d,v/d


def nonstar():
    # Two deep capacities =1. Coupling g12=12/5, surface couplings g01=9/5,g02=2/5.
    # Dirichlet interior K=[[21/5,-12/5],[-12/5,14/5]], rates 1,6.
    # Constant-vector mass projected onto rational unit eigenvectors is 49/25,1/25.
    g,h12,k1,k2=F(9,5),F(12,5),F(21,5),F(14,5)
    gg=F(2,5)
    lam,c0=F(3,2),F(7,5)
    for w in map(F,(1,4,7)):
        U=cp((k1,w),(k2,w))
        det=(U[0]-h12*h12,U[1])
        u1=cq((g*k2+h12*gg,g*w),det)
        u2=cq((gg*k1+h12*g,gg*w),det)
        direct=lam+g+gg-g*u1[0]-gg*u2[0],w*c0-g*u1[1]-gg*u2[1]
        assert direct==Y([(F(1),F(49,25)),(F(6),F(1,25))],lam,c0,w)
    print("reciprocally coupled nonstar exact nodal phasor matches spectral star")


if __name__=="__main__":
    two_wave()
    one_wave()
    nonstar()
    print("PASS")
