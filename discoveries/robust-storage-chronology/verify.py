"""Offline rational checks for exact continuum coarse-bin reserve certificates."""
from fractions import Fraction as Q
from itertools import combinations, product


def arrangement(bins):
    T=bins[-1][1]
    X={0,T}
    Y={0,T,2*T}
    D={0,T}
    for a,b,h in bins:
        ell=b-a-h
        assert 0<=h<=b-a
        X.update([a,b,b-ell])
        Y.update([a,b,a+ell,a+T,b+T,a+T+ell])
        D.update([ell,T-(b-a)+ell])
    V=set(product(X,Y))
    V.update((x,x+d) for x,d in product(X,D))
    V.update((y-d,y) for y,d in product(Y,D))
    return {(t,u) for t,u in V if 0<=t<=T and 0<=u-t<=T}


def w(t,a,b,base):
    # If base=e^rho>1, this is rho times the weighted integral over [a,b].
    # For rho=0 (base=1) return the ordinary integral.
    return b-a if base==1 else base**(b-t)-base**(a-t)


def rearranged(t,u,bins,base,c,d):
    T=bins[-1][1]
    selected=Q(0)
    for a,b,h in bins:
        ell=b-a-h
        pieces=[(max(t,a+shift),min(u,b+shift)) for shift in (0,T)]
        for p,q in reversed(pieces):
            length=min(ell,max(q-p,0))
            if length:
                selected+=w(t,q-length,q,base)
                ell-=length
    return c*w(t,t,u,base)-(c+d)*selected


def patterns(bins,c,d):
    cases=[]
    for a,b,h in bins:
        cases.append([tuple(c if i in inds else -d for i in range(b-a))
                      for inds in combinations(range(b-a),h)])
    return [sum(parts,()) for parts in product(*cases)]


def brute_arc(t,u,F,base):
    mult=base-1 if base>1 else 1
    T=len(F)
    return mult*sum(base**(s-t)*F[s%T] for s in range(t,u))


def simulate(F,phase,Z,base,steps):
    # Independent forward reflected ODE over each constant unit-input step.
    # State is rho*X for positive rho and X for zero rho.
    value=Z
    for i in range(steps):
        f=F[(phase+i)%len(F)]
        potential=(value+(base-1)*f)/base if base>1 else value+f
        value=min(Z,potential)
        if value<0:
            return False
    return True


def check_case(bins,base,c,d):
    T=bins[-1][1]
    V=arrangement(bins)
    assert all(isinstance(t,int) and isinstance(u,int) for t,u in V)
    allF=patterns(bins,c,d)
    cut=lambda t,u:rearranged(t,u,bins,base,c,d)
    for t,u in V:
        assert cut(t,u)==min(brute_arc(t,u,F,base) for F in allF)
    cyc=min(cut(t,u) for t,u in V if u-t==T)
    reserve=max(-cut(t,u) for t,u in V)
    assert cyc==min(brute_arc(t,t+T,F,base) for t in range(T) for F in allF)
    assert reserve==max(-brute_arc(t,u,F,base) for F in allF
                        for t in range(T) for u in range(t,t+T+1))
    if cyc>=0:
        assert all(simulate(F,t,reserve,base,3*T) for F in allF for t in range(T))
        if reserve>0:
            assert not all(simulate(F,t,reserve-Q(1,10**8),base,T)
                           for F in allF for t in range(T))
    return V,cyc,reserve


def gross_meter():
    T=8
    base=Q(6,5)
    c,d=Q(1),Q(1)
    assert 2*base**5-base**8-1<0  # net-only coarse contract
    allF=[]
    for pos in combinations(range(T),3):
        for neg in set(range(T))-set(pos):
            allF.append(tuple(c if i in pos else -d if i==neg else Q(0)
                              for i in range(T)))
    V=arrangement(((0,T,3),))|arrangement(((0,T,7),))
    # Form cross intersections of all breakpoint lines: all are integer in
    # this control, so checking all integer cuts covers every vertex.
    V.update((t,u) for t in range(T+1) for u in range(t,t+T+1))
    def gross(t,u):
        M_pos=rearranged(t,u,((0,T,3),),base,Q(1),Q(0))
        M_neg=rearranged(t,u,((0,T,7),),base,Q(0),Q(1))
        return M_pos+M_neg
    assert all(gross(t,u)==min(brute_arc(t,u,F,base) for F in allF) for t,u in V)
    cyc=min(gross(t,u) for t,u in V if u-t==T)
    assert cyc==base**3-1-base**8+base**7>0
    res=max(-gross(t,u) for t,u in V)
    assert res==base**5-base**4==Q(1296,3125)
    assert all(simulate(F,0,res,base,3*T) for F in allF)
    print("separate gross records rescue net-only infinite case; exact rho*E=1296/3125")


def main():
    # Two bins, positive coarse mean throughout. Worst full-cycle cut starts inside.
    bins=((0,3,2),(3,6,2))
    base=Q(11,10)
    V,cyc,res=check_case(bins,base,Q(1),Q(1))
    assert cyc==Q(156639,10**6) and res==Q(21,100)
    assert min(rearranged(t,t+6,bins,base,Q(1),Q(1)) for t in (0,3))==Q(207459,10**6)
    assert (1,7) in V and rearranged(1,7,bins,base,Q(1),Q(1))==cyc
    assert brute_arc(2,4,(1,1,-1,-1,1,1),base)==-res
    print("two-bin exact all-history reserve rho*E=21/100; worst full-cycle is interior")

    # A bin-boundary-only robust feasibility test can even certify the wrong sign.
    _,cyc,_=check_case(bins,base,Q(3,5),Q(1))
    boundary=min(rearranged(t,t+6,bins,base,Q(3,5),Q(1)) for t in (0,3))
    assert boundary==Q(11655,10**6)>0 and cyc==Q(-29001,10**6)<0
    print("boundary-only cycle check false-positive: +11655/10^6 vs interior -29001/10^6")

    # Equality at leaky feasibility threshold with asymmetric rates and physical caps.
    # e^rho=2, h=1/2 period, endpoint rates +2,-1.
    V,cyc,res=check_case(((0,2,1),),Q(2),Q(2),Q(1))
    assert cyc==0 and res==1
    print("one-bin equality case: weighted full-cycle=0 yet exact reserve rho*E=1")

    # Converter clipping can coexist with coarse raw-net totals; check an alternate instance.
    # P+=2, P_c=1, eta_c=3/4 => c=3/4; P-=1, eta_d=4/5 => d=5/4.
    bins2=((0,3,2),(3,7,3))
    V,cyc,res=check_case(bins2,Q(21,20),Q(3,4),Q(5,4))
    assert cyc>0 and res>0
    # Lossless boundary of the same theorem.
    V,cyc,res=check_case(bins,Q(1),Q(1),Q(1))
    assert cyc==2 and res==2
    print("asymmetric clipped-conversion case and lossless limits cross-check")

    # Explicit repeating trace: averaged profile p=1/100, E=0 is feasible;
    # positive 101 then negative 99 has negative no-spill periodic charge at start.
    base=Q(1001,1000)
    total=2*base**101-base**200-1
    assert total<-Q(1,1000)
    assert Q(101-99,200)==Q(1,100)>0
    assert total/(base**200-1)<0
    print("zero-capacity coarse prediction vs no finite feasible reserve at 1% mean surplus")
    gross_meter()
    print("PASS")


if __name__=="__main__":
    main()
