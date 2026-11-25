#!/usr/bin/env python3
"""Exact independent probability/phase-label replay; no optimizer or author-code import."""
from fractions import Fraction as Q
from itertools import combinations, product

X=((0,0),(1,0),(0,1),(1,1))

def as_p(num,den):
    z=tuple(Q(n,den) for n in num)
    assert sum(z)==1 and min(z)>=0
    return z

def step(x,i,s):
    # output bit, physical Kraus label, conditional mass.
    if not x[i]:
        return [(0,0,Q(1))]
    return [(1,0,s),(0,1,1-s)]

def sector_outputs(p,s,order=(0,1),reset=None):
    # We propagate squared amplitudes because each joint tagged-basis
    # output has a unique initial x within its Kraus sector. Accumulated
    # integer phase labels are attached to the remaining sensor bits.
    i,j=order
    sectors={L:[] for L in X}
    jumped=Q(0)
    for x,w in zip(X,p):
        for u,L1,A in step(x,i,s[i]):
            W=w*A
            assert u^x[i]==L1
            if reset is not None and L1:
                jumped+=W
                continue
            for v,L2,B in step(x,j,s[j]):
                assert v^x[j]==L2
                sectors[(L1,L2)].append(((u,v,x[i],x[j]),u+v,W*B))
    if reset is not None:
        for b,fb in ((0,1-reset),(1,reset)):
            for v,L2,B in step((b,),0,s[j]):
                assert v^b==L2
                sectors[(1,L2)].append(((0,v,1,b),v,jumped*fb*B))
    seen=set()
    for L,data in sectors.items():
        for y,k,W in data:
            assert W>=0 and k==y[0]+y[1]
            if W:
                assert (y[0]^y[2],y[1]^y[3])==L and y not in seen
                seen.add(y)
    assert sum(W for data in sectors.values() for _,_,W in data)==1
    return sectors

def qfi(sectors):
    # For a local optimal +/- measurement in each pure sector:
    # p_plus=p_minus=T/2, slope magnitudes T sqrt(variance).
    # The unnormalized variance numerator is the two-copy pair sum.
    val=Q(0)
    for row in sectors.values():
        T=sum(w for _,_,w in row)
        if T:
            Z=sum(w*v*(k-l)**2 for (_,k,w),(_,l,v) in combinations(row,2))
            val+=4*Z/T
    return val

def upper(s,means):
    # Trial posterior loss means, computed by the generic Kraus path,
    # retaining excitation count on input for branch differences.
    scores=[]
    for x in X:
        c=Q(0)
        for _,a,A in step(x,0,s[0]):
            for _,b,B in step(x,1,s[1]):
                c+=A*B*(sum(x)-means[(a,b)])**2
        scores.append(c)
    return 4*max(scores)

def U(a,b,c):
    return {(0,0):a,(1,0):b,(0,1):c,(1,1):2}

def ghz_control():
    s=(Q(3,5),Q(3,5)); p=as_p((3,0,0,5),8)
    assert qfi(sector_outputs(p,s))==upper(s,U(Q(3,4),2,2))==Q(9,4)
    assert qfi(sector_outputs(p,s,reset=Q(9,16)))==Q(2421,992)
    # Independent countercontrol: no tags for the same GHZ state.
    # Phase-carrying 00/11 block includes double-loss mass in its diagonal.
    a,b=s; r,_,_,t=p; c2=r*a*b*t
    assert 16*c2/(r+(1-a)*(1-b)*t+a*b*t)==Q(27,14)<Q(9,4)
    s=(Q(16,25),Q(81,100)); p=as_p((18,0,0,25),43)
    P=qfi(sector_outputs(p,s))
    assert P==upper(s,U(Q(36,43),2,2))==Q(5184,1849)
    A=qfi(sector_outputs(p,s,reset=Q(10,19)))
    B=qfi(sector_outputs(p,s,(1,0),reset=Q(5,9)))
    assert A==Q(1996812,667489) and A-P==Q(2916,15523)
    assert B-P==Q(304,3483) and A>B
    print("Matched diagonal and unequal GHZ exact controls; untagged diagonal GHZ loses QFI (27/14 < 9/4).")

def high_loss_controls():
    s=(Q(1,4),Q(1,4)); p=as_p((4,3,3,12),22)
    assert upper(s,U(Q(491,1000),Q(743,500),Q(743,500)))==Q(967669,1000000)
    assert qfi(sector_outputs(p,s,reset=Q(2,3)))==Q(1087,1100)
    s=(Q(4,25),Q(9,25)); p=as_p((33,41,19,107),200)
    P=qfi(sector_outputs(p,s))
    assert P==qfi(sector_outputs(p,s,(1,0)))==Q(1173318793,1224558300)
    A=qfi(sector_outputs(p,s,reset=Q(5,8)))
    B=qfi(sector_outputs(p,s,(1,0),reset=Q(5,7)))
    cap=upper(s,U(Q(489,1000),Q(1485,1000),Q(1477,1000)))
    assert cap==Q(149825737,156250000)
    assert A==Q(3360283141,3449460000)<B==Q(9944872779,10204652500)
    assert B-A==Q(101407,256452000)
    assert B-cap==Q(9986716030943,637790781250000)>Q(1565,100000)
    assert A>cap and P<=cap
    print("Unequal high loss: cap",cap,"noisy first",A,"less noisy first",B,"order difference",B-A,"certified gap",B-cap)

def boundary_controls():
    s=(Q(9,16),Q(1)); p=as_p((3,0,0,4),7)
    assert qfi(sector_outputs(p,s))==upper(s,U(Q(6,7),2,2))==Q(144,49)
    assert qfi(sector_outputs(p,s,reset=Q(1,2)))==Q(625,196)
    assert qfi(sector_outputs(p,s,(1,0),reset=Q(4,7)))==Q(144,49)
    # Replacer before or after the only useful call: choose its input ground.
    s=(Q(9,16),Q(0)); p=as_p((3,4,0,0),7)
    assert qfi(sector_outputs(p,s))==qfi(sector_outputs(p,s,(1,0)))==Q(36,49)
    # Both noiseless.
    s=(Q(1),Q(1)); p=as_p((1,0,0,1),2)
    assert qfi(sector_outputs(p,s))==4
    print("Noiseless/damped edge: 144/49 parallel, 625/196 noisy-first; replacer and two-unitary controls pass.")

if __name__=="__main__":
    ghz_control()
    high_loss_controls()
    boundary_controls()
