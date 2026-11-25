#!/usr/bin/env python3
"""Exact finite identities used by the 22/23-root continuum proof."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product

def require(test, message):
    if not test:
        raise AssertionError(message)

def roots(d):
    return {tuple(a if j==i else b if j==k else 0 for j in range(d))
            for i,k in combinations(range(d),2) for a,b in product((-1,1),repeat=2)}

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def basis(d):
    return tuple(tuple(int(i==j) for j in range(d)) for i in range(d))

def neg(x):
    return tuple(-i for i in x)

def vec_add(x,y):
    return tuple(a+b for a,b in zip(x,y))

def sum_plus_sqrt2_sign(a,b):
    # Sign of a+b*sqrt(2) without evaluating sqrt(2).
    s=lambda x: int(x>0)-int(x<0)
    if not a:
        return s(b)
    if not b or s(a)==s(b):
        return s(a)
    return s(a)*s(a*a-2*b*b)

def qdot(x,y):
    a=sum(p*r+2*q*s for (p,q),(r,s) in zip(x,y))
    b=sum(p*s+q*r for (p,q),(r,s) in zip(x,y))
    return a,b

def lift(x,t):
    return tuple((F(z),F(0)) for z in x)+((F(t),F(0)),)

def triality():
    D=roots(4)
    E=basis(4)
    A={v for e in E for v in (e,neg(e))}
    H=[{tuple(F(s,2) for s in a) for a in product((-1,1),repeat=4)
        if sum(s<0 for s in a)%2==j} for j in range(2)]
    frames=(A,*H)
    full=A|H[0]|H[1]
    require(len(full)==24 and all(dot(u,u)==1 for u in full),"exact polar candidates")
    require(all(dot(u,r)<=1 for u in full for r in D),"polar feasibility")
    require(all((dot(x,y) in (-1,0)) for V in frames for x,y in combinations(V,2)),
            "orthoplex within each frame")
    require(all(abs(dot(x,y))==F(1,2) for S,T in combinations(frames,2)
                for x,y in product(S,T)),"incompatible mixed frames")
    require(all(D <= {vec_add(x,y) for x,y in combinations(V,2) if dot(x,y)==0}
                for V in frames),"every root splits in every frame")
    T=tuple(tuple(F(s,2) for s in a) for a in
            ((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)))
    standard=lambda x:tuple(dot(r,x) for r in T)
    flip=lambda x:(-x[0],*x[1:])
    require({standard(x) for x in D}==D and {flip(x) for x in D}==D,"equator symmetry")
    for trans, permutation in ((standard,(1,0,2)),(flip,(0,2,1))):
        require(all({trans(x) for x in frames[i]}==frames[j] for i,j in enumerate(permutation)),
                "triality permutation")
    holes=Counter(dot(r,s) for r,s in combinations(D,2))
    require(holes=={-2:12,-1:96,0:72,1:96},"all pair relations")
    pairset=({(0,1),(2,3)},{(0,2),(1,3)},{(0,3),(1,2)})
    for r,s in combinations(D,2):
        blocked={tuple(i for i,x in enumerate(v) if x) for v in (r,s)}
        require(any(not P&blocked for P in pairset),"intact matching")
        for a,b in ((r,s),(s,r)):
            if dot(a,b)<=0:
                require(b not in {v for v in D if dot(a,v)==1},"retained cone inequalities")
    print("24 polar directions and all 3 frames; missing pair counts",sorted(holes.items()))
    return D,frames

def cone(D):
    r,s=(1,1,0,0),(1,0,1,0)
    e1,e2,*_=basis(4)
    old=[lift(r,0),lift(e1,1),lift(e2,1)]
    old+=[lift((F(1,2),F(1,2),F(a,2),F(b,2)),1) for a,b in product((-1,1),repeat=2)]
    W=lift((1,F(1,2),F(1,2),0),0)[:-1]+((F(0),F(1,2)),)
    qnorm=lambda x:qdot(x,x)
    require(all(qnorm(x)==(2,0) for x in [*old,W]),"sphere generator norms")
    allpairs=[qdot(x,y) for x,y in combinations(old,2)]
    require(min(a for a,b in allpairs)==1 and all(b==0 for a,b in allpairs),"original cone")
    require(all(sum_plus_sqrt2_sign(a-1,b)>0 for a,b in (qdot(W,x) for x in old)),
            "additional generator endpoints")
    require(qdot(W,old[0])==(F(3,2),0),"two possible single-ray lows")
    require({p for p in D if sum_plus_sqrt2_sign(qdot(W,lift(p,0))[0]-1,0)>0}=={r,s},
            "positive test for genuine double violation")
    require(sum_plus_sqrt2_sign(-F(3,2),F(1,2))<0,
            "mixed generator fails original single-cone height criterion")
    print("seven-ray and eight-ray compatibility; exact mixed low violates only the adjacent two roots")

def sharp_controls(D,frames):
    r={(*v,0) for v in D}
    lower={(*u,-1) for u in frames[0]}
    d5=r|lower|{(*u,1) for u in frames[0]}
    l5=r|lower|{(*u,1) for u in frames[1]}
    require(d5==roots(5),"D5 construction")
    for label,S,antipodes in (("D5",d5,40),("L5",l5,24)):
        require(len(S)==40 and all(dot(x,x)==2 for x in S),"witness norm or size")
        require(max(dot(x,y) for x,y in combinations(S,2))==1,"witness separation")
        require(sum(neg(x) in S for x in S)==antipodes,"distinct isometry classes")
        print(label,"exact 40-point witness;",antipodes,"points with an antipode")

if __name__=="__main__":
    D,frames=triality()
    cone(D)
    sharp_controls(D,frames)
