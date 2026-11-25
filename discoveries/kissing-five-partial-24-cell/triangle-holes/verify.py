#!/usr/bin/env python3
"""Exact finite controls for the triangular-hole continuum theorem."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations, product

def must(ok, name):
    if not ok:
        raise AssertionError(name)

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def add(x,y):
    return tuple(a+b for a,b in zip(x,y))

def scale(x,m):
    return tuple(m*a for a in x)

def supp(x):
    return {i for i,a in enumerate(x) if a}

def root(i,j,a,b):
    return tuple(a if k==i else b if k==j else 0 for k in range(4))

roots={root(i,j,a,b) for i,j in combinations(range(4),2)
       for a,b in product((-1,1),repeat=2)}
ax={tuple(int(i==j)*s for i in range(4)) for j in range(4) for s in (-1,1)}
cube={scale(s,Q(1,2)) for s in product((-1,1),repeat=4)}
even={x for x in cube if sum(a<0 for a in x)%2==0}
F=(ax,even,cube-even)
polar=set().union(*F)
mat=((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
H=lambda x:tuple(Q(dot(row,x),2) for row in mat)
J=lambda x:(x[0],x[1],x[2],-x[3])

def structure_and_cones():
    triangles={frozenset(M) for M in combinations(roots,3)
               if all(dot(r,s)==1 for r,s in combinations(M,2))}
    must(len(triangles)==96,"all triangles")
    def face_kind(M):
        S=set().union(*(supp(x) for x in M))
        p=set.intersection(*(supp(x) for x in M))
        must(len(S)==3 or len(p)==1,"pairwise adjacent support type")
        return len(S)==3
    transforms=(lambda x:x,H,lambda x:H(J(x)))
    for T,V in zip(transforms,F):
        must({T(x) for x in roots}==roots and {T(x) for x in V}==ax,
             "triality conjugation")
    for M in triangles:
        must(sum(face_kind({T(x) for x in M}) for T in transforms)==1,
             "unique face-relative frame")
        for r in M:
            old=[(r,0,0)]+[(v,1,0) for v in polar if dot(v,r)==1]
            extra=[(scale(add(r,s),Q(1,2)),0,Q(1,2)) for s in M-{r}]
            C=old+extra
            must(len(old)==7,"old ray count")
            for U,V in combinations(C,2):
                a,b=qd(U,V,2)
                must(sign(a-1,b,2)>=0,"compatible ray pair")
                if U[1]==V[1]==0:
                    must(sign(a-1,b,2)>0,"strict low single-ray pairing")
            for w in polar:
                if dot(w,r)!=1: continue
                for V in old:
                    a,b=qd((w,1,0),V,2)
                    must(sign(a-1,b,2)>=0,"matching polar high/base")
                for V in extra:
                    a,b=qd((w,1,0),V,2)
                    must(sign(a-1,b,2)>0,"matching polar high/extra")
    print("All 96 hole triangles: unique face frame, both W rays and polar matching")
    return triangles

def sign(a,b,n):
    f=lambda x:(x>0)-(x<0)
    if not a: return f(b)
    if not b or f(a)==f(b): return f(a)
    return f(a)*f(a*a-n*b*b)

def qd(U,V,n):
    # (u,a,b) means horizontal u and height a+b sqrt(n).
    u,a,b=U
    v,c,d=V
    return dot(u,v)+a*c+n*b*d,a*d+b*c

def qcode(V,n,N,M):
    must(len(V)==N,"code size")
    must(all(qd(U,U,n)==(2,0) for U in V),"squared norm")
    must(all(sign(a-1,b,n)<=0 for a,b in (qd(U,W,n) for U,W in combinations(V,2))),
         "all pair separations")
    must({u for u,a,b in V if not a and not b and u in roots}==roots-M,
         "exact named intersection")

def layer_controls():
    profiles=Counter()
    for S in combinations(sorted(polar),5):
        if any(scale(u,-1) in S for u in S): continue
        if any(dot(u,v)>0 for u,v in combinations(S,2)): continue
        chunks=sorted(([u for u in S if u in T] for T in F),key=len)
        profiles[tuple(map(len,chunks))]+=1
        must(tuple(map(len,chunks))==(0,1,4),"five-polar profile")
        total=tuple(sum(v[i] for v in chunks[2]) for i in range(4))
        must(chunks[1][0]==scale(total,Q(-1,2)),"positive dependence")
    must(profiles=={(0,1,4):48},"count of five-polar cases")
    r,s,t=(1,1,0,0),(1,0,1,0),(0,1,1,0)
    T={r,s,t}
    examples=(((Q(20,25),Q(12,25),Q(9,25),0),{r,s},{r}),
              ((Q(6,11),Q(6,11),Q(7,11),0),T,T))
    for w,bad,block in examples:
        must(dot(w,w)==1 and {v for v in roots if dot(w,v)>1}==bad,
             "exact nonpolar high violation profile")
        for rr in block:
            cones=[(rr,0,0)]+[(v,1,0) for v in polar if dot(v,rr)==1]
            cones +=[(scale(add(rr,ss),Q(1,2)),0,Q(1,2)) for ss in T-{rr}]
            for Y in cones:
                a,b=qd((w,1,0),Y,2)
                must(sign(a-1,b,2)>0,"nonpolar high blocks full cone")
    print("Five polar highs without an antipode: 48 cases; exact mixed/triple high-ray controls")

def lower_and_negatives(triangles):
    r,s,t=(1,1,0,0),(1,0,1,0),(0,1,1,0)
    M={r,s,t}
    eq=[(u,0,0) for u in roots-M]
    top=[(w,1,0) for w in ax]
    both=top+[(w,-1,0) for w in ax]
    around=((1,Q(3,5),Q(4,5),0),0,0)
    canonical=eq+both+[around]
    qcode(canonical,2,38,M)
    images={}
    for T in (lambda x:x,H):
        for p,sig in product(permutations(range(4)),product((-1,1),repeat=4)):
            def act(x): 
                z=T(x)
                return tuple(z[i]*a for i,a in zip(p,sig))
            K=frozenset(act(x) for x in M)
            images.setdefault(K,[(act(u),a,b) for u,a,b in canonical])
    must(set(images)==triangles,"all triangle configurations transported")
    for K,V in images.items():
        qcode(V,2,38,K)
        must(sum(a==1 for u,a,b in V)==sum(a==-1 for u,a,b in V)==8,
             "both high layers full")
    print("Transported exact 38-point both-full code for each of 96 hole masks")
    u=(Q(3,5),Q(3,5),Q(3,5),0)
    pair=[(u,0,Q(1,5)),(u,0,Q(-1,5))]
    qcode(eq+pair,23,23,M)
    must(qd(*pair,23)==(Q(4,25),0),"opposite triple dot")
    x=((Q(4,5),Q(4,5),Q(4,5),Q(1,5)),Q(1,5),0)
    y=((Q(8,15),Q(8,15),Q(8,15),Q(-7,15)),0,Q(-1,15))
    qcode(eq+top+[x,y],209,31,M)
    must(qd(x,y,209)==(Q(89,75),Q(-1,75)),"one-full opposition dot")
    must(sign(-1,Q(1,15),209)<0 and sign(-Q(8,15),Q(1,15),209)>0,"other low height")
    must(all(r0=={z for z in roots if dot(z,U[0])>1}
             for r0,U in zip((M,M),(x,y))),"the two lows are triple")
    print("Exact two-sided opposition controls: one full cap can coexist with two triple lows")

if __name__=="__main__":
    Ts=structure_and_cones()
    layer_controls()
    lower_and_negatives(Ts)
