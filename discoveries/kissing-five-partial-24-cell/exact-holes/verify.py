#!/usr/bin/env python3
"""Original, dependency-free exact controls for the companion theorem."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product

def need(ok, reason):
    if not ok:
        raise AssertionError(reason)

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def add(x,y):
    return tuple(a+b for a,b in zip(x,y))

def neg(x):
    return tuple(-a for a in x)

def scaled(x,k):
    return tuple(k*a for a in x)

def root(i,j,a=1,b=1):
    return tuple(a if t==i else b if t==j else 0 for t in range(4))

D=tuple(sorted(root(i,j,a,b) for i,j in combinations(range(4),2)
               for a,b in product((-1,1),repeat=2)))
axes=tuple(tuple(int(i==j)*s for i in range(4)) for j in range(4) for s in (-1,1))
cubes=tuple(product((-1,1),repeat=4))
# Coordinates are doubled: use dot=4 as the unit-squared-norm scale.
frames=(tuple(scaled(x,2) for x in axes),
        tuple(x for x in cubes if x.count(1)%2==0),
        tuple(x for x in cubes if x.count(1)%2==1))
F=tuple(sorted(sum(frames,())))
RSET=frozenset(D)

def mixed_sign(a,b):
    # a+b sqrt(2), exact with rationals.
    f=lambda z:(z>0)-(z<0)
    if not a: return f(b)
    if not b or f(a)==f(b): return f(a)
    return f(a)*f(a*a-2*b*b)

def rays(r,M):
    # (2u, a, b) encodes 2t=a+b sqrt(2).
    base=[(x,2,0) for x in F if dot(r,x)==2]
    need(len(base)==6,"base-ray count")
    apex=(scaled(r,2),0,0)
    plus=[(add(r,s),0,1) for s in M if dot(r,s)==1]
    return [apex,*base,*plus]

def ray_dot(z,w):
    u,a,b=z
    v,c,d=w
    return dot(u,v)+a*c+2*b*d,a*d+b*c

def cone_controls(M):
    for r in M:
        G=rays(r,M)
        need(all(ray_dot(g,g)==(8,0) for g in G),"ray norm")
        for g,h in combinations(G,2):
            a,b=ray_dot(g,h)
            need(mixed_sign(a-4,b)>=0,"acute generator pairing")
            if not g[1] and not h[1]:
                need(mixed_sign(a-4,b)>0,"strict mixed/apex pairing")

def adj_number(M):
    return sum(dot(r,s)==1 for r,s in combinations(M,2))

def high_controls():
    need(all(abs(dot(x,y))==2
             for A,B in combinations(frames,2) for x,y in product(A,B)),
         "triality cross-frame dots")
    need(all(len(S)==8 and all(dot(x,y) in (-4,0) for x,y in combinations(S,2))
             for S in frames),"three frames")
    def independent(S,k,pre=()):
        if not k:
            yield pre
        for i,x in enumerate(S):
            if len(S)-i<k: break
            yield from independent(tuple(y for y in S[i+1:] if dot(x,y)<=0),
                                   k-1,pre+(x,))
    data=[]
    for k,n in ((6,84),(7,24),(8,3)):
        ct=0
        for S in independent(F,k):
            ct+=1
            need(any(set(S)<=set(T) for T in frames),"large polar layer has a single frame")
        need(ct==n,"polar subset count")
        data.append(ct)
    triple=[S for S in combinations(F,3) if all(dot(x,y)==-2 for x,y in combinations(S,2))]
    need(triple and all(not any(all(dot(x,y)==0 for y in S) for x in F) for S in triple),
         "a singular polar triple has no orthogonal polar vector")
    print("6/7/8 polar layers:",*data,"; singular polar triples:",len(triple))

def high_ray_controls():
    r,s,t=root(0,1),root(0,2),root(1,2,1,-1)
    M=(r,s,t)
    need(adj_number(M)==2,"test has path holes")
    # Distinct rational controls: first sole violation of r; next double and
    # a further obeyed neighbor t. All vectors are squared norm one.
    U=[(Fraction(3,5),Fraction(4,5),0,0),
       (Fraction(20,25),Fraction(12,25),Fraction(9,25),0)]
    for w,expected in zip(U,({r},{r,s})):
        need(dot(w,w)==1 and {x for x in D if dot(x,w)>1}==expected,
             "unit high violation profile")
        need(dot(r,w)>max(dot(s,w),dot(t,w)),"maximal assignment")
        for Y in rays(r,M):
            a,b=ray_dot((scaled(w,2),2,0),Y)
            need(mixed_sign(a-4,b)>0,"nonpolar high conflicts with every type-ray")
    v=(0,0,Fraction(5,13),Fraction(12,13))
    rr=root(2,3)
    nonpolar=(axes[0],axes[1],axes[2],axes[3],(0,0,-1,0),(0,0,0,-1),v)
    need(all(dot(x,x)==1 for x in nonpolar) and
         all(dot(x,y)<=0 for x,y in combinations(nonpolar,2)),"genuine nonpolar 7-layer")
    need({x for x in D if dot(x,v)>1}=={rr},"its sole violation")
    print("Nonpolar 7-high control and W conflicts incl. obeyed adjacent hole pass")

def code_check(V, size, holes):
    eq={(*scaled(r,10),0) for r in D}
    expect={(*scaled(r,10),0) for r in holes}
    need(len(V)==size and (V&eq)==eq-expect,"exact number/template intersection")
    need(all(dot(x,x)==200 for x in V),"all scaled norms")
    need(all(dot(x,y)<=100 for x,y in combinations(V,2)),"all scaled separations")

def do_count():
    E={(*scaled(x,10),0) for x in D}
    plain={(*scaled(x,5),t) for x,t in product(frames[0],(-10,10))}
    leech={(*scaled(x,5),10) for x in frames[1]}|{(*scaled(x,5),-10) for x in frames[0]}
    D5=E|plain
    L5=E|leech
    cnt=0
    for M in combinations(D,2):
        cone_controls(M)
        rem={(*scaled(x,10),0) for x in M}
        code_check(D5-rem,38,M)
        code_check(L5-rem,38,M)
        cnt+=1
    need(cnt==276,"hole pair count")
    print("Exact 38-point witnesses for all 276 pairs, both aligned families")
    counts=Counter()
    triangle=[]
    for M in combinations(D,3):
        e=adj_number(M)
        counts[e]+=1
        if e==3:
            triangle.append(frozenset(M))
            continue
        cone_controls(M)
        rem={(*scaled(x,10),0) for x in M}
        code_check(D5-rem,37,M)
    need(counts=={0:488,1:1056,2:384,3:96},"adjacency counts")
    print("Exact 37-point witnesses for all",sum(counts[e] for e in (0,1,2)),
          "nontriangles;",sorted(counts.items()))
    return D5,set(triangle)

def triangle_transport(D5,expected):
    H=((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
    base=frozenset((root(0,1),root(0,2),root(1,2)))
    initial=(D5-{(*scaled(x,10),0) for x in base})|{(10,6,8,0,0)}
    def had(x):
        d=tuple(dot(h,x[:4]) for h in H)
        need(all(a%2==0 for a in d),"integral scaled Hadamard transport")
        return tuple(a//2 for a in d)+(x[4],)
    star=frozenset(had((*x,0))[:4] for x in base)
    need(adj_number(star)==3 and {tuple(i for i,a in enumerate(x) if a) for x in star}
         =={(0,1),(0,2),(0,3)},"Hadamard takes face to star")
    V={base:initial,star:{had(x) for x in initial}}
    found={}
    for M,code in V.items():
        for perm,signs in product(permutations(range(4)),product((-1,1),repeat=4)):
            def act(x):
                return tuple(x[i]*sg for i,sg in zip(perm,signs))
            U=frozenset(act(x) for x in M)
            if U not in found:
                found[U]={(*act(x[:4]),x[4]) for x in code}
    need(set(found)==expected,"two transported triangle orbits cover every hole triangle")
    for M,V in found.items():
        code_check(V,38,M)
    print("Transported exact 38-point constructions for all",len(found),"triangular masks")

if __name__=="__main__":
    high_controls()
    high_ray_controls()
    triangle_transport(*do_count())
