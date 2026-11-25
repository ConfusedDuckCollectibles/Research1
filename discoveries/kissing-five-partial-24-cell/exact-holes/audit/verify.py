"""Independent exact finite controls for named 24-cell hole bounds."""
from fractions import Fraction as Q
from itertools import combinations,product,permutations
from collections import Counter

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def plus(a,b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
axis={tuple(s if i==j else 0 for i in range(4)) for j in range(4) for s in (-1,1)}
D={plus(a,b) for a,b in combinations(axis,2) if dot(a,b)==0}
H=set(product((Q(-1,2),Q(1,2)),repeat=4))
even={x for x in H if sum(a<0 for a in x)%2==0}
frames=(axis,even,H-even)
F=set().union(*frames)
assert len(D)==len(F)==24
for i,j in combinations(range(3),2):
    assert all(abs(dot(x,y))==Q(1,2) for x in frames[i] for y in frames[j])
assert all({u for u in all_frame if dot(u,r)==1} for r in D for all_frame in frames)
cap_signs=lambda frame,h:{x+(h,) for x in frame}

def check(K,N):
    assert len(K)==N and all(dot(x,x)==2 for x in K)
    assert all(dot(x,y)<=1 for x,y in combinations(K,2))

equator={r+(0,) for r in D}
D5=equator|cap_signs(axis,1)|cap_signs(axis,-1)
check(D5,40)
for size in (2,3):
    for holes in combinations(D,size):
        K=D5-{r+(0,) for r in holes}
        assert len(K)==40-size and len(K&equator)==24-size
print("Exact sharp deletion witnesses: all 276 pairs and all 2024 triples")

def lows_rays(r,M):
    # Each height is u+v sqrt(2); all ray horizontal projections are rational.
    apex=(r,Q(0),Q(0))
    base=[(b,Q(1),Q(0)) for b in F if dot(r,b)==1 and (b in axis or b in H)]
    assert len(base)==6
    ws=[(tuple(Q(u+v,2) for u,v in zip(r,t)),Q(0),Q(1,2))
        for t in M-{r} if dot(r,t)==1]
    return [apex]+base+ws

def raydot(a,b):
    p,u,v=a
    q,x,y=b
    return dot(p,q)+u*x+2*v*y,u*y+v*x
def pos_pair(a,b,strict=False): # a+b sqrt(2)>=0 or >0; here b>=0
    assert b>=0
    if a>=0:
        return (a>0 or b>0) if strict else True
    return 2*b*b>a*a if strict else 2*b*b>=a*a

triple_count=Counter()
for holes in combinations(D,3):
    M=set(holes)
    ed=sum(dot(r,t)==1 for r,t in combinations(M,2))
    triple_count[ed]+=1
    for r,t in combinations(M,2):
        if dot(r,t)==0:
            # A fully retained pair with the same sum certifies simultaneous
            # strict violation impossible, independent of heights.
            assert any(plus(p,q)==plus(r,t) for p,q in combinations(D-M,2))
    if ed<3:
        for r in M:
            rays=lows_rays(r,M)
            assert all(raydot(p,p)==(2,0) for p in rays)
            assert all(pos_pair(a-1,b) for a,b in
                       (raydot(p,q) for p,q in combinations(rays,2)))
            single_lows=[p for p in rays if p[1]==0]
            assert all(pos_pair(a-1,b,True) for a,b in
                       (raydot(p,q) for p,q in combinations(single_lows,2)))
assert triple_count=={0:488,1:1056,2:384,3:96}
print("All 2024 triple masks by adjacency edge count:",dict(sorted(triple_count.items())))
print("Nonadjacent zero-dot never co-violates with <=3 holes; all 1928 nontriangle three-hole ray systems pass")

def all_nonpositive(cand,k,current=()):
    if k==0:
        yield current
        return
    for i,w in enumerate(cand):
        if i+k>len(cand): return
        yield from all_nonpositive([u for u in cand[i+1:] if dot(u,w)<=0],k-1,current+(w,))
for k in (6,7,8):
    C=list(all_nonpositive(sorted(F),k))
    assert C
    assert all(any(neg(w) in c for w in c) for c in C)
    assert all(any(set(c)<=f for f in frames) for c in C)
    print("Exact nonpositive subsets of polar F at size",k,":",len(C),"; all within one frame")

r=(0,0,1,1); s=(1,1,0,0)
layer={w+(1,) for w in [(1,0,0,0),(-1,0,0,0),(0,1,0,0),(0,-1,0,0),
                         (0,0,-1,0),(0,0,0,-1),(0,0,Q(3,5),Q(4,5))]}
check((equator-{r+(0,),s+(0,)})|layer,29)
check((equator-{r+(0,)})|layer,30)
assert {v for v in D if dot(v,(0,0,Q(3,5),Q(4,5)))>1}=={r}
print("Nonpolar seven-high layer exists; no assumption that all high are polar")
r,s,t=(1,1,0,0),(1,0,1,0),(0,1,1,0)
example=(1,Q(3,5),Q(4,5),0,0)
T=(D5-{a+(0,) for a in (r,s,t)})|{example}
check(T,38); assert len(T&equator)==21
assert {a for a in D if dot(a,example[:4])>1}=={r,s,t}
# A second same-chamber point without caps, x_h=(1/10,0,1,0), t_x=sqrt(99)/10.
xh=(Q(1,10),0,1,0)
assert dot(xh,xh)+Q(99,100)==2
assert {a for a in D if dot(a,xh)>1}=={s}
assert max((dot(a,example[:4]),a) for a in D)[1]==s
assert dot(xh,example[:4])==Q(9,10)
print("Triangle: exact 38-point counterexample to universal 37; per-type low conflict also fails at dot 9/10")
