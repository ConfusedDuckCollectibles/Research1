#!/usr/bin/env python3
"""Independent exact finite controls using Fraction, no solver and no patient data."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from random import Random

def path_words(n, arcs, cap=3):
    id = {a:i for i,a in enumerate(arcs)}
    for k in range(1,cap+1):
        for w in permutations(range(n),k):
            a = tuple(id.get((v,u)) for v,u in zip(w,w[1:]))
            if None not in a:
                yield w,a

def payoff(path, bits):
    val = 1 if path else 0
    for a in path:
        if not bits[a]:
            break
        val += 1
    return val

def cut_case(n, arcs, expected):
    full = [a for w,a in path_words(n,arcs) if len(w)==3]
    H = set(tuple(sorted(a)) for a in full)
    assert len(H)==len(full)
    k,m = len(H),len(arcs)
    all_splits = list(product((0,1),repeat=m))
    cuts = {x:sum(x[i]!=x[j] for i,j in H) for x in all_splits}
    x = max(cuts,key=cuts.__getitem__)
    neg = tuple(1-b for b in x)
    mincut = F(3,2)+F(k-cuts[x],2*k)
    direct = F(sum(payoff(a,x)+payoff(a,neg) for a in full),2*k)
    assert mincut == direct == expected
    if k%2:
        shifted = [t for t in all_splits if cuts[t]==k-1]
        # Uniform maximum cuts and their orientations give matching upper
        # bound on *each* full chain in all two odd-cycle controls.
        assert all(F(sum(payoff(a,t) for t in shifted),len(shifted))==expected for a in full)
    else:
        assert all(F(payoff(a,x)+payoff(a,neg),2)==expected for a in full)
    print("derived H",sorted(H),"exact game",expected)

def inv(A):
    n=len(A)
    t=[[F(x) for x in A[j]]+[F(i==j) for i in range(n)] for j in range(n)]
    for i in range(n):
        k=next((k for k in range(i,n) if t[k][i]),None)
        if k is None:
            return None
        t[i],t[k]=t[k],t[i]
        p=t[i][i]
        t[i]=[v/p for v in t[i]]
        for j in range(n):
            if i!=j:
                p=t[j][i]
                t[j]=[x-p*y for x,y in zip(t[j],t[i])]
    return [row[n:] for row in t]

def all_vertex_triangle_laws(states):
    rows=list(zip(*[(1,*m) for m in states]))
    basis=[]
    for inds in combinations(range(8),4):
        Ai=inv([[row[j] for j in inds] for row in rows])
        if Ai:
            basis.append((inds,Ai))
    assert len(basis)==58
    def at(p):
        rhs=(F(1),p,p,p)
        for inds,Ai in basis:
            x=[sum(v*b for v,b in zip(row,rhs)) for row in Ai]
            if min(x)>=0:
                q=[F(0)]*8
                for i,xi in zip(inds,x):
                    q[i]=xi
                assert sum(q)==1
                assert all(sum(t[i]*p0 for t,p0 in zip(states,q))==p for i in range(3))
                yield q
    return at

def triangle():
    states=list(product((0,1),repeat=3))
    vertices=all_vertex_triangle_laws(states)
    rng=Random(1873)
    for p in {F(i,60) for i in range(61)} | {F(9,20),F(11,20)}:
        states_h=[(x[0]*x[1],x[1]*x[2],x[2]*x[0]) for x in states]
        ext = list(vertices(p))
        h_laws=[tuple(sum(q*h[i] for q,h in zip(law,states_h)) for i in range(3)) for law in ext]
        f=max(F(0),2*p-1)
        g=max(F(0),3*p-1-3*f)
        for w0 in [(1,1,1),(3,0,0),tuple(rng.randint(1,40) for j in range(3))]:
            w=[F(x,sum(w0)) for x in w0]
            worst=min(sum(a*b for a,b in zip(w,h)) for h in h_laws)
            assert worst==f+g*min(w)
    p=F(1,2)
    w=[F(7,15),F(4,15),F(4,15)]
    R=[]
    for law in vertices(p):
        R.append([sum(q*(w[i]+w[(i-1)%3]*bits[(i-1)%3]+
                            w[(i-2)%3]*bits[(i-2)%3]*bits[(i-1)%3])
                      for q,bits in zip(law,states)) for i in range(3)])
    guarantees=tuple(min(r[i] for r in R) for i in range(3))
    val=min(sum(r) for r in R)
    assert guarantees==(F(3,5),F(1,2),F(2,5)) and val==F(49,30)
    print("triangle p-grid exact 61 points x three tests; priority",guarantees,val)

def two_chain_control():
    # Use two disjoint sources: length 3, length 2, and disjoint length-2 cycle.
    arcs = [(0,1),(1,2),(3,4),(5,6),(6,5)]
    total,both = 0,[0]*len(arcs)
    for order in permutations(range(7)):
        rank = {x:i for i,x in enumerate(order)}
        x=tuple(rank[u]<rank[v] for u,v in arcs)
        total+=payoff((0,1),x)+payoff((2,),x)+2*(x[3] and x[4])
        both=[a+int(v) for a,v in zip(both,x)]
    assert set(both)=={2520} and F(total,5040)==F(19,6)
    coupled=[]
    for a,b,c in product((0,1),repeat=3):
        x=(a,1-a,b,c,1-c)
        coupled.append(payoff((0,1),x)+payoff((2,),x)+2*(x[3] and x[4]))
    assert F(sum(coupled),len(coupled))==F(3)
    print("two-chain + atomic cycle: common rank law 19/6, pure adverse coupling 3")

if __name__=="__main__":
    cut_case(3,[(0,1),(1,2),(2,0)],F(5,3))
    cut_case(5,[(0,1),(1,2),(2,3),(3,4),(1,3)],F(8,5))
    cut_case(4,[(0,1),(1,2),(2,3),(3,0)],F(3,2))
    triangle()
    two_chain_control()
