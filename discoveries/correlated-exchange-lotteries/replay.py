#!/usr/bin/env python3
"""Exact certificates for synthetic common-law prefix exchanges."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, combinations, product
from math import comb, factorial
from pathlib import Path

def candidates(E,n,L=3):
    # Any first recipient is reached via a certain source edge.
    return [v for v in permutations(range(n),L) if all(e in E for e in zip(v,v[1:]))]

def maxcut_gaps(E,P,w):
    opts={e:i for i,e in enumerate(E)}
    props=[tuple(opts[e] for e in zip(v,v[1:])) for v in P]
    assert all(len(v)==2 and v[0]!=v[1] for v in props)
    assert len(w)==len(P) and sum(w)==1
    vals=[]
    for mode in product(range(2),repeat=len(E)):
        vals.append(sum(x for (a,b),x in zip(props,w) if mode[a]==mode[b]))
    return F(1,2)*min(vals),props

def order_law(n,E):
    D=Counter()
    for order in permutations(range(n)):
        rank={v:i for i,v in enumerate(order)}
        mask=sum(1<<j for j,(u,v) in enumerate(E) if rank[u]<rank[v])
        D[mask]+=1
    return {z:F(c,factorial(n)) for z,c in D.items()}

def pmarg(law,m):
    assert all(p>=0 for p in law.values()) and sum(law.values())==1
    return tuple(sum(p for z,p in law.items() if z>>i&1) for i in range(m))

def risk_curve(p):
    f=max(F(0),2*p-1)
    g=max(F(0),3*p-1-3*f)
    return 1+p+f,1+p+f+g/3
def triangle_law(p,i):
    # Unequal lottery response for the conjunction of e_i and e_(i+1)%3.
    if p<=F(1,3):return {0:1-3*p,1:p,2:p,4:p}
    if p>=F(2,3):return {7:3*p-2,3:1-p,5:1-p,6:1-p}
    h=[max(F(0),2*p-1)]*3;h[i]+=3*p-1-sum(h)
    law={(1<<j)|(1<<((j+1)%3)):h[j] for j in range(3)}
    law.update({1<<j:p-h[j]-h[(j-1)%3] for j in range(3)})
    assert pmarg(law,3)==(p,p,p)
    return law
def triangle_read(law,r):
    vec=[sum(p for z,p in law.items() if z>>i&1) for i in range(3)]
    return [r[i]+vec[(i-1)%3]*r[(i-1)%3]+r[(i-2)%3]*sum(p for z,p in law.items() if (z>>((i-1)%3)&1) and (z>>((i-2)%3)&1)) for i in range(3)]
def exact_bench():
    cases=[
      ("triangle",3,[(0,1),(1,2),(2,0)],F(1,6)),
      ("acyclic shortcut",5,[(0,1),(1,2),(2,3),(3,4),(1,3)],F(1,10)),
      ("even circle",4,[(0,1),(1,2),(2,3),(3,0)],F(0)),
      ("simple line",4,[(0,1),(1,2),(2,3)],F(0))]
    for lab,n,E,bon in cases:
        P=candidates(E,n)
        w=[F(1,len(P))]*len(P)
        b,props=maxcut_gaps(E,P,w)
        assert b==bon
        order=order_law(n,E)
        assert pmarg(order,len(E))==(F(1,2),)*len(E)
        for i,j in props:assert sum(p for z,p in order.items() if z>>i&1 and z>>j&1)==F(1,6)
        print(lab,"valid full chains",P,"fragility pairs",props,"exact uniform robust",F(3,2)+b,"ratio",1+b/F(3,2))
    for p in [F(i,20) for i in range(21)]:
        f,q=risk_curve(p)
        for i in range(3):
            law=triangle_law(p,i)
            assert pmarg(law,3)==(p,p,p)
            assert sum(triangle_read(law,[F(1,3)]*3))==q
            assert sum(triangle_read(law,[F(j==((i+1)%3)) for j in range(3)]))==f
        if F(9,20)<=p<=F(11,20):
            assert q/f>=F(21,20)
            print("triangle p",p,"pure/mixed",f,q,"factor",q/f)
    w=[F(7,15),F(4,15),F(4,15)]
    receipts=[triangle_read(triangle_law(F(1,2),j),w) for j in range(3)]
    p=[min(v) for v in zip(*receipts)]
    assert p==[F(3,5),F(1,2),F(2,5)]
    q=min(sum(v) for v in receipts)
    assert q==F(49,30) and q/F(3,2)>=F(21,20)
    print("priority lottery",w,"individual robust guarantees",p,"unit robust",q,"unconstrained pure benchmark",F(3,2))

def complete_table():
    lines=iter(Path("complete.txt").read_text().splitlines())
    for header in lines:
        n,m,fac,count=map(int,header.split())
        vals=[];raw=[]
        E=[(i,j) for i in range(n) for j in range(n) if i!=j]
        assert fac==factorial(n) and count==2**len(E)
        P=[tuple(E.index(e) for e in zip(v,v[1:])) for v in permutations(range(n))]
        for k in range(m+1):
            nn,score,mask,c=map(int,next(lines).split())
            assert nn==k==mask.bit_count() and c==comb(m,k)
            actual=0
            for path in P:
                actual+=1
                for e in path:
                    if not (mask>>e&1):break
                    actual+=1
            assert score==actual
            vals.append(F(score,fac));raw.append(mask)
        x=m//2
        v=vals[x];which=(x,x)
        for i in range(x):
            for j in range(x+1,m+1):
                q=F((j-x),j-i)*vals[i]+F((x-i),j-i)*vals[j]
                if q<v:v,which=q,(i,j)
        cap=sum(F(1,factorial(j)) for j in range(1,n+1))
        assert v=={3:F(5,3),4:F(5,3),5:F(17,10)}[n] and v<=cap
        assert which==(x,x)
        mask=raw[x];law=[e for i,e in enumerate(E) if mask>>i&1]
        print("complete",n,"all edge states",count,"value",v,"rank-law upper",cap,"min adverse seed before relabeling",law)

if __name__=="__main__":
    exact_bench();complete_table()
