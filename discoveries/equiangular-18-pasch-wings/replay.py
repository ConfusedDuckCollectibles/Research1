#!/usr/bin/env python3
"""Exact standalone proof replay. Python standard library + C++20 compiler.
Run in this directory: python3 replay.py
Avoid Python -O and do not define NDEBUG when compiling the C++ helper.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import lcm
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import argparse

# Each R^6 block has coordinate sum zero; drop its sixth coordinate.
# Each R^2 block has coordinate sum zero; drop its second coordinate.
# The metric below reconstructs those missing coordinates.
def dp(a,b):
    return sum(x*y for x,y in zip(a,b))

def weights(a):
    u=tuple(a[i]+sum(a[(i//5)*5:(i//5)*5+5]) for i in range(15))
    return u+tuple(2*x for x in a[15:])

def metric(a,b):
    return dp(weights(a),b)

# Inputs store six times the physical vectors. Primary paper §5.
rho=(0,)*18+(6,)
wings=tuple(tuple((tuple(map(int,I)),sign*s) for I,s in M)
            for M in (
                (('135',-1),('146',-1),('236',-1),('245',-1),('123',1),('345',1),('156',1)),
                (('135',-1),('146',-1),('236',-1),('245',-1),('126',1),('234',1),('456',1)),
            ) for sign in (-1,1))
def vec(i,I,e):
    a=[0]*19
    a[5*i:5*i+5]=[3-6*(j in I) for j in (1,2,3,4,5)]
    a[15+(i+1)%3]=3*e
    a[15+(i+2)%3]=3
    a[18]=3
    return tuple(a)

def selected(ch):
    return tuple(vec(i,I,e) for i,c in enumerate(ch) for I,e in wings[c])

def orbits():
    """Concurrent block permutations and independent signs on the first
    three A1 coordinates, then signs of individual lines and S6^3."""
    names={}
    for i,M in enumerate(wings):
        for p in permutations(range(1,7)):
            tag=tuple(sorted((tuple(sorted(p[j-1] for j in I)),s) for I,s in M))
            assert tag not in names or names[tag]==i
            names[tag]=i
    assert len(names)==960
    def image(ch,p,signs):
        out=[None]*3
        for i,c in enumerate(ch):
            j=p[i]
            raw=[]
            for I,e in wings[c]:
                a=[0]*3
                u,v=p[(i+1)%3],p[(i+2)%3]
                a[u],a[v]=e*signs[u],signs[v]
                s,t=a[(j+1)%3],a[(j+2)%3]
                J=I if t==1 else tuple(x for x in range(1,7) if x not in I)
                raw.append((J,s*t))
            out[j]=names[tuple(sorted(raw))]
        return tuple(out)
    remain=set(product(range(4),repeat=3))
    result=[]
    while remain:
        rep=min(remain)
        got={image(rep,p,s) for p in permutations(range(3))
             for s in product((-1,1),repeat=3)}
        assert rep in got and got<=remain
        remain-=got
        result.append((rep,len(got)))
    assert result==[((0,0,0),8),((0,0,1),8),((0,0,2),24),((0,0,3),24)]
    return result

def independent(rows):
    result=[]
    E={}
    for row in rows:
        a=list(map(F,row))
        for j,b in sorted(E.items()):
            t=a[j]
            if t:a=[x-t*y for x,y in zip(a,b)]
        piv=next((j for j,t in enumerate(a) if t),None)
        if piv is not None:
            d=a[piv]
            E[piv]=[t/d for t in a]
            result.append(row)
    return result

def inv(A):
    n=len(A)
    T=[[F(a) for a in row]+[F(j==i) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        i=next(j for j in range(k,n) if T[j][k])
        T[i],T[k]=T[k],T[i]
        d=T[k][k]
        T[k]=[x/d for x in T[k]]
        for j in range(n):
            if j==k:continue
            t=T[j][k]
            if t:T[j]=[x-t*y for x,y in zip(T[j],T[k])]
    return [a[n:] for a in T]

def candidates(ch,exe):
    W=selected(ch)
    assert len(set(W))==21
    assert all(metric(w,w)==108 and metric(w,rho)==36 for w in W)
    assert all(metric(v,w) in (0,36) for v,w in combinations(W,2))
    B=independent((rho,)+W)
    assert len(B)==19 and B[0]==rho
    G=[]
    for a in B:
        row=[metric(a,b) for b in B]
        assert all(x%36==0 for x in row)
        G.append([x//36 for x in row])
    U=inv(G)
    # v = B * G^{-1} * z in 6x coordinates. Store columns in C.
    C=[[sum(B[j][i]*U[j][k] for j in range(19)) for i in range(19)] for k in range(19)]
    mul=lcm(*(x.denominator for a in C for x in a))
    d=6*mul
    T=[d,*[int(mul*t) for a in C for t in a],*[t for a in W for t in a]]
    output=run([exe],input=' '.join(map(str,T)),text=True,check=True,capture_output=True).stdout
    V=[tuple(map(int,row.split())) for row in output.splitlines()]
    assert len(V)==len(set(V))
    assert all(metric(v,v)==3*d*d and 2*v[-1]==d for v in V)
    return d,V,W

def make_graphs(V,d):
    m=[weights(v) for v in V]
    n=len(V)
    A=[1<<i for i in range(n)] # norm 3 is integral: include self in A.
    B=[0]*n
    for i in range(n):
        for j in range(i):
            z=dp(m[i],V[j])
            if z%(d*d)==0:
                A[i]|=1<<j;A[j]|=1<<i
                if z in (0,d*d):
                    B[i]|=1<<j;B[j]|=1<<i
    return A,B

def branches(A):
    """Group true twins (equal closed neighborhoods), then enumerate
    maximal cliques of that integrality quotient by Bron--Kerbosch."""
    groups={}
    for i,row in enumerate(A): groups.setdefault(row,[]).append(i)
    C=list(groups.values())
    R=[{j for j,W in enumerate(C) if i!=j and A[T[0]]>>W[0]&1} for i,T in enumerate(C)]
    def rec(S,P,X):
        if not(P or X):
            yield tuple(i for k in S for i in C[k])
            return
        p=max(P|X,key=lambda j:(len(P&R[j]),-j))
        for i in sorted(P-R[p]):
            yield from rec(S+(i,),P&R[i],X&R[i])
            P.remove(i);X.add(i)
    return list(rec((),set(range(len(C))),set()))

def slices(V,A,T):
    """Supply at most 36 independent groups in equiangularity graph."""
    for i,j in combinations(range(3),2):
        for mode in range(4):
            C=defaultdict(list)
            for k in T:
                a=V[k][5*i:5*i+5];b=V[k][5*j:5*j+5]
                aa=tuple(-x for x in a);bb=tuple(-x for x in b)
                key=((a,b),min((a,b),(aa,bb)),(min(a,aa),min(b,bb)),(min(a,aa),b))[mode]
                C[key].append(k)
            if len(C)<=36 and all(not (A[k]>>l&1) for G in C.values() for k,l in combinations(G,2)):
                return len(C)
    return None

def complement_components(A,T):
    """Independence numbers in the connected components of the complement."""
    adj=[0]*len(T)
    for i,x in enumerate(T):
        for j in range(i):
            if not (A[x]>>T[j]&1):
                adj[i]|=1<<j
                adj[j]|=1<<i
    todo=(1<<len(T))-1
    comp=[]
    @lru_cache(maxsize=None)
    def alpha(R):
        if not R:return 0
        q=R&-R
        v=q.bit_length()-1
        return max(alpha(R^q),1+alpha((R^q)&~adj[v]))
    while todo:
        pending=todo&-todo
        seen=0
        while pending:
            q=pending&-pending
            seen|=q
            v=q.bit_length()-1
            pending=(pending^q)|(adj[v]&~seen)
        todo&=~seen
        assert seen.bit_count()<=15
        comp.append((seen.bit_count(),alpha(seen)))
    return Counter(comp)


def control():
    # One published construction: use the cyclic Latin square with the
    # first representative of the three wing configurations.
    W=list(selected((0,0,0)))
    for i,j in product(range(1,7),repeat=2):
        k=1+(i+j-2)%6
        W.append(tuple(1-6*(x==u) for x in (i,j,k) for u in range(1,6))+(0,0,0,3))
    assert len(set(W))==57
    assert len(independent((rho,)+tuple(W)))==19
    assert all(metric(v,v)==108 and metric(v,rho)==36 for v in W)
    assert all(metric(v,u) in (0,36) for v,u in combinations(W,2))

def audit(exe):
    control()
    out=[]
    for ch,deg in orbits():
        d,V,W=candidates(ch,exe)
        A,B=make_graphs(V,d)
        hist=Counter()
        for T in branches(A):
            hist[len(T)]+=1
            if len(T)==216:
                assert slices(V,B,T)==36
                continue
            assert len(T)==93
            assert complement_components(B,T)=={(8,4):6,(15,4):3}
        row=(ch,deg,d,len(V),dict(sorted(hist.items())))
        print(row,flush=True)
        out.append(row)
    assert [(r[2],r[3],r[4]) for r in out]==[
         (36,864,{216:4}),(12,1053,{93:36,216:4}),
         (36,864,{216:4}),(12,1053,{93:36,216:4})]
    print('PASS: all three-wing extensions have at most 57; known 57 positive control checked.')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--cxx',default='c++')
    p.add_argument('--binary',type=Path)
    args=p.parse_args()
    if args.binary: audit(str(args.binary.resolve()))
    else:
        with TemporaryDirectory() as tf:
            exe=str(Path(tf)/'enumerate')
            source=Path(__file__).with_name('enumerate.cpp')
            run([args.cxx,'-O3','-std=c++20',str(source),'-o',exe],check=True)
            audit(exe)
