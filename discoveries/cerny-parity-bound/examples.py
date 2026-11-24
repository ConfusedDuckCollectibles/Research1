"""Construct, then independently validate plain transition/reset certificates."""
from collections import Counter, deque
from itertools import combinations, permutations
from fractions import Fraction
from pathlib import Path
import json, random, sys
from bound import parity_collision, theoretical, profile

CFG={
  'A':((2,0,0,2),(0,1)),
  'B':((2,2,0,0),(0,1)),
  'C':((2,0,0,2),(1,1)),
}

def graph_cycle_lengths(f):
    cycles=set()
    for i in range(len(f)):
        track={}
        while i not in track:
            track[i]=len(track); i=f[i]
        cycle=[q for q in track if track[q]>=track[i]]
        cycles.add(tuple(sorted(cycle)))
    return sorted(map(len,cycles))

def accepts_weights(T,excluded):
    n=len(T[0]);p=None
    AA=[Counter(f) for f in T]
    for j in range(n):
        if j in excluded:continue
        a,b=AA[0][j],AA[1][j]
        if a==b:
            if a!=1:return False
        else:
            z=Fraction(1-b,a-b)
            if not 0<z<1 or (p is not None and z!=p):return False
            p=z
    return True

def structural_gap(T,kind):
    n=len(T[0])
    C=[graph_cycle_lengths(f) for f in T]
    if any(len(L)<2 for L in C) or all(max(L)<=1 for L in C):
        return None
    if any(len(set(f))>n-2 for f in T):return None
    if accepts_weights(T,set()):return None
    AA=[Counter(f) for f in T]
    # For every candidate sole exceptional x, some node outside
    # x and its two outputs forces the letter probability to 1/2.
    if any(not any((AA[0][j],AA[1][j]) in ((2,0),(0,2))
                   for j in range(n) if j not in {x,T[0][x],T[1][x]})
           for x in range(n)):
        return None
    if any(accepts_weights(T,{x,s}) and all(f[x] in {x,s} for f in T)
           for x,s in permutations(range(n),2)):
        return None
    found=set();covered=set()
    for x,y,s in permutations(range(n),3):
        if not all(f[q] in {x,y,s} for q in (x,y) for f in T):continue
        if not accepts_weights(T,{x,y,s}):continue
        to_local={x:0,y:1,s:2}
        bd=tuple(sum(f[q]==z for f in T for q in range(n) if q not in {x,y}) for z in (x,y))
        types=[]
        for a,b in (T,T[::-1]):
            loc=tuple(to_local[f[q]] for q,f in ((x,a),(x,b),(y,a),(y,b)))
            types += [k for k,cv in CFG.items() if cv==(loc,bd)]
        witness=(frozenset({x,y}),s)
        found.add(witness)
        if types:covered.add(witness)
    if found!=covered:return None
    return C,len(found)

def strong(T):
    n=len(T[0])
    for q in range(n):
        C=[q]; S={q}
        for x in C:
            for f in T:
                y=f[x]
                if y not in S:S.add(y);C.append(y)
        if len(C)!=n:return False
    return True

def reset_word(T):
    n=len(T[0])
    P=list(combinations(range(n),2))
    R={q:[] for q in P}
    merging={}
    for q in P:
        for l,f in enumerate(T):
            a,b=q
            if f[a]==f[b]:merging[q]=str(l)
            else:R[tuple(sorted((f[a],f[b])))].append((l,q))
    Q=deque(merging)
    while Q:
        z=Q.popleft()
        for l,q in R[z]:
            if q not in merging:merging[q]=str(l)+merging[z];Q.append(q)
    if len(merging)!=len(P):return None
    S=set(range(n));U=''
    while len(S)>1:
        u=min(merging[q] for q in combinations(sorted(S),2))
        for l in map(int,u):S={T[l][q] for q in S}
        U+=u
    return U

def check(obj):
    T=tuple(tuple(f) for f in obj['T'])
    n=len(T[0]);kind=obj['kind']
    u=obj['reset']
    assert all(len(f)==n and all(0<=j<n for j in f) for f in T)
    assert strong(T)
    x,bd=CFG[kind]
    assert (T[0][0],T[1][0],T[0][1],T[1][1])==x
    counts=Counter(j for a in T for j in a[2:])
    assert (counts[0],counts[1])==bd
    assert counts[2]==2-sum(bd) and all(counts[j]==2 for j in range(3,n))
    S=[]
    for q in range(n):
        for j in u:q=T[int(j)][q]
        S.append(q)
    assert len(set(S))==1
    classes,nwit=structural_gap(T,kind)
    tau=parity_collision(T,kind)
    # Independent exhaustive layer expansion, not the 2n-node BFS distances.
    L={(0,1),(1,0),(2,0)}
    for t in range(1,tau+1):
        L|={(f[q],c^p) for q,c in L for f,p in zip(T,(1, int(kind=='B')))}
        assert (len({q for q,c in L})!=len(L))==(t==tau)
    bound,*_=theoretical(n,kind,tau)
    refined,*_=profile(T,kind)
    assert refined<=bound<=(n-1)**2
    print(kind,n,'tau',tau,'bound',bound,'profile',refined,'Cerny',(n-1)**2,
          'reset length',len(u),'cycle lengths',classes,'ranks',[len(set(f)) for f in T],
          'unordered QE witnesses',nwit)

if __name__=='__main__':
    path=Path('examples.json')
    if '--make' in sys.argv:
        n=20
        rng=random.Random(82250309)
        data=[]
        for k in CFG:
            local,indeg=CFG[k]
            for _ in range(6000):
                L=[j for j in range(3,n) for _ in range(2)]+[0]*indeg[0]+[1]*indeg[1]+[2]*(2-sum(indeg))
                rng.shuffle(L)
                T=(local[0],local[2],*L[::2]),(local[1],local[3],*L[1::2])
                tau=parity_collision(T,k)
                if tau is None or theoretical(n,k,tau)[0]>(n-1)**2:continue
                if not strong(T) or not structural_gap(T,k):continue
                u=reset_word(T)
                if u is None:continue
                data.append({'kind':k,'T':T,'reset':u})
                break
        assert len(data)==3
        path.write_text(json.dumps(data)+'\n')
    for obj in json.loads(path.read_text()):check(obj)
