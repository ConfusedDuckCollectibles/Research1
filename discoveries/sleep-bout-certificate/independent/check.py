#!/usr/bin/env python3
"""Independent synthetic-only replay, reconstructed after workspace reset."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product
from random import Random

N=960
half=Q(1,2)

def track(n, intervals):
    x=[0]*n
    for a,b in intervals:
        assert 0<=a<b<=n and not any(x[a:b])
        x[a:b]=[1]*(b-a)
    return tuple(x)

def runs(x):
    i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[j]==x[i]: j+=1
        yield x[i],i,j
        i=j

def dwell(x):
    return Counter((v,b-a) for v,a,b in runs(x))

def bouts(x,minimum):
    return [(a,b) for v,a,b in runs(x) if v and b-a>=minimum]

def hist(rs,ps,k):
    s=tuple(2*r+p for r,p in zip(rs,ps,strict=True))
    return tuple(Counter(s[i:i+j] for i in range(len(s)-j+1))
                 for j in range(1,k+1))

def iou(r,p):
    a,b=r; c,d=p
    h=max(0,min(b,d)-max(a,c))
    return Q(h,b-a+d-c-h)

def stats(r,p,minimum=10):
    r,p=bouts(r,minimum),bouts(p,minimum)
    graph=[[j for j,pp in enumerate(p) if iou(rr,pp)>=half] for rr in r]
    owner={}
    def augment(i,seen):
        for j in graph[i]:
            if j in seen: continue
            seen.add(j)
            if j not in owner or augment(owner[j],seen):
                owner[j]=i
                return True
        return False
    hits=sum(augment(i,set()) for i in range(len(r)))
    assert hits==sum(map(len,graph))
    return hits,len(r),len(p)

def same(a,b,K,full_dwell,ref_fixed=False):
    assert {len(x) for x in a+b}=={len(a[0])}
    assert all(x[0]==x[-1]==0 for x in a+b)
    assert hist(*a,K)==hist(*b,K)
    if ref_fixed: assert a[0]==b[0]
    for x,y in zip(a,b,strict=True):
        first,second=dwell(x),dwell(y)
        if not full_dwell:
            first={k:v for k,v in first.items() if k[0]}
            second={k:v for k,v in second.items() if k[0]}
        assert first==second

def confusion(fixture):
    c=hist(*fixture,1)[0]
    return tuple(c[(i,)] for i in range(4))

def shift(x,L,fixed):
    s=72 if fixed else 90
    r,p=[],[]
    for z in x:
        r.append((s,s+L))
        p.append((s+L-z,s+2*L-z))
        s+=100 if fixed else 2*L-z+90
    return track(N,r),track(N,p)

def ceildiv(x,y): return -((-x)//y)

def shift_hits(B,S,a,b,q):
    if not B*a<=S<=B*b: return set()
    if q<=a: return {B}
    if q>b: return {0}
    lo=max(0,ceildiv(S-B*(q-1),b-q+1))
    hi=min(B,(S-B*a)//(q-a))
    return set(range(lo,hi+1))

def test_shift():
    fa,fb=shift((10,)*6,16,False),shift((11,)*5+(5,),16,False)
    same(fa,fb,5,False)
    assert (stats(*fa),stats(*fb))==((0,6,6),(5,6,6))
    assert confusion(fa)==(828,36,36,60)
    assert shift_hits(6,60,4,12,11)==set(range(6))
    for B in range(1,6):
        for a,b in ((1,1),(0,3),(2,5)):
            for q in range(a-1,b+3):
                reference=defaultdict(set)
                for x in product(range(a,b+1),repeat=B):
                    reference[sum(x)].add(sum(v>=q for v in x))
                for S in range(B*a-1,B*b+2):
                    assert shift_hits(B,S,a,b,q)==reference[S]
    print('Stage 1: matched K=5 histogram, W dwells, confusion (828,36,36,60); 0/6 vs 5/6; exact hit set 0..5.')

    x=(3,6,5,6,6,6,6,6)
    y=(3,3,4,7,7,7,7,6)
    fa,fb=shift(x,10,True),shift(y,10,True)
    same(fa,fb,4,True,True)
    assert (stats(*fa),stats(*fb))==((0,8,8),(4,8,8))
    assert confusion(fa)==(844,36,36,44)
    assert dwell(fa[1])==Counter({(1,10):8,(0,90):4,(0,87):1,(0,91):1,(0,89):1,(0,79):1,(0,174):1})

    def key(x):
        return (x[0],x[-1],sum(x),tuple(sorted(v-u for u,v in zip(x,x[1:]))))

    selected=key(x); counts=defaultdict(set); witnesses=[]
    for z in product(range(3,8),repeat=8):
        S=sum(z)
        if 2*S<80 or 100*(N-2*(80-S))<88*N: continue
        k=key(z); h=sum(v>=7 for v in z)
        counts[k].add(h)
        if k==selected:
            fz=shift(z,10,True)
            same(fa,fz,4,True,True)
            assert stats(*fz)[0]==h
            witnesses.append(h)
    assert len(counts)==31752
    assert len(witnesses)==11 and set(witnesses)==set(range(5))
    assert max(max(h)-min(h) for h in counts.values() if min(h)==0 and max(h)>=4)==4
    print('Stage 2: fixed reference, all dwells, K=4 histogram, confusion (844,36,36,44); 0/8 vs 4/8.')
    print('Stage 2: scanned all 390625 absolute states; 31752 eligible summary classes, 11 selected trajectories, exact hit set 0..4.')

def trade(m,T):
    shape=Counter(u+v for u in range(m) for v in range(m))
    a=[(T+k,T+k+2) for k,d in shape.items() for _ in range(d)]
    b=[(T+k+1,T+k+1) for k,d in shape.items() for _ in range(d)]
    b.remove((T+m,T+m))
    b.append((T,T+2*m))
    order=lambda w:(sum(w),w[0])
    a.sort(key=order); b.sort(key=order)
    assert len(a)==len(b)==m*m
    assert tuple(map(sum,a))==tuple(map(sum,b))
    assert Counter(k for w in a for k in w)==Counter(k for w in b for k in w)
    return a,b

def layout(rows,K,mirror,target=0):
    gap=max(20,K+2)
    r,p=[],[]
    s=gap
    for side in range(1+mirror):
        for a,b in rows:
            single=[(s,s+a+b+1)]
            double=[(s,s+a),(s+a+1,s+a+b+1)]
            u,v=(single,double) if side==0 else (double,single)
            r.extend(u);p.extend(v)
            s+=a+b+1+gap
    assert not target or s<=target
    return track(max(s,target),r),track(max(s,target),p)

def test_trade():
    for K,m,T in ((1,2,10),(11,1,10),(11,2,10),(12,3,11),(17,2,16),(4,5,10),(6,12,18)):
        a,b=trade(m,T)
        for mirror in (False,True):
            f,g=layout(a,K,mirror),layout(b,K,mirror)
            same(f,g,K,True,not mirror)
            h,r,p=stats(*f)
            l,r1,p1=stats(*g)
            assert (h,l)==((1+mirror)*m*m,1+mirror)
            assert (r,p)==(r1,p1)
            wake_stack=m*m*(2*T+2*m)
            c=confusion(f)
            assert c[1:]==(m*m*mirror,m*m,wake_stack*(1+mirror))
            if mirror:
                assert r==p==3*m*m
                assert Q(c[3],c[2]+c[3])==1-Q(1,4*T+4*m+1)
    a,b=trade(2,10)
    f,g=layout(a,11,True,960),layout(b,11,True,960)
    assert (stats(*f),stats(*g))==((8,12,12),(2,12,12))
    assert confusion(f)==(760,4,4,192) and f[0]!=g[0]
    print('Any K controls through K=T+1 pass, including m=1 zero-gap. At N=960: (760,4,4,192), 8/12 vs 2/12, both W=196.')

def local_count(fixture,U,minimum):
    M=4*U//3
    H=hist(*fixture,M+2); cumulative=0
    for width in range(3,M+3):
        for word,n in H[width-1].items():
            r=tuple(x>>1 for x in word);p=tuple(x&1 for x in word)
            cumulative+=n*sum(min(rr[0],pp[0])==1 and max(rr[1],pp[1])==width-1 and iou(rr,pp)>=half
                              for rr in bouts(r,minimum) for pp in bouts(p,minimum))
    return cumulative

def test_local():
    generator=Random(3611)
    fixture=(track(32,[(1,13)]),track(32,[(5,17)]))
    assert iou((1,13),(5,17))==half
    examples=[fixture]
    for _ in range(150):
        n=generator.randrange(6,40)
        examples.append(tuple((0,)+tuple(generator.getrandbits(1) for _ in range(n))+(0,) for _ in range(2)))
    for f in examples:
        U=max((b-a for x in f for a,b in bouts(x,1)),default=1)
        for minimum in (1,2,3,5):
            assert local_count(f,U,minimum)==stats(*f,minimum)[0]
    for tau in (Q(2,5),Q(49,100),Q(499,1000)):
        L=max(10,ceildiv(tau.numerator,tau.denominator-2*tau.numerator))
        assert all(tau<=iou((1,2*L+2),pair)<half for pair in ((1,L+1),(L+2,2*L+2)))
    print('Independent graph matching and reconstruction through floor(4U/3)+2: 604 exact controls, including inclusive equality and lower-threshold ambiguity.')

if __name__=='__main__':
    test_shift()
    test_trade()
    test_local()
