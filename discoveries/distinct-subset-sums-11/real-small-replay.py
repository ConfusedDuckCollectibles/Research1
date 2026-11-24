"""Independent exact replay for the finite real n=9 and n=10 cover cases."""
from collections import Counter
from itertools import combinations
from math import comb
import argparse

if not __debug__:
    raise RuntimeError('Run without -O.')


def see(n):
    m=sum(comb(j,j//2) for j in range(n))
    Q=1<<n
    k=4 if n==9 else 5
    d=k-1
    low=n-k
    order=[sum(1<<i for i in s)
           for j in range(n+1) for s in combinations(range(n),j)]
    rank={s:i for i,s in enumerate(order,1)}
    N=0
    P=[]
    for p,s in enumerate(order,1):
        N|=(1<<s)|sum(1<<(s^(1<<i)) for i in range(n))
        assert N.bit_count()-p<=m
        if N.bit_count()-p==m:
            assert not s>>(n-1)
            assert rank[s|(1<<(n-1))]==p+m
            P.append(p)
    L=P[0]; W=P[-1]-L
    assert L+P[-1]+m==Q+1
    assert (n,m,W,len(P),L) in ((9,148,49,16,158),(10,274,175,32,288))
    masks=[x>>low for x in order]
    gap_witness=[]
    for a in range(d):
        j=low+a+1
        _,p,q=min((q-p,p,q) for p in P for q in P if q>p
                      and any(s>>(j-1)==0 for s in order[p:q]))
        gap_witness.append((j,p,q,q-p-1))
    D=[x[-1] for x in gap_witness]
    assert D==([20,5,1] if n==9 else [76,20,5,1])
    print('profiles',n,m,L,W,len(P),'bounds, via rank rises',gap_witness)
    p=1<<k
    F=[0]*p
    for shift in (0,m):
        for a,b in zip(P,P[1:]):
            q,r=a+shift,b+shift
            width=r-q
            lo=((1<<width)-1)<<(q-L)
            hi=lo<<1
            for g in range(k+1):
                d=sum(e<(1<<g) for e in masks[q:r])
                c=sum(((p-1)^e)<(1<<g) for e in masks[q-1:r-1])
                for e in range(p):
                    eligible=e<(1<<g)
                    if (d==0 and eligible) or (d==width and not eligible):
                        F[e]|=hi
                    eligible=((p-1)^e)<(1<<g)
                    if (c==0 and eligible) or (c==width and not eligible):
                        F[e]|=lo
    assert 0<=63<m
    tests=list(range(W+64))+list(range(m,m+W+1))
    assert len(tests)==len(set(tests)) and max(tests)==m+W
    return m,D,F,tests


def least(b):
    return (b&-b).bit_length()-1


def left(y,k):
    return y<<k if k>=0 else y>>(-k)


def prove(m,t,F,tests):
    K=len(t)
    E=[sum(m-t[i] for i in range(K) if mask>>i&1) for mask in range(1<<K)]
    if len(E)!=len(set(E)):
        return 'duplicate',0
    z=(K+1)*m
    U=sum(1<<(p+z) for p in tests)
    A=0
    for e in E:
        A|=U>>e
    for e,f in zip(E,F,strict=True):
        A &= ~(f<<(z-e))
    deltas={a-b for a in E for b in E if a!=b}
    conflicts=sum(1<<(d+z) for d in deltas)
    D=0
    def rec(A,T):
        nonlocal D
        while True:
            extra=U&~T
            if not extra:return True
            B=0; C=0
            for e in E:
                X=(A<<e)&extra
                C|=X&B
                B|=X
            if extra&~B:
                return False
            lone=B&~C
            if not lone:
                break
            p=least(lone)
            x,=tuple(p-e for e in E if (A>>(p-e))&1)
            T |= sum(1<<(x+e) for e in E)
            A &= ~(left(conflicts,x-z)|(1<<x))
        best=None
        while extra:
            p=least(extra)
            extra&=extra-1
            found=tuple(p-e for e in E if (A>>(p-e))&1)
            if best is None or len(found)<len(best):
                best=found
                if len(found)==2:break
        assert best and len(best)>1
        for x in best:
            D+=1
            if D>100:
                raise RuntimeError('cap, no conclusion')
            if rec(A&~(left(conflicts,x-z)|(1<<x)),T|sum(1<<(x+e) for e in E)):
                return True
        return False
    if rec(A,0):
        return 'witness',D
    return ('branch' if D else 'unit'),D


def controls():
    # Upper sums from (1,2) are 0,1,2,3 and a single lower point 0
    # is feasible. Banning all four representations of rank 0 is not.
    assert prove(2,(1,0),[0]*4,list(range(4)))[0]=='witness'
    assert prove(2,(1,0),[1]*4,list(range(4)))==('unit',0)
    print('Exact-cover controls passed.')


def main(n):
    m,b,F,Q=see(n)
    if n==9:
        T=[(a,c,1,0) for c in range(5,1,-1) for a in range(20,c,-1)]
    else:
        T=[(a,b,c,1,0) for c in range(5,1,-1) for b in range(20,c,-1) for a in range(76,b,-1)]
    result=Counter(); nodes=peak=0
    for t in T:
        label,d=prove(m,t,F,Q)
        result[label]+=1
        nodes+=d
        peak=max(peak,d)
    expected={9:{'duplicate':4,'branch':3,'unit':59},10:{'duplicate':455,'branch':262,'unit':3493}}
    assert result==expected[n]
    print('CHECKED',n,len(T),dict(result),'children',nodes,'peak',peak)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('n',nargs='?',type=int,choices=(9,10))
    opt=ap.parse_args()
    controls()
    for n in ((9,10) if opt.n is None else (opt.n,)):
        main(n)
