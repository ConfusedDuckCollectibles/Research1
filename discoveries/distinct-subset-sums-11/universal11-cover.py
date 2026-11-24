"""Build exact certificates for the dimension-independent cover obstruction."""
from math import comb
from itertools import combinations
from collections import Counter
import argparse
import json

SCALE=2048

if not __debug__:
    raise RuntimeError('Run without -O.')


def rank_profile(n,m):
    """Build simplicial cube prefix profiles without external data."""
    V=1<<n
    B=0
    C=[(0,)*(n+1)]
    cut=[]
    O=[sum(1<<i for i in S) for k in range(n+1) for S in combinations(range(n),k)]
    for q,S in enumerate(O,1):
        B|=(1<<S)|sum(1<<(S^(1<<i)) for i in range(n))
        C.append(tuple(c+(S<(1<<j)) for j,c in enumerate(C[-1])))
        size=B.bit_count()-q
        if size>m:
            raise ValueError('not an equality profile')
        if size==m:
            cut.append(q)
    assert cut and cut[-1]-cut[0]<m
    P=sorted(set(cut+[x+m for x in cut]))
    L,R=P[0],P[-1]
    assert L+R==V+1
    return L,R,P,C


def move(x,S,rad):
    return S<< (x-rad) if x>=rad else S>>(rad-x)


def cover(E,A,band):
    seen=multi=0
    for e in E:
        p=A<<e
        multi|=seen&p
        seen|=p
    return band&~seen,band&seen&~multi


def pick(x):
    return (x&-x).bit_length()-1

def local_bans():
    # Suffix template (5 transversal pairs followed by the maximum vertex);
    # grade masks are read using the last 5 bits.  The reflected interval
    # [p,q) prescribes counts for the conjugate grade (tail superset).
    n=11; H=sum(comb(i,i//2) for i in range(n))
    lo,hi,P,C=rank_profile(n,H)
    A=[sum(1<<i for i in s)>>(n-5) for k in range(n+1) for s in combinations(range(n),k)]
    P0=P[:len(P)//2]
    assert all(x-lo<=175 for x in P0)
    ans=[]
    for h in range(2):
        table=[0]*176
        for p,q in zip(P0,P0[1:]):
            p0,q0=p+h*H,q+h*H
            for k in range(6):
                g=(1<<k)-1
                hits=sum(x&~g==0 for x in A[p0:q0])
                mirrorhits=sum(g|x==31 for x in A[p0-1:q0-1])
                D=0; R=0
                for x in range(32):
                    if (hits==0 and x&~g==0) or (hits==q-p and x&~g!=0):
                        D|=1<<x
                    if (mirrorhits==0 and g|x==31) or (mirrorhits==q-p and g|x!=31):
                        R|=1<<x
                for u in range(p-lo+1,q-lo+1):
                    table[u]|=D
                for u in range(p-lo,q-lo):
                    table[u]|=R
        ans.append(table)
    return [[sum(1<<q for q,v in enumerate(tab) if v>>i&1)
             for i in range(32)] for tab in ans]


def prove(t,T,limit,K,upper=True):
    # e from high tail weights SCALE minus deficit (0 to 102).
    E=[sum(SCALE-d for i,d in enumerate(t) if x>>i&1) for x in range(32)]
    if len(E)!=len(set(E)):return ['D'],0
    r=6*SCALE
    B=((1<<(K+1))-1)<<r
    if upper:
        B|=((1<<176)-1)<<(r+SCALE)
    A=0
    for e in E:A|=B>>e
    for i,e in enumerate(E):
        for k,tab in enumerate(T):
            F=tab[i]
            A&=~(F<<(r+k*SCALE-e))
    D={p-q for p in E for q in E if p!=q}
    D=sum(1<<(q+r) for q in D)
    exp=0
    def visit(A,X):
        nonlocal exp
        steps=[]
        while True:
            z,unit=cover(E,A,B)
            if z:return [steps,'H',pick(z)-r]
            todo=0
            while unit:
                q=pick(unit);unit&=unit-1
                a=next(q-e for e in E if (A>>(q-e))&1)
                if not (X>>a)&1:
                    todo=1;break
            if not todo:break
            rem=move(a,D,r)
            assert not rem&X
            A&=~rem; X|=1<<a
            steps.append(q-r)
        todo=B
        for e in E:todo&=~(X<<e)
        if not todo:return None
        opts=None
        while todo:
            q=pick(todo);todo&=todo-1
            O=[q-e for e in E if (A>>(q-e))&1]
            if opts is None or len(O)<len(opts):
                opts=O
                mark=q
                if len(O)==2:break
        pf=[]
        for a in opts:
            if exp>=limit:return None
            exp+=1
            rem=move(a,D,r)
            assert not rem&X
            k=visit(A&~rem,X|(1<<a))
            if k is None:return None
            pf.append(k)
        return [steps,'B',mark-r,pf]
    p=visit(A,0)
    return p,exp

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='universal11.cert.jsonl')
    args=ap.parse_args()
    K,cap,up=238,100,True
    T=local_bans()
    tails=[(a,b,c,1,0) for c in range(5,1,-1) for b in range(20,c,-1) for a in range(76,b,-1)]
    counts=Counter(); maxi=0
    with open(args.out,'w') as f:
        for i,t in enumerate(tails,1):
            p,e=prove(t,T,cap,K,up)
            f.write(json.dumps([t,p,e],separators=(',',':'))+'\n')
            counts['unknown' if p is None else 'duplicate' if p==['D'] else 'branch' if e else 'unit']+=1
            maxi=max(maxi,e)
            if i%250==0:print(i,counts,maxi)
    assert counts=={'unit':3487,'duplicate':455,'branch':268}
    print('FINISHED',dict(counts),'largest tree edges',maxi)
