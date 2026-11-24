"""Exact certificate producer for the f_Z(11)=526 equality test.

All arithmetic is integer. Rows are ruled out by an exact-cover hole,
possibly after forced single coverage and consequent collision blocks.
A separate direct, non-bitset checker validates each certificate.
"""
from collections import Counter
from itertools import combinations
import argparse
import json

M=526
N=11

if not __debug__:
    raise RuntimeError('Run without -O, assertions are checks.')


def sums(a):
    h=[0]
    for w in a:
        h.extend([w+x for x in h])
    return h


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


def candidates(profile):
    _,_,P,C=profile
    # Each elementary lower bound uses C_j(p)-C_(j-1)(p)<C_(j-1)(q).
    triples=[(7,1249,800),(8,1340,835),(9,1365,845),(10,1372,848)]
    assert all(p in P and q in P and C[p][j]-C[p][j-1]<C[q][j-1] for j,p,q in triples)
    assert [p-q+1 for j,p,q in triples]==[450,506,521,525]
    for c in range(521,525):
        for b in range(506,c):
            for a in range(450,b):
                yield (a,b,c,525,526)


def setup(t,profile,n=N):
    band_lo,band_hi,P,C=profile
    j0=n-len(t)
    h=[sums(t[:k]) for k in range(len(t)+1)]
    E=h[-1]
    if len(E)!=len(set(E)):
        return
    V=1<<n
    tau=V+1-sum(t)
    width=j0*t[0]-j0*(j0+1)//2
    L=(tau-width+1)//2
    R=(tau+width)//2
    # B contains disallowed indices, U linear order L,...,R.
    B=0
    def clip_mask(l,r):
        l=max(L,l)
        r=min(R,r)
        return 0 if l>r else ((1<<(r-l+1))-1)<<(l-L)
    def cut(a,b,H):
        nonlocal B
        for e in H:
            left,right=a+1-e,b-e
            B|=clip_mask(left,right)|clip_mask(tau-right,tau-left)
    for j in range(j0,n+1):
        Hj=h[j-j0]
        points=[(L-1,0),*((p,C[p][j]) for p in P),(R+sum(t),1<<j)]
        assert points==sorted(points)
        nonj=set(E)-set(Hj)
        for (a,c),(b,d) in zip(points,points[1:]):
            if d==c:
                cut(a,b,Hj)
            if d-c==b-a:
                cut(a,b,nonj)
    # Missing midpoint is implied by complement distinctness of low subsets.
    if tau%2==0:
        B|=clip_mask(tau//2,tau//2)
    D={a-b for a in E for b in E if a!=b}
    for d in D:
        # If x and its distinct reflection differ by a top subset difference,
        # selecting the reflection orbit creates a sum collision.
        if (tau-d)%2==0:
            B|=clip_mask((tau-d)//2,(tau-d)//2)
    W=((1<<(R-L+1))-1)
    A=W^B
    # Small dimensions may have low positions to the right of the band start;
    # clamp band mask to nonnegative shifted positions.
    band=((1<<(band_hi-max(band_lo,L)+1))-1)<<(max(band_lo,L)-L)
    S=sum(1<<(d+sum(t)) for d in D)
    return L,R,A,tau,sorted(E),band,S


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


def prove(t,profile,n=N,branch_limit=500):
    tmp=setup(t,profile,n)
    if tmp is None:
        return ['D'],0
    L,R,A,tau,E,band,S=tmp
    rad=sum(t)
    width=R-L+1
    target=1<<(n-len(t))
    explored=0
    def visit(A,chosen):
        nonlocal explored
        steps=[]
        while True:
            if chosen.bit_count()>target or A.bit_count()<target:
                return [steps,'K']
            hole,unit=cover(E,A,band)
            if hole:
                return [steps,'H',pick(hole)+L]
            todo=0
            while unit:
                k=pick(unit)
                unit&=unit-1
                a=next(k-e for e in E if k>=e and (A>>(k-e))&1)
                if not (chosen>>a)&1:
                    todo=1
                    break
            if not todo:
                break
            z=tau-2*L-a
            # By symmetry z should still be active.
            assert 0<=a<width and 0<=z<width and (A>>z)&1
            assert a!=z
            C=(1<<a)|(1<<z)
            bad=move(a,S,rad)|move(z,S,rad)
            assert not bad&C and not bad&chosen
            A&=~bad
            chosen|=C
            steps.append(k+L)
        # For first pilot feasibility: branch on a smallest uncovered target
        # row. Branches are proved explicitly, never treated as absent.
        if not branch_limit or explored>=branch_limit:
            return None
        todo=band
        for e in E:
            todo&=~(chosen<<e)
        if not todo:
            return None
        opts=None
        while todo:
            k=pick(todo)
            todo&=todo-1
            cur=[k-e for e in E if k>=e and (A>>(k-e))&1]
            if opts is None or len(cur)<len(opts):
                q=k
                opts=cur
                if len(cur)==2:
                    break
        # Every branch from the exact cover row uses a distinct orbit.
        assert opts and len(opts)>=2
        kids=[]
        for a in opts:
            if explored>=branch_limit:
                return None
            explored+=1
            z=tau-2*L-a
            chosen2=chosen|(1<<a)|(1<<z)
            bad=move(a,S,rad)|move(z,S,rad)
            assert not bad&chosen2
            v=visit(A&~bad,chosen2)
            if v is None:
                return None
            kids.append(v)
        return [steps,'B',q+L,kids]
    p=visit(A,0)
    return p,explored


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--limit',type=int,default=500)
    ap.add_argument('--prefix',type=int,default=5000)
    ap.add_argument('--out',default='occupancy-full.cert.jsonl')
    a=ap.parse_args()
    profile=rank_profile(N,M)
    assert profile[:2]==(674,1375)
    T=list(candidates(profile))
    h=Counter()
    max_ex=0
    with open(a.out,'w') as f:
        for i,t in enumerate(T[:a.prefix]):
            cert,e=prove(t,profile,branch_limit=a.limit)
            max_ex=max(max_ex,e)
            if cert is None:
                f.write(json.dumps([t,None,e],separators=(',',':'))+'\n')
                h['unresolved']+=1
            else:
                f.write(json.dumps([t,cert,e],separators=(',',':'))+'\n')
                h['duplicate' if cert==['D'] else 'unit-only' if e==0 else 'branched']+=1
            if i%200==199:
                print(i+1,dict(h),'max branch nodes',max_ex,flush=True)
    print('FINAL',sum(h.values()),dict(h),'max branch nodes',max_ex)

if __name__=='__main__':
    main()
