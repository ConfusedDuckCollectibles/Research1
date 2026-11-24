"""Direct coordinate-pair verifier; does not import the bitset generator."""
from collections import Counter
from itertools import combinations, product
from math import comb
import argparse
import json
from pathlib import Path

if not __debug__:
    raise RuntimeError('Run without -O.')

Point=tuple[int,int]
W=2048
K=238
MAXUP=175


def plus(p,q):
    return (p[0]+q[0],p[1]+q[1])


def minus(p,q):
    return (p[0]-q[0],p[1]-q[1])


def template():
    n=11
    H=sum(comb(i,i//2) for i in range(n)); assert H==526
    # Build the suffix transversals directly; no boundary generator used.
    sets=list(product(*(range(j,j+2) for j in range(1,11,2))))
    assert len(sets)==32
    lo=tuple(sets[0]); hi=tuple(sets[-1])
    order=[s for k in range(n+1) for s in combinations(range(1,n+1),k)]
    rank={s:i+1 for i,s in enumerate(order)}
    P=sorted(rank[s] for s in sets)
    L=rank[lo]; U=rank[hi]; assert (L,U)==(674,849) and U-L==MAXUP
    assert all(rank[s+(11,)]==rank[s]+H for s in sets)
    # Independently verify all four required grade-count rises.
    cmp=[(7,723,800),(8,814,835),(9,839,845),(10,846,848)]
    assert all(q in P and p in P and q<U for j,p,q in cmp)
    for j,p,q in cmp:
        assert any(max(s)<j for s in order[p:q])
    assert [q-p-1 for j,p,q in cmp]==[76,20,5,1]
    masks=[sum(1<<(i-7) for i in s if i>=7) for s in order]
    bans=[set() for _ in range(32)]
    for h in (0,1):
        for a,b in zip(P,P[1:]):
            p,q=a+h*H,b+h*H
            nsites=q-p
            assert q-p>0
            # Use explicit images of the two different half-open intervals.
            direct=[masks[s-1] for s in range(p+1,q+1)]
            mirror=[31^masks[s-1] for s in range(p,q)]
            for k in range(6):
                prefix=set(range(1<<k))
                d=sum(x in prefix for x in direct)
                v=sum(x in prefix for x in mirror)
                for z in range(32):
                    if d==0 and z in prefix or d==nsites and z not in prefix:
                        bans[z].update((h,s-L) for s in range(a+1,b+1))
                    if v==0 and (31^z) in prefix or v==nsites and (31^z) not in prefix:
                        bans[z].update((h,s-L) for s in range(a,b))
    requ={(0,q) for q in range(K+1)}|{(1,q) for q in range(MAXUP+1)}
    return bans,requ


class Instance:
    def __init__(self,shifts,required,bans):
        self.shifts=tuple(shifts) # certificate children follow upper-mask index
        assert len(self.shifts)==len(set(self.shifts))
        self.required=required
        candidates={minus(q,e) for q in required for e in shifts}
        for e,F in zip(shifts,bans,strict=True):
            candidates.difference_update(minus(q,e) for q in F)
        self.candidates=candidates
        self.diff={minus(e,f) for e in shifts for f in shifts if e!=f}

    def check(self,tree,allowed=None,chosen=frozenset()):
        if allowed is None:
            allowed=self.candidates
        def choices(v):
            q=divmod(v,W); assert q in self.required
            return [minus(q,e) for e in self.shifts if minus(q,e) in allowed]
        def select(x):
            assert x in allowed and x not in chosen
            A=allowed-{plus(x,d) for d in self.diff}
            X=chosen|{x}
            assert X<=A
            return A,X
        assert isinstance(tree,list) and len(tree)>=2
        forced,action,*rest=tree
        assert isinstance(forced,list)
        for f in forced:
            cand=choices(f)
            assert len(cand)==1
            allowed,chosen=select(cand[0])
        if action=='H':
            assert len(rest)==1 and not choices(rest[0])
            return 0,1,len(forced)
        assert action=='B' and len(rest)==2
        cand=choices(rest[0]); down=rest[1]
        assert len(cand)>=2 and len(down)==len(cand) and not chosen.intersection(cand)
        edges,leaves,units=len(cand),0,len(forced)
        for x,T in zip(cand,down,strict=True):
            e,l,u=self.check(T,*select(x))
            edges+=e;leaves+=l;units+=u
        return edges,leaves,units


def controls():
    # A feasible translated cover: {0,3}+{0,1} covers {0,1,3,4}.
    good=Instance([(0,0),(0,1)],{(0,i) for i in (0,1,3,4)},[set(),set()])
    witness={(0,0),(0,3)}
    assert witness<=good.candidates
    assert {plus(x,e) for x in witness for e in good.shifts}==good.required
    for bad in ([[],'H',0],[[],'B',1,[[[],'H',4]]]):
        try:
            good.check(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('accepted a false certificate')
    # A one-row certificate with both shifts explicitly forbidden.
    no=Instance([(0,0),(0,1)],{(0,0)},[{(0,0)},{(0,0)}])
    assert no.check([[],'H',0])==(0,1,0)
    print('Feasible, infeasible and malformed-certificate controls passed.')


def main(path):
    controls()
    bans,required=template()
    tails=[(a,b,c,1,0) for c in range(5,1,-1) for b in range(20,c,-1) for a in range(76,b,-1)]
    lines=Path(path).read_text().splitlines()
    assert len(tails)==len(lines)==4210
    total=Counter(); audit=Counter()
    for t,s in zip(tails,lines,strict=True):
        tail,cert,nodes=json.loads(s)
        assert tuple(tail)==t
        shifts=[(m.bit_count(),-sum(t[j] for j in range(5) if m>>j&1)) for m in range(32)]
        if len(set(shifts))<32:
            assert cert==['D'] and nodes==0
            total['duplicates']+=1
            continue
        edges,leaves,units=Instance(shifts,required,bans).check(cert)
        assert edges==nodes
        total['branched' if nodes else 'unit-only']+=1
        audit['edges']+=edges
        audit['leaves']+=leaves
        audit['units']+=units
        audit['max edges']=max(audit['max edges'],edges)
    assert total=={'duplicates':455,'unit-only':3487,'branched':268}
    print('UNIVERSAL CERTIFICATE CHECKED',dict(total),dict(audit))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('certificate',nargs='?',default='universal11.cert.jsonl')
    main(p.parse_args().certificate)
