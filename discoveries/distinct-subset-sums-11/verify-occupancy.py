"""Independent direct-set verifier for the f_Z(11)>=527 certificate.

No import of the producer or the previous rank-tree code. Requires only
Python standard library. A proof node covers all supported choices in
its stated exact-cover row or ends at an explicit contradiction.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import argparse
import json


if not __debug__:
    raise RuntimeError('Run the verifier without -O, its assertions are checks.')


def simplicial(n, m):
    """Compute count tables using tuples and ordinary cube sets."""
    order = [x for k in range(n+1) for x in combinations(range(1,n+1),k)]
    vertices=[frozenset(x) for x in order]
    included = set()
    boundary = set()
    saturation = []
    for p,x in enumerate(order,1):
        s=vertices[p-1]
        included.add(s)
        boundary.discard(s)
        for j in range(1,n+1):
            other=s ^ {j}
            if other not in included:
                boundary.add(other)
        assert len(boundary)<=m
        if len(boundary)==m:
            assert included|boundary==set(vertices[:p+m])
            saturation.append(p)
    assert saturation
    band=(saturation[0],saturation[-1]+m)
    assert max(saturation)-min(saturation)<m
    assert sum(band)==(1<<n)+1
    P=sorted(set(saturation)|{p+m for p in saturation})
    C={p:[sum(not x or max(x)<=j for x in order[:p]) for j in range(n+1)]
       for p in P}
    return band,P,C


def upper(t):
    return [(sum(t[i] for i in range(len(t)) if mask & (1<<i)), mask.bit_length())
            for mask in range(1<<len(t))]


class State:
    def __init__(self, n, t, profile):
        self.t = t
        self.low = n-len(t)
        self.target = 1<<self.low
        self.central, P, C = profile
        self.tau = (1<<n)+1-sum(t)
        raw = upper(t)
        self.E = sorted(x for x,_ in raw)
        assert len(set(self.E))==len(self.E)
        # Since the low factors are distinct integers less than t[0], their
        # total is at most low*t[0]-low*(low+1)//2.
        width = sum(t[0]-i for i in range(1,self.low+1))
        a=(self.tau-width+1)//2
        b=(self.tau+width)//2
        assert self.tau-a==b
        F=[set() for _ in range(len(t)+1)]
        for j in range(self.low,n+1):
            edges=[(a-1,0)]+[(p,C[p][j]) for p in P]+[(b+sum(t),1<<j)]
            assert edges==sorted(edges)
            for (p,c),(q,d) in zip(edges,edges[1:]):
                if d==c:
                    for g in range(j-self.low+1):
                        F[g].update(range(p+1,q+1))
                if d-c==q-p:
                    for g in range(j-self.low+1,len(t)+1):
                        F[g].update(range(p+1,q+1))
        impossible=set()
        for e,g in raw:
            impossible.update(z-e for z in F[g])
        window=set(range(a,b+1))
        window-=impossible
        self.allowed=window & {self.tau-q for q in window}
        self.allowed={q for q in self.allowed
                      if q*2 != self.tau
                      and len({r+e for r in (q,self.tau-q) for e in self.E})==2*len(self.E)}
        self.forced=set()

    def choices(self,q):
        lo,hi=self.central
        assert lo<=q<=hi
        return [q-e for e in self.E if q-e in self.allowed]

    def select(self,x):
        assert x not in self.forced
        couple={x,self.tau-x}
        assert len(couple)==2 and couple<=self.allowed
        fixed={q+e for q in couple for e in self.E}
        assert len(fixed)==len(couple)*len(self.E)
        assert not fixed.intersection(q+e for q in self.forced for e in self.E)
        blocked={w-e for w in fixed for e in self.E}-couple
        allowed=self.allowed-blocked
        forced=self.forced|couple
        assert forced<=allowed
        s=object.__new__(State)
        s.__dict__.update(self.__dict__)
        s.allowed,s.forced=allowed,forced
        return s


def proof(s,tree):
    assert isinstance(tree,list) and len(tree)>=2
    steps,*v=tree
    assert isinstance(steps,list)
    for q in steps:
        k=s.choices(q)
        assert len(k)==1
        s=s.select(k[0])
    if v[0]=='H':
        assert len(v)==2 and not s.choices(v[1])
        return 0,1,len(steps)
    if v[0]=='K':
        assert len(v)==1
        assert len(s.forced)>s.target or len(s.allowed)<s.target
        return 0,1,len(steps)
    assert len(v)==3 and v[0]=='B'
    k=s.choices(v[1])
    assert len(k)>=2 and not (set(k)&s.forced)
    assert isinstance(v[2],list) and len(k)==len(v[2])
    nodes=len(k)
    leaves=0
    units=len(steps)
    for x,sub in zip(k,v[2],strict=True):
        n,l,u=proof(s.select(x),sub)
        nodes+=n
        leaves+=l
        units+=u
    return nodes,leaves,units


def controls():
    # Cases of equality in the bandwidth lower bound. Check the actual
    # witness's full low occupancy survives the relaxation, and that a
    # deliberately false missing-rank claim cannot be verified.
    for a in ((2,3,4),(3,5,6,7),(6,9,11,12,13)):
        n=len(a)
        feature=simplicial(n,a[-1])
        A=sorted(x for x,_ in upper(a))
        assert len(A)==len(set(A))
        (lo,hi),P,C=feature
        u=A[lo-1]-lo
        assert all(A[q-1]==q+u for q in range(lo,hi+1))
        low=a[:-2]
        high=a[-2:]
        T={x-u for x,_ in upper(low)}
        S=State(n,high,feature)
        assert T<=S.allowed
        for p in P:
            for j in range(n+1):
                assert sum(x<=p+u for x,_ in upper(a[:j]))==C[p][j]
        for q in range(lo,hi+1):
            assert len(set(S.choices(q))&T)==1
        try:
            proof(S,[[],'H',lo])
        except AssertionError:
            pass
        else:
            raise AssertionError('accepted an invalid missing-rank proof')
    print('Three positive equality controls passed; false certificates rejected.')


def main(p):
    controls()
    profile=simplicial(11,526)
    band,P,C=profile
    assert band==(674,1375) and len(P)==64
    def lower(j,p,q):
        assert C[p][j]-C[p][j-1]<C[q][j-1]
        return p-q+1
    a_min=lower(7,1249,800)
    b_min=lower(8,1340,835)
    c_min=lower(9,1365,845)
    d_min=lower(10,1372,848)
    assert (a_min,b_min,c_min,d_min)==(450,506,521,525)
    tails=[(a,b,c,525,526) for c in range(c_min,d_min)
           for b in range(b_min,c) for a in range(a_min,b)]
    lines=Path(p).read_text().splitlines()
    assert len(tails)==4210==len(lines)
    out=Counter()
    marked=Counter()
    for t,line in zip(tails,lines,strict=True):
        t_cert,cert,e=json.loads(line)
        assert tuple(t_cert)==t
        if len({x for x,_ in upper(t)})<32:
            assert cert==['D'] and e==0
            out['duplicate']+=1
            continue
        assert cert is not None
        n,l,u=proof(State(11,t,profile),cert)
        assert n==e
        out['branched' if n else 'unit-only']+=1
        marked['branch nodes']+=n
        marked['leaves']+=l
        marked['unit selections']+=u
        marked['largest tree'] = max(marked['largest tree'],n)
    assert out=={'unit-only':3639,'duplicate':455,'branched':116}
    print('ALL CHECKED',dict(out),dict(marked))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('certificate',nargs='?',default='occupancy-full.cert.jsonl')
    a=ap.parse_args()
    main(a.certificate)
