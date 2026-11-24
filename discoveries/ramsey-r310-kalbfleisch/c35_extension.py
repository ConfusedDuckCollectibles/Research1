#!/usr/bin/env python3
"""Exact three-terminal-type reduction for induced C35 in a hypothetical R(3,10,40)."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import argparse,json
HERE=Path(__file__).parent
N=35
D=(1,7,11,16)
SHIFTS={x%N for a in D for x in (a,-a)}
core=[sum(1<<((i+d)%N) for d in SHIFTS) for i in range(N)]
EXTERIORS={
    'star':tuple((u,4) for u in range(4)),
    'bipartite':tuple((u,v) for u in range(2) for v in range(2,5)),
    'pentagon':tuple(sorted(tuple(sorted((i,(i+1)%5))) for i in range(5))),
}
MASK=lambda n: (1<<n)-1

def inds(adj,universe):
    out={i:[] for i in range(len(adj)+1)}
    def go(S,U):
        if U==0:
            return
        while U:
            x=U&-U;U^=x;i=x.bit_length()-1
            T=S|x
            out[T.bit_count()].append(T)
            go(T,U&~adj[i])
    out[0].append(0)
    go(0,universe)
    return out

def mask_members(x):
    while x:
        y=x&-x;x^=y;yield y.bit_length()-1

def laws(g):
    assert all(v==sum(1<<j for j in range(len(g)) if g[j]>>i&1) and not (v>>i&1)
               for i,v in enumerate(g))
    assert all(not (g[i]&g[j]) for i,v in enumerate(g) for j in mask_members(v))
    for k,L in inds(g,MASK(len(g))).items():
        if L: big=k
    return big

def verify_baseline():
    m=(HERE/'exoo39.txt').read_text().split()
    assert len(m)==39 and all(len(s)==39 and set(s)<={'0','1'} for s in m)
    assert all(m[i][j]==m[j][i] and m[i][i]=='0' for i in range(39) for j in range(39))
    adj=[sum(1<<j for j,b in enumerate(r) if i!=j and b=='0') for i,r in enumerate(m)]
    assert laws(adj)==9
    assert sum(a.bit_count() for a in adj)==2*167
    assert laws(core)==8
    deletions=[]
    for B in combinations(range(39),4):
        if 35 not in B:
            continue # degree six cannot be inside an 8-regular induced graph
        T=sum(1<<i for i in B)
        if all((a&~T).bit_count()==8 for i,a in enumerate(adj) if i not in B):
            deletions.append(B)
    assert len(deletions)==1
    B=deletions[0]
    return {'edges':167,'alpha':9,'core_candidate_exterior':B,
            'exterior_edges':[(a,b) for a,b in combinations(B,2) if adj[a]>>b&1],
            'old_degree_histogram':dict(Counter(v.bit_count() for v in adj))}

def check_three_types():
    U=tuple(combinations(range(5),2))
    buckets=Counter()
    def graph(E):
        return [sum(1<<v if u==i else 1<<u if v==i else 0 for u,v in E) for i in range(5)]
    for h in range(1<<10):
        E=tuple(U[i] for i in mask_members(h))
        a=graph(E)
        if any(a[i]&a[j] for i,j in E):
            continue
        if any(not (a[i]&a[j]) for i,j in U if (i,j) not in E):
            continue
        deg=tuple(sorted(x.bit_count() for x in a))
        buckets[deg]+=1
    assert buckets=={(1,1,1,1,4):5,(2,2,2,3,3):10,(2,2,2,2,2):12}
    return {'-'.join(map(str,d)):n for d,n in buckets.items()}

class CNF:
    def __init__(self):
        self.cl=[]
        self.n=175
    def new(self):
        self.n+=1
        return self.n
    def add(self,*cl):
        if any(type(x) is bool and x is True for x in cl):
            return
        lits=dict.fromkeys(x for x in cl if type(x) is int)
        if any(-x in lits for x in lits):
            return
        self.cl.append(tuple(lits))
    @staticmethod
    def inv(x):
        return not x if type(x) is bool else -x
    def bounds(self,L,lo,hi):
        P={0:True}
        for i,x in enumerate(L,1):
            Q={0:True}
            for j in range(1,min(i,hi+1)+1):
                A=P.get(j,False); B=P.get(j-1,False)
                s=self.new()
                Q[j]=s
                # s <-> A or (B and x)
                self.add(self.inv(A),s)
                self.add(self.inv(B),-x,s)
                self.add(-s,A,B)
                self.add(-s,A,x)
            P=Q
        self.add(P.get(lo,False))
        self.add(self.inv(P.get(hi+1,False)))
    def output(self,path):
        with path.open('w') as f:
            f.write(f'p cnf {self.n} {len(self.cl)}\n')
            for c in self.cl:
                f.write(' '.join(map(str,c))+' 0\n')

def term(u,h):
    return 5*u+h+1

def build(F, I, anchor=True, cross_saturate=True):
    T=CNF()
    out=[sum(1<<v if u==i else 1<<u if v==i else 0 for u,v in F) for i in range(5)]
    J=inds(out,31)
    counts={}
    for x in range(N):
        for u,v in combinations(range(5),2):
            T.add(-term(x,u),-term(x,v))
    for x in range(N):
        for y in mask_members(core[x]):
            if y>x:
                for j in range(5):
                    T.add(-term(x,j),-term(y,j))
    for a,A in J.items():
        k=10-a
        if not A or k>8:continue
        cnt=0
        for i in I[k]:
            xs=tuple(mask_members(i))
            for b in A:
                ys=tuple(mask_members(b))
                T.add(*(term(x,h) for x in xs for h in ys))
                cnt+=1
        counts[a]=cnt
    for h in range(5):
        d=out[h].bit_count()
        # Minimum degree >=4 follows from the published R(3,9)=36.
        T.bounds([term(x,h) for x in range(N)],max(0,4-d),9-d)
    if cross_saturate:
        for x in range(N):
            for h in range(5):
                T.add(term(x,h),*(term(y,h) for y in mask_members(core[x])),
                    *(term(x,j) for j in mask_members(out[h])))
    if anchor:
        assert out[0].bit_count()<4 # hence at least one core link by the min degree bound
        T.add(term(0,0))
    return T,counts

def write_cubes(skip_exoo=False):
    I=inds(core,MASK(N))
    assert max(i for i,s in I.items() if s)==8
    out={'baselines':{} if skip_exoo else verify_baseline(),'labeled_allowed_exteriors':check_three_types(),
         'core_independence_counts':{i:len(L) for i,L in I.items() if L},'types':{}}
    for name,F in EXTERIORS.items():
        c,a=build(F,I)
        path=HERE/f'{name}.cnf'
        c.output(path)
        out['types'][name]={'vars':c.n,'clauses':len(c.cl),'cover_counts':a,'exterior_edges':F}
    (HERE/'model-audit.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--skip-exoo-baseline',action='store_true',help='No third-party data required.')
    write_cubes(p.parse_args().skip_exoo_baseline)
if __name__=='__main__':
    main()
