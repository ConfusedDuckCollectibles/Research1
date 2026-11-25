#!/usr/bin/env python3
"""Independent exact integer enumeration. Synthetic stock ranks, no primary data."""
from collections import Counter, deque
from itertools import permutations, product
from math import factorial
from random import Random

def neighbors(world, a, b):
    return sum(1<<i for i,j in enumerate(world) if i<a and j<b)

def min_union_tables():
    examined=0
    for n in range(5):
        types=list(product(range(1,n+1),repeat=2))
        num=1<<len(types)
        cut=[n+1]*num
        lower=[0]*num
        for m in range(1,num):
            bit=m&-m
            a,b=types[bit.bit_length()-1]
            lower[m]=max(lower[m^bit],a+b-n)
        for w in permutations(range(n)):
            edges=[neighbors(w,*q) for q in types]
            unions=[0]*num
            cut[0]=0
            for m in range(1,num):
                bit=m&-m
                unions[m]=unions[m^bit]|edges[bit.bit_length()-1]
                cut[m]=min(cut[m],unions[m].bit_count())
        assert cut==lower
        examined+=num
    print(f'Exact minimum neighbor union against every rank permutation for {examined} scopes, N=0..4.')

def matched_bits(n, edges):
    """Enumerate feasible matched request subsets, from item-by-item assignment.
    Independent of all cutoff formulae; state is only the served request bitset.
    """
    states={0}
    for i in range(n):
        eligible=[j for j,e in enumerate(edges) if e>>i&1]
        states|={m|1<<j for m in states for j in eligible if not m>>j&1}
    return states

def gf_rank(columns):
    basis={}
    for value in columns:
        v=value
        while v:
            key=v.bit_length()-1
            if key in basis: v^=basis[key]
            else:
                basis[key]=v
                break
    return len(basis)

def robust_formula(n, g, rank):
    return min(z+rank([i for i,s in enumerate(g) if s>z]) for z in range(n+1))

def virtual_mask(g,mask):
    s=sorted(g[j] for j in range(len(g)) if mask>>j&1)
    return all(x>=i for i,x in enumerate(s,1))

def audit_small():
    G=Random(6729841)
    checked=0
    for n in range(6):
        Worlds=tuple(permutations(range(n)))
        for _ in range(28):
            q=[(G.randrange(n+1),G.randrange(n+1)) for _ in range(G.randrange(8))]
            tags=[G.randrange(3) for _ in q]
            cols=[G.randrange(8) for _ in q]
            m=len(q)
            g=[max(0,a+b-n) for a,b in q]
            parts=[[i for i,t in enumerate(tags) if t==j] for j in range(3)]
            limits=tuple(product(*(range(len(P)+1) for P in parts)))
            feasible={v for v in limits}
            univ=set(range(1<<m))
            rcache={b:gf_rank(cols[i] for i in range(m) if b>>i&1) for b in univ}
            robust_plain=n+1; robust_G=n+1
            for w in Worlds:
                choices=matched_bits(n,[neighbors(w,*qj) for qj in q])
                univ &= choices
                counts={tuple(sum(mask>>i&1 for i in P) for P in parts) for mask in choices}
                feasible &= counts  # each family downward closed
                robust_plain=min(robust_plain,max(mask.bit_count() for mask in choices))
                robust_G=min(robust_G,max(rcache[b] for b in choices if rcache[b]==b.bit_count()))
                checked+=1
            assert univ=={b for b in range(1<<m) if virtual_mask(g,b)}
            assert robust_plain==robust_formula(n,g,len)
            assert robust_G==robust_formula(n,g,lambda ix:gf_rank(cols[i] for i in ix))
            pred={v for v in limits if all(sum(max(0,v[j]-sum(g[i]>z for i in P)) for j,P in enumerate(parts))<=z
                                           for z in range(n+1))}
            assert feasible==pred
    print(f'All-completion assignment, fixed robust request IDs, partition quotas and binary matroids: {checked} graph checks.')

def cutoff_menu():
    for n in range(1,9):
        families=Counter()
        for w in permutations(range(1,n+1)):
            for z in range(n):
                inside=[i for i,j in enumerate(w,1) if i+j<=n+z]
                assert len(inside)>=z
                if len(inside)==z:
                    families[z]+=1
                    assert inside==list(range(1,z+1))
                    assert w[n-1]==z+1
                    assert all(w[i-1]==n+z+1-i for i in range(z+1,n+1))
        assert families=={z:factorial(z) for z in range(n)}
    print('Unique reversed tails and disjoint minimum scenario families exhaustively checked through N=8.')

def maxflow(edges,s,t):
    adj={i:{} for e in edges for i in e[:2]}
    for u,v,c in edges:
        adj[u][v]=adj[u].get(v,0)+c
        adj[v][u]=adj[v].get(u,0)
    count=0
    while True:
        src={s:None}; q=deque([s])
        while q and t not in src:
            v=q.popleft()
            for w,cap in adj[v].items():
                if cap>0 and w not in src:
                    src[w]=v;q.append(w)
        if t not in src: return count
        v=t; f=10**12
        while src[v] is not None:
            f=min(f,adj[src[v]][v]); v=src[v]
        v=t
        while src[v] is not None:
            u=src[v]
            adj[u][v]-=f;adj[v][u]+=f
            v=u
        count+=f

def full_case(strong, demand, group):
    counts=(9,2,9)
    supply={(p,e): int(strong[e-1] if p==1 else counts[e-1]-strong[e-1])
            for p,e in product((0,1),range(1,4))}
    edges=[('s',('G',j),x) for j,x in enumerate(group)]
    for i,(p,e,n,j) in enumerate(demand):
        edges.append((('G',j),('q',i),n))
        edges.extend((('q',i),('u',p0,e0),20) for p0,e0 in supply if p0>=p and e0>=e)
    edges.extend((('u',p0,e0),'t',c) for (p0,e0),c in supply.items())
    return maxflow(edges,'s','t')

def fixture20():
    demand=((1,1,9,0),(1,2,1,0),(0,2,1,2),(0,3,9,1))
    stock=[x for x in product(range(10),range(3),range(10)) if sum(x)==10]
    vals=Counter(full_case(v,demand,(10,9,1)) for v in stock)
    assert len(stock)==28 and sorted(vals)==list(range(11,21))
    quotas={(i,j) for i in range(11) for j in range(10)
            if all(full_case(v,demand,(i,j,0))==i+j for v in stock)}
    assert quotas=={(i,j) for i in range(11) for j in range(10) if i+j<=10}
    assert full_case((9,1,0),demand,(10,9,1))==20
    assert full_case((0,1,9),demand,(10,9,1))==11
    print('N=20: all 28 latent tables, exact totals',sorted(vals.items()),'; quota pair region i+j<=10 verified.')

def adaptive_controls():
    # Two worlds illustrate why the objective is not (a) an ordinal blind plan,
    # (b) the union of independently worst request neighborhoods.
    direct=[(0,1),(1,0)]
    qs=[(1,2),(2,1),(1,1),(2,2)]
    sets=[matched_bits(2,[neighbors(w,*q) for q in qs]) for w in direct]
    assert 0b0011 not in sets[0] and 0b0011 in sets[1] # complementary rows
    assert 0b0100 in sets[0] and 0b0100 not in sets[1] # high/high singleton opposite
    assert 0b1001 in sets[0] and 0b1001 in sets[1]
    hereditary={0,1,2,4,8,3,12} # downward union of two 2-element faces, not matroid
    assert min(max(x.bit_count() for x in (M&hereditary)) for M in sets)==2
    assert max(x.bit_count() for x in (sets[0]&sets[1]&hereditary))==1
    for high in product(range(2),repeat=3):
        req=[(1<<high[i])&(1<<high[j]) for i,j in ((0,1),(0,2),(1,2))]
        assert max(m.bit_count() for m in matched_bits(2,req))==1
    print('Negative controls pass: no universal coupling, nonmatroid adaptive gap, and 3-axis obstruction.')

def both_margins():
    rng=Random(786)
    exam=0
    for _ in range(26):
        n=rng.randrange(1,5)
        margins=[]
        for _ in range(2):
            l=rng.randrange(1,4)
            # Request difficulties represented by eligibility counts, in decreasing difficulty order.
            margins.append(([rng.randrange(n+1) for _ in range(l)],[rng.randrange(n+1) for _ in range(l)]))
        demandWorlds=[list(set(tuple(zip(A,B)) for B in permutations(b))) for A,b in margins]
        L=[]
        for j,(A,B) in enumerate(margins):
            lj=[]
            for z in range(n+1):
                lower=max([0]+[sum(x>=a for x in A)+sum(y>=b for y in B)-len(A)
                               for a,b in product(A,B) if a+b>n+z])
                assert lower==min(sum(a+b>n+z for a,b in dw) for dw in demandWorlds[j])
                lj.append(lower)
            L.append(lj)
        quotas=set(product(*(range(len(A)+1) for A,B in margins)))
        robust=n+1
        for w1,w2 in product(*demandWorlds):
            q=w1+w2
            for stock in permutations(range(n)):
                choices=matched_bits(n,[neighbors(stock,*p) for p in q])
                robust=min(robust,max(x.bit_count() for x in choices))
                counts={( (x&((1<<len(w1))-1)).bit_count(), (x>>len(w1)).bit_count() ) for x in choices}
                quotas &= counts
                exam+=1
        assert robust==min(z+L[0][z]+L[1][z] for z in range(n+1))
        pred={tuple(v) for v in product(*(range(len(A)+1) for A,B in margins))
              if all(sum(max(0,v[j]-L[j][z]) for j in range(2))<=z for z in range(n+1))}
        assert quotas==pred
    print(f'Both stock and demand margins hidden: {exam} direct joint-world assignment checks.')

def no_holes_and_best():
    rng=Random(30079)
    cases=0
    for n in range(5):
        for _ in range(45):
            q=[(rng.randrange(n+1),rng.randrange(n+1)) for _ in range(rng.randrange(1,7))]
            worlds=list(permutations(range(n)))
            sizes=set()
            optimum=0
            # Check arbitrary hereditary families generated by a few random maximal sets.
            maxima=[rng.randrange(1<<len(q)) for _ in range(3)]
            for world in worlds:
                choices=matched_bits(n,[neighbors(world,*qj) for qj in q])
                sizes.add(max(x.bit_count() for x in choices if any(x&~c==0 for c in maxima)))
                optimum=max(optimum,max(x.bit_count() for x in choices))
            assert sorted(sizes)==list(range(min(sizes),max(sizes)+1))
            edges=[('s',('A',i),1) for i in range(n)]
            edges.extend((('B',j),'t',1) for j in range(n))
            edges += [('s','t',0)]
            for i,(a,b) in enumerate(q):
                edges.append((('Qin',i),('Qout',i),1))
                edges.extend((('A',j),('Qin',i),1) for j in range(a))
                edges.extend((('Qout',i),('B',j),1) for j in range(b))
            assert maxflow(edges,'s','t')==optimum
            cases+=1
    print(f'Best coupling path-flow and no integer holes under hereditary constraints: {cases} batches.')

if __name__=='__main__':
    min_union_tables()
    audit_small()
    cutoff_menu()
    fixture20()
    adaptive_controls()
    both_margins()
    no_holes_and_best()
