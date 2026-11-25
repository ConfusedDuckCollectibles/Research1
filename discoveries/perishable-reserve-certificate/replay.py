#!/usr/bin/env python3
"""Exact small-integer replay. Only synthetic data."""
from collections import Counter
from itertools import permutations, product
from random import Random

def gscore(r,s,D):
    return [max(0,sum(a>=p for a in r)+sum(b>=t for b in s)-len(r)) for p,t in D]

def perm(z,n):
    return list(range(z))+list(range(n-1,z-1,-1))

def canonical(r,s,z):
    return list(zip(sorted(r,reverse=True),[sorted(s,reverse=True)[i] for i in perm(z,len(r))]))

def formula(gs,N,rank=None):
    D=range(len(gs))
    def v(z):
        E=tuple(i for i in D if gs[i]>z)
        return z+(len(E) if rank is None else rank(E))
    return min(v(z) for z in range(N+1))

def reserve(gs,N,groups,quota):
    return all(sum(max(0,k-sum(a==j and g>z for a,g in zip(groups,gs)))
                   for j,k in enumerate(quota)) <= z for z in range(N+1))

def supports(S,D):
    return [sum((1<<j) for j,(a,b) in enumerate(S) if a>=p and b>=t) for p,t in D]

def match_subsets(S,D):
    # Independent brute unit assignments: propagate both request and supply masks.
    states={(0,0)}
    for i,e in enumerate(supports(S,D)):
        new=set(states)
        for demand,used in states:
            open_=e&~used
            while open_:
                y=open_&-open_; open_^=y
                new.add((demand|1<<i,used|y))
        states=new
    return {d for d,u in states}

def matching(S,D):
    # Independent augmenting-path search for the 20-slot witnesses.
    M=[-1]*len(S)
    links=[[j for j,(a,b) in enumerate(S) if a>=p and b>=t] for p,t in D]
    def augment(i,seen):
        for j in links[i]:
            if j not in seen:
                seen.add(j)
                if M[j]<0 or augment(M[j],seen):
                    M[j]=i
                    return True
        return False
    return sum(augment(i,set()) for i in range(len(D)))

def gf2(X,I):
    basis={}
    for i in I:
        p=X[i]
        while p:
            a=p.bit_length()-1
            if a in basis: p^=basis[a]
            else:
                basis[a]=p
                break
    return len(basis)

def find(mask,K): return tuple(i for i in range(K) if mask>>i&1)

def all_stock(r,s):
    return [tuple(zip(r,sp)) for sp in set(permutations(s))]

def best_axis_admits(a,d):
    return len(d)<=len(a) and all(x>=y for x,y in zip(sorted(a,reverse=True),sorted(d,reverse=True)))

def lower_easy(r,s,dm,z):
    # dm=(first req margin, second req margin), one independently linkable pool.
    a,b=dm
    K=len(a)
    return max([0]+[sum(x<=p for x in a)+sum(y<=t for y in b)-K
                    for p in set(a) for t in set(b) if gscore(r,s,[(p,t)])[0]>z])

def systematic():
    rng=Random(667)
    shapes=couplings=matroid=quota_checks=0
    for _ in range(240):
        N=rng.randint(0,5); K=rng.randint(0,6)
        r=tuple(sorted(rng.randrange(3) for _ in range(N)))
        s=tuple(sorted(rng.randrange(4) for _ in range(N)))
        D=[(rng.randrange(4),rng.randrange(5)) for _ in range(K)]
        gs=gscore(r,s,D)
        X=[rng.randrange(8) for _ in D]
        ix=[find(u,K) for u in range(1<<K)]
        indep={u for u,ii in enumerate(ix) if gf2(X,ii)==len(ii)}
        rows=all_stock(r,s)
        couple=[match_subsets(S,D) for S in rows]
        opt=[max(map(int.bit_count,J)) for J in couple]
        optm=[max(map(int.bit_count,J&indep)) for J in couple]
        assert min(opt)==formula(gs,N)
        assert min(optm)==formula(gs,N,lambda A:gf2(X,A))
        for arr in (opt,optm): assert set(arr)==set(range(min(arr),max(arr)+1))
        assert max(opt)==max(len(I) for I in ix if
                      best_axis_admits(r,[D[i][0] for i in I]) and
                      best_axis_admits(s,[D[i][1] for i in I]))
        U=set.intersection(*couple)
        virtual={u for u,I in enumerate(ix)
                 if all(g>=j for j,g in enumerate(sorted(gs[i] for i in I),1))}
        assert U==virtual
        group=[rng.randrange(3) for q in D]
        for floor in product(*(range(group.count(j)+1) for j in range(3))):
            def satisfied(I):
                return all(sum(group[i]==j for i in I)>=k for j,k in enumerate(floor))
            assert reserve(gs,N,group,floor)==all(any(satisfied(ix[u]) for u in J) for J in couple)
            quota_checks+=1
        for z in range(N+1):
            stock=canonical(r,s,z)
            for i,su in enumerate(supports(stock,D)):
                if gs[i]<=z: assert su < (1<<z)
        shapes+=1; couplings+=len(rows); matroid+=len(rows)
    print('Tiny known-demand instances/couplings/matroid comparisons/floor vectors:',
          shapes,couplings,matroid,quota_checks)

def both_hidden():
    rng=Random(712)
    checks=hw=0
    for _ in range(60):
        N=rng.randint(1,4)
        r=[rng.randrange(3) for _ in range(N)]
        s=[rng.randrange(3) for _ in range(N)]
        dm=[([rng.randrange(3) for _ in range(k)],[rng.randrange(3) for _ in range(k)])
            for k in (rng.randint(1,3),rng.randint(1,3))]
        sets=[all_stock(*ab) for ab in dm]
        rows=all_stock(r,s); K=sum(len(x[0]) for x in dm)
        idx=range(len(dm[0][0]))
        lefts=sum(1<<i for i in idx)
        best=[]; freqs=[]
        ls=[[lower_easy(r,s,ab,z) for ab in dm] for z in range(N+1)]
        for c,d in product(*sets):
            D=c+d
            assert all(ls[z][0]<=sum(g>z for g in gscore(r,s,c)) and
                       ls[z][1]<=sum(g>z for g in gscore(r,s,d)) for z in range(N+1))
            for S in rows:
                F=match_subsets(S,D)
                best.append(max(map(int.bit_count,F)))
                freqs.append({((u&lefts).bit_count(),(u&~lefts).bit_count()) for u in F})
                hw+=1
        assert min(best)==min(z+sum(a) for z,a in enumerate(ls))
        for ka in range(len(dm[0][0])+1):
            for kb in range(len(dm[1][0])+1):
                gt=all(any(a>=ka and b>=kb for a,b in F) for F in freqs)
                assert gt==all(max(0,ka-aa)+max(0,kb-bb)<=z for z,(aa,bb) in enumerate(ls))
                checks+=1
    print('Both-unlinked joint completions / checked floor vectors:',hw,checks)

def motif():
    r=[2]*10+[1]*10;s=[1]*9+[2]*2+[3]*9
    D=[(2,1)]*9+[(2,2),(1,2)]+[(1,3)]*9
    labels=[0]*10+[2]+[1]*9
    gs=gscore(r,s,D)
    assert Counter(gs)==Counter({10:9,1:1,11:1,9:9})
    seen=[]
    for a in range(10):
        for b in range(3):
            c=10-a-b
            if not 0<=c<=9: continue
            S=[(p,t) for p,row in ((2,(a,b,c)),(1,(9-a,2-b,9-c)))
               for t,m in enumerate(row,1) for _ in range(m)]
            seen.append(matching(S,D))
    assert set(seen)==set(range(11,21))
    accepted={(ka,kb) for ka in range(11) for kb in range(10) if
                  reserve(gs,20,labels,(ka,kb,0))}
    assert accepted=={(ka,kb) for ka in range(11) for kb in range(10) if ka+kb<=10}
    assert matching(canonical(r,s,10),D)==11
    print('Frozen motif, exact supported total / 2x3 tables / robust protected mix:',
          sorted(set(seen)),len(seen),'A<=10, B<=9, A+B<=10')

def tier_control():
    n_cases=0
    for N in range(1,6):
        R=range(N)
        for pt in product(R,repeat=N):
            for dt in (tuple(R),tuple(reversed(R)),tuple(q//2 for q in R)):
                D=tuple(zip(pt,dt))
                if not (best_axis_admits(R,pt) and best_axis_admits(R,dt)): continue
                floor=max(max(Counter(pt).values()),max(Counter(dt).values()))
                assert formula(gscore(R,R,D),N)>=floor
                n_cases+=1
    for k in range(1,8):
        for b in range(1,8):
            N=k*b
            s=[i for i in range(k) for _ in range(b)]
            D=[(i,k-1-i) for i in s]
            assert formula(gscore(s,s,D),N)==b
            assert matching(list(zip(s,reversed(s))),D)==N
    print('Tier/wave bound margin checks and sharp examples:',n_cases,49)

def menu_and_dimension():
    cnt=0
    for n in range(2,8):
        minimizers=[0]*n
        for p in permutations(range(1,n+1)):
            at=[sum(i+b<=n+z for i,b in enumerate(p,1)) for z in range(n)]
            assert all(v>=z for z,v in enumerate(at))
            winners=[z for z,v in enumerate(at) if v==z]
            assert len(winners)<=1
            if winners: minimizers[winners[0]]+=1
            cnt+=1
        assert all(minimizers)
        for z in range(n):
            p=[i+1 for i in perm(z,n)]
            assert sum(i+b<=n+z for i,b in enumerate(p,1))==z
    dd=[(1,1,0),(1,0,1),(0,1,1)]
    opt=[]
    for r in product((0,1),repeat=3):
        S=[r,tuple(1-x for x in r)]
        edges=[{u for u,s in enumerate(S) if all(a>=b for a,b in zip(s,q))} for q in dd]
        opt.append(max(sum(u>=0 for u in a) for a in product((-1,0,1),repeat=3)
                       if len({u for u in a if u>=0})==sum(u>=0 for u in a) and
                       all(u==-1 or u in e for u,e in zip(a,edges))))
    assert min(opt)==1
    print('Scenario-menu permutation controls:',cnt,'; three-axis worst=1, not 0')

if __name__=='__main__':
    motif(); systematic(); both_hidden(); tier_control(); menu_and_dimension()
