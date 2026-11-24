#!/usr/bin/env python3
"""Independent exact CSS checks using low pivots and physical meet-in-the-middle."""
from collections import defaultdict,Counter
from itertools import combinations
import argparse
import json
from pathlib import Path

def dot(a,b):
    return (a&b).bit_count()%2

def closure(words):
    s={0}
    for w in words:
        if w not in s:
            s.update([v^w for v in s])
    return frozenset(s)

def low_basis(rows):
    ans={}
    for v in rows:
        while v:
            j=(v&-v).bit_length()-1
            if j in ans:v^=ans[j]
            else:
                ans[j]=v
                break
    return ans

def rank(rows):
    return len(low_basis(rows))

def orth_basis(rows,n):
    bz=low_basis(rows)
    answer=[]
    for j in range(n):
        if j in bz: continue
        v=1<<j
        for q,y in sorted(bz.items(),reverse=True):
            if dot(y,v):v|=1<<q
        answer.append(v)
    assert len(answer)==n-len(bz) and all(dot(y,v)==0 for y in rows for v in answer)
    return answer

def half_words(cols):
    sig=0;last=0
    for q in range(1<<len(cols)):
        word=q^(q>>1)
        flip=word^last
        if flip:sig^=cols[flip.bit_length()-1]
        last=word
        yield word,sig

def mitm(checks,trivial,n):
    normals=orth_basis(trivial,n)
    allrows=checks+normals
    code=[sum((1<<j) for j,r in enumerate(allrows) if r>>q&1) for q in range(n)]
    m=n//2
    low=(1<<len(checks))-1
    groups=defaultdict(dict)
    for x,key in half_words(code[:m]):
        d=groups[key&low]
        lg=key>>len(checks)
        opt=(x.bit_count(),x,lg)
        if lg not in d or opt<d[lg]:d[lg]=opt
    neighbors={s:sorted(g.values())[:2] for s,g in groups.items()}
    answer=(n+1,0)
    for x,key in half_words(code[m:]):
        lg=key>>len(checks)
        for w,v,k in neighbors.get(key&low,()):
            if k!=lg:
                answer=min(answer,(w+x.bit_count(),v|(x<<m)))
                break
    assert answer[1] and all(not dot(x,answer[1]) for x in checks)
    assert any(dot(x,answer[1]) for x in normals)
    return answer

def classical_lift(X,Z,widths):
    # Physical numbering is block-contiguous, but local parity is checked independently.
    start=[sum(widths[:i]) for i in range(len(widths))]
    pieces=[((1<<w)-1)<<s for s,w in zip(start,widths)]
    xx=[3<<q for s,w in zip(start,widths) for q in range(s,s+w-1)]
    xx += [sum(1<<s for i,s in enumerate(start) if m>>i&1) for m in X]
    zz = [sum(b for i,b in enumerate(pieces) if m>>i&1) for m in Z]
    return xx,zz

def compressed_from_duals(X,Z,n):
    def sig(x,dual):return sum(dot(x,v)<<i for i,v in enumerate(dual))
    dx,dz=orth_basis(X,n),orth_basis(Z,n)
    gs=defaultdict(int)
    for i in range(n):gs[sig(1<<i,dx)]|=1<<i
    masks=[m for f,m in gs.items() if f and sig(m,dz)]
    a=[row_view(x,masks,True) for x in X]
    b=[row_view(z,masks,False) for z in Z]
    return masks,a,b

def check_region(n,X,Z,D,dX):
    masks,a,b=compressed_from_duals(X,Z,n)
    m=len(masks)
    L=sum(1<<i for i,g in enumerate(masks) if g.bit_count()>=D)
    S=((1<<m)-1)^L
    physS=sum(g for i,g in enumerate(masks) if S>>i&1)
    physL=sum(g for i,g in enumerate(masks) if L>>i&1)
    # Solve the actual restricted physical kernel, then take its image
    # modulo the original physical Z stabilizers. Deleted factors do not
    # need to be included in this support calculation.
    cols=[i for i in range(n) if physS>>i&1]
    small=[sum((x>>i&1)<<j for j,i in enumerate(cols)) for x in X]
    cut=orth_basis(small,len(cols))
    embed=[sum((x>>j&1)<<i for j,i in enumerate(cols)) for x in cut]
    ell=rank(Z+embed)-rank(Z)
    pa=rank(x&S for x in a)
    assert ell==S.bit_count()-pa-rank(b)+rank(x&L for x in b)
    local=None
    present=all((x&L).bit_count()<=2 for x in b)
    atom2=[v for v in closure(x&L for x in b) if 0<v.bit_count()<=2]
    intrinsic=rank(atom2)==rank(x&L for x in b)
    assert not present or intrinsic
    n0=n-physS.bit_count()-physL.bit_count()
    k=n-rank(X)-rank(Z)
    if intrinsic:
        components=exact_graph(L,atom2,a)
        row_a=closure(a)
        eligible=[r for r in components if r not in row_a]
        assert rank(list(a)+eligible)-rank(a)==k-ell
        local=min((x.bit_count() for x in eligible),default=None)
        dregs=[x.bit_count() for x in kernels(closure(b),m)-row_a if not x&S]
        assert min(dregs,default=None)==local
        assert physL.bit_count()>=D*dX*len(eligible)>=D*dX*(k-ell)
    return dict(n0=n0,nS=physS.bit_count(),nL=physL.bit_count(),ell=ell,rS=pa,
                presentation=present,intrinsic=intrinsic,local=local,
                bound=n0+physS.bit_count()+D*dX*(k-ell))

def larger(path):
    with path.open() as f: source=json.load(f)
    controls={x['control']:x for x in source['controls']}
    X1,Z1=classical_lift([], [3,6], [6]*3)
    steane=[int('1010101',2),int('1100110',2),int('1111000',2)]
    X2,Z2=classical_lift(steane,steane,[2]*7)
    # On Steane, coordinates representing F2^3 vectors x,y,x^y give a
    # weight-three word, not a check (all nonzero checks weigh four).
    word2=sum(1<<(2*i) for i in (0,1,2))
    # Each of 3 repeated blocks from the all-equal outer X word.
    word1=(1<<0)|(1<<6)|(1<<12)
    xs=X1+[x<<18 for x in X2];zs=Z1+[z<<18 for z in Z2]
    xr=xs[:];xr[0]^=xr[-1]
    xa,za=classical_lift([], [7], [3]*3)
    xone,zone=classical_lift([], [], [6])
    actual={
       'pair_only_asymmetric':(18,X1,Z1),
       'repeated_steane_bias_control':(14,X2,Z2),
       'mixed_check_budget_exact_saturation':(32,xs,zs),
       'mixed_X_links_clear_small_logicals':(32,xs+[word1|(word2<<18)],zs),
       'basis_change_intrinsic_pair_invariance':(32,xr,zs),
       'three_block_parity_boundary':(9,xa,za),
       'one_sector_protected':(6,xone,zone),
       'pinned_and_fixed_factors':(21,X1+[1<<18,3<<19],Z1+[3<<19]),
    }
    for name,(n,X,Z) in actual.items():
        assert not any(dot(x,z) for x in X for z in Z)
        c=controls[name]
        read=lambda name:[sum(1<<q for q in supp) for supp in c['presented_generators'][name]]
        XC,ZC=read('X'),read('Z')
        # All eight physical presentations have the same spaces as the
        # reported inputs, though the coupled example may choose a
        # different term if the producer's minimum-word tie breaks differ.
        assert rank(X+XC)==rank(X)==rank(XC)
        assert rank(Z+ZC)==rank(Z)==rank(ZC)
        dX,x=mitm(Z,X,n)
        dZ,z=mitm(X,Z,n)
        assert (n-rank(X)-rank(Z),dX,dZ)==(c['k'],c['d_X'],c['d_Z'])
        col=lambda R:max(sum(v>>i&1 for v in R) for i in range(n))
        row=lambda R:max((v.bit_count() for v in R),default=0)
        assert [row(XC),row(ZC)]==c['physical_generator_max_row']
        assert [col(XC),col(ZC)]==c['physical_generator_max_col']
        q=check_region(n,XC,ZC,c['target_D'],dX)
        expected=c['certificate']
        assert (q['nS'],q['nL'],q['ell'],q['rS'])==(
            expected['n_small'],expected['n_large'],expected['ell_Z_on_small'],expected['rank_A_on_small'])
        assert q['presentation']==expected['applies_actual_large_check_projection']
        if q['presentation']:
            assert q['local']==expected['exact_dX_support_large']
            # The producer's number is the post-removal lower bound.
            assert q['bound']==expected['n_bound']
        if name=='three_block_parity_boundary':
            D=c['target_D']
            assert c['physical_generator_max_row'][1]==3*D
            assert n<c['k']*dX*D
        print(name, 'n,k,dX,dZ',n,c['k'],dX,dZ,'physical witnesses',hex(x),hex(z),q)
    X,Z=classical_lift([], [3,6,12], [6]*4)
    Z[0]^=Z[-1]
    assert max(z.bit_count() for z in Z)==24
    assert (mitm(Z,X,24)[0],mitm(X,Z,24)[0])==(4,6)
    q=check_region(24,X,Z,6,4)
    assert not q['presentation'] and q['intrinsic'] and q['bound']==24
    print('hidden projected low-weight basis: [[24,1,(4,6)]], W_Z=24,',q)

def subspaces(n):
    all={frozenset({0})}
    layer=set(all)
    for _ in range(n):
        layer={frozenset(tuple(s)+tuple(y^w for y in s))
               for s in layer for w in range(1<<n) if w not in s}
        all|=layer
    return sorted(all,key=lambda s:(len(s),tuple(sorted(s))))

def greedy(S):
    got={0}; rows=[]
    for x in sorted(S,key=lambda x:(x.bit_count(),x)):
        if x not in got:
            rows.append(x)
            got.update([v^x for v in got])
    return rows

def row_view(v,masks,parity):
    return sum((dot(v,b) if parity else int(bool(v&b)))<<i for i,b in enumerate(masks))

def from_signatures(U,V,n):
    # Signatures in F2^n/U. Use least member of each exact coset as key,
    # rather than union-find over pair checks.
    groups=defaultdict(int)
    for q in range(n):
        key=min((1<<q)^v for v in U)
        groups[key]|=1<<q
    masks=[g for k,g in groups.items() if k and g not in V]
    A=frozenset(row_view(x,masks,True) for x in U)
    B=frozenset(row_view(z,masks,False) for z in V)
    return masks,A,B

def dim(S):
    n=len(S)
    assert n and n&(n-1)==0
    return n.bit_length()-1

def kernels(A,m):
    return frozenset(x for x in range(1<<m) if all(not dot(x,a) for a in A))

def exact_graph(large, Sz, A):
    rem=set(i for i in range(large.bit_length()) if large>>i&1)
    connected=defaultdict(set)
    forced=set()
    for z in Sz:
        ls=[i for i in rem if z>>i&1]
        assert len(ls)<=2
        if len(ls)==1:forced.update(ls)
        if len(ls)==2:
            a,b=ls;connected[a].add(b);connected[b].add(a)
    indicators=[]
    while rem:
        w=rem.pop(); part={w};front={w}
        while front:
            shell=set().union(*(connected[v] for v in front))&rem
            rem-=shell;part|=shell;front=shell
        if not part&forced:indicators.append(sum(1<<i for i in part))
    return indicators

def small(nmin,nmax):
    assert nmin>=1 and nmax<=6
    dims=Counter()
    for n in range(nmin,nmax+1):
        allsp=subspaces(n)
        for U in allsp:
            orthU=frozenset(x for x in range(1<<n) if all(not dot(x,u) for u in U))
            for V in allsp:
                if len(V)>=len(orthU):break
                if not V<=orthU:continue
                dims[n,'positive-CSS']+=1
                k=n-dim(U)-dim(V)
                x0=kernels(V,n)-U;z0=orthU-V
                dX=min(x.bit_count() for x in x0)
                dZ=min(z.bit_count() for z in z0)
                masks,A,B=from_signatures(U,V,n); M=len(masks); T=(1<<M)-1
                W=[s.bit_count() for s in masks]
                assert M-dim(A)-dim(B)==k
                assert all(z & g in (0,g) for g in masks for z in V)
                assert all(v.bit_count()>=3 for v in A if v)
                assert all((1<<i) not in B for i in range(M))
                X=kernels(B,M)-A; Z=kernels(A,M)-B
                assert min(v.bit_count() for v in X)==dX
                assert min(sum(w for i,w in enumerate(W) if v>>i&1) for v in Z)==dZ
                PX=greedy(U);PZ=greedy(V)
                Wx=max((x.bit_count() for x in PX),default=0)
                Wz=max((z.bit_count() for z in PZ),default=0)
                for D in range(1,dZ+1):
                    dims[n,'D-targets']+=1
                    L=sum(1<<i for i,w in enumerate(W) if w>=D); S=T^L
                    physicalS=sum(masks[i] for i in range(M) if S>>i&1)
                    physicalL=sum(masks[i] for i in range(M) if L>>i&1)
                    ell=dim(closure(list(V)+[x for x in z0 if x&~physicalS==0]))-dim(V)
                    rS=dim({u&S for u in A})
                    assert ell==S.bit_count()-rS-dim(B)+dim({z&L for z in B})
                    assert 0<=ell<=k and ell<=S.bit_count()-rS<=(Wx-1)*rS
                    assert dim(closure(list(U)+[x for x in x0 if x&~physicalL==0]))-dim(U)==k-ell
                    if physicalS.bit_count()<D: assert ell==0
                    if Wz<2*D and dX>=2: assert ell==k
                    if all(x.bit_count()<=2 for x in PX):
                        assert A=={0} and not S
                    viewed=[row_view(z,masks,False)&L for z in PZ]
                    if any(z.bit_count()>2 for z in viewed):
                        assert Wz>=3*D
                        continue
                    dims[n,'graph-qualifying']+=1
                    gens=exact_graph(L,viewed,A)
                    local0=[v.bit_count() for v in x0 if v&~physicalL==0]
                    candidates=[z for z in gens if z not in A]
                    active=rank(list(A)+candidates)-dim(A)
                    assert active==k-ell and physicalL.bit_count()>=D*dX*active
                    assert min(local0,default=None)==min((v.bit_count() for v in candidates),default=None)
        print('exhaustive n=',n,', positive CSS pairs =',dims[n,'positive-CSS'],
              ', D targets =',dims[n,'D-targets'],', graph-qualifying =',dims[n,'graph-qualifying'],flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=Path(__file__).with_name('checks.json'))
    p.add_argument('--nmax',type=int,default=5)
    p.add_argument('--small-only',action='store_true')
    args=p.parse_args()
    if not args.small_only:larger(args.source)
    if args.nmax:small(1,args.nmax)
    print('PASS',flush=True)

if __name__=='__main__':
    main()
