"""Independent binary-signature enumerator; no HNF, no author imports."""
from itertools import combinations,product,permutations
from pathlib import Path
from collections import Counter,defaultdict
from math import lcm
import argparse,json
import numpy as np
from sympy import Matrix

MET=np.eye(18,dtype=np.int64)
for i in range(3):
    MET[5*i:5*i+5,5*i:5*i+5]+=1
for i in (15,16,17):
    MET[i,i]=2

pasch=((1,3,5),(1,4,6),(2,3,6),(2,4,5))
C=[tuple((I,-sgn) for I in pasch)+tuple((I,sgn) for I in rem)
   for rem in (((1,2,3),(3,4,5),(1,5,6)),((1,2,6),(2,3,4),(4,5,6)))
   for sgn in (-1,1)]

def wing(i,I,s):
    p=[0]*18
    p[5*i:5*i+5]=[1-2*(j in I) for j in range(1,6)]
    p[15+(i+1)%3]=s
    p[15+(i+2)%3]=1
    return p

def entries(ch):
    W=np.array([wing(i,I,s) for i,k in enumerate(ch) for I,s in C[k]],dtype=np.int64)
    assert np.all(np.diag(W@MET@W.T)==10)
    assert all(abs(int((W@MET@W.T)[i,j]))==2 for i,j in combinations(range(21),2))
    A=Matrix((W@MET).tolist())
    h=A.T.rref()[1]
    assert len(h)==18
    As=A[list(h),:]
    B=2*As.inv()
    D=lcm(*(int(x.q) for x in B))
    K=np.array((D*B).tolist(),dtype=np.int64)
    assert As*(D*B)==2*D*Matrix.eye(18)
    Q=max(sum(abs(int(z)) for z in row) for row in K)
    # Every vector coordinate abs<=Q before screening. Norm and residual
    # products cannot overflow signed int64. No floating operations occur.
    assert 96*max(Q,1)**2 + 10*D*D < 2**63
    rest=np.array([list(A.row(j)) for j in range(21) if j not in h],dtype=np.int64)
    keep=set()
    for batches in range(1<<5):
        words=np.arange(batches*(1<<12),(batches+1)*(1<<12),dtype=np.int64)
        sig=np.concatenate((np.ones((len(words),1),dtype=np.int64),
                            1-2*((words[:,None] >> np.arange(17)) & 1)),axis=1)
        v=sig@K.T
        v=v[np.all(np.abs(v@rest.T)==2*D,axis=1)]
        v=v[np.sum((v@MET)*v,axis=1)==10*D*D]
        for row in v:
            p=tuple(map(int,row))
            q=tuple(-a for a in p)
            assert p!=q
            u=min(p,q)
            assert u not in keep
            keep.add(u)
    # Check all 36 obvious Latin-square control points appear.
    latin=[(a,b,(a+b)%6) for a,b in product(range(6),repeat=2)]
    assert all(sum(x==y for x,y in zip(a,b))<=1 for a,b in combinations(latin,2))
    for ids in latin:
        # z=2h: in each block h has entries 1/6 and one entry -5/6.
        assert D%3==0
        p=tuple(D//3*(1-6*(j==i)) for i in ids for j in range(5))+(0,0,0)
        assert min(p,tuple(-z for z in p)) in keep
    return D,np.array(sorted(keep),dtype=np.int64)

def int_signature_branches(P,D):
    Q=2*D*D
    # v=(P/(2D),height 1/2 A1 coefficient), so v·v'=(P M P'+2D²)/(4D²).
    # Row-profile equality may merge original classes but is exact for compatibility.
    dots=P@MET@P.T
    M=(dots+Q)%(2*Q)==0
    groups=defaultdict(list)
    for i,r in enumerate(M):
        groups[r.tobytes()].append(i)
    classes=list(groups.values())
    S=[set(j for j,e in enumerate(classes) if i!=j and M[c[0],e[0]])
       for i,c in enumerate(classes)]
    def branch(C,F,X):
        if not F and not X:
            yield C
            return
        p=max(F|X,key=lambda a:len(F&S[a]))
        for v in sorted(F-S[p]):
            yield from branch(C+[v], F&S[v], X&S[v])
            F.remove(v);X.add(v)
    B=list(branch([],set(range(len(S))),set()))
    return dots,[[i for j in b for i in classes[j]] for b in B],len(groups)

def cert(P,D,dots,V):
    V=np.array(V,dtype=np.int64)
    graph=abs(dots[np.ix_(V,V)])==2*D*D
    if len(V)==216:
        for i,j in combinations(range(3),2):
            # Independently computed projection-fibre colouring.
            cols=[]
            for x in P[V]:
                a=tuple(x[5*i:5*i+5])
                b=tuple(x[5*j:5*j+5])
                cols.append((min(a,tuple(-z for z in a)),min(b,tuple(-z for z in b))))
            if len(set(cols))<=36 and all(cols[i]!=cols[j] for i,j in zip(*graph.nonzero())):
                return "36-color"
    # On the 93-vertex branches the complement breaks into tiny
    # components. Exact maximum independent set by include/exclude.
    n=len(V)
    assert n==93
    masks=[sum(1<<int(j) for j in np.flatnonzero(~row) if j!=i)
           for i,row in enumerate(graph)]
    todo=(1<<n)-1
    comps=[]
    while todo:
        queue=todo&-todo
        todo^=queue
        comp=queue
        while queue:
            b=queue&-queue;queue^=b
            q=masks[b.bit_length()-1]&todo
            todo^=q;queue|=q;comp|=q
        comps.append(comp)
    assert sorted(c.bit_count() for c in comps)==[8]*6+[15]*3
    from functools import lru_cache
    @lru_cache(None)
    def alpha(P):
        if not P:return 0
        c=P&-P
        return max(alpha(P^c),1+alpha((P^c)&~masks[c.bit_length()-1]))
    maxima=Counter((c.bit_count(),alpha(c)) for c in comps)
    assert sum(a*k for (_,a),k in maxima.items())<=36
    return ",".join(f"{m}:{a}x{k}" for (m,a),k in sorted(maxima.items()))

def check_orbits():
    choices=list(product(range(4),repeat=3))
    canon={}
    for k,Ck in enumerate(C):
        for S in permutations(range(1,7)):
            c=tuple(sorted((tuple(sorted(S[i-1] for i in I)),s) for I,s in Ck))
            assert c not in canon or canon[c]==k
            canon[c]=k
    # Independently generate every seven-clique of the published 40-vertex
    # normalized wing graph by extending only compatible higher labels.
    vlist=list(product(combinations(range(1,7),3),(-1,1)))
    nbr=[{j for j,(J,t) in enumerate(vlist) if j>i and
          3-len(set(I)&set(J))+(s!=t) in (2,3)}
         for i,(I,s) in enumerate(vlist)]
    def enumerate_wings(q,p):
        if len(q)==7:
            yield tuple(sorted(vlist[i] for i in q))
        else:
            for i in sorted(p):
                yield from enumerate_wings(q+[i],p&nbr[i])
    assert set(enumerate_wings([],set(range(40))))==set(canon)
    assert Counter(canon.values())=={k:240 for k in range(4)}
    sizes=[]
    seen=set()
    for last in range(4):
        images=set()
        for perm in permutations(range(3)):
            for flips in product((-1,1),repeat=3):
                decoded=[]
                for i,k in enumerate((0,0,last)):
                    arm=tuple(wing(i,I,s) for I,s in C[k])
                    j=perm[i]
                    blob=[]
                    for v in arm:
                        # Apply the direct coordinate isometry first.
                        out=[0]*18
                        for k in range(3):
                            out[5*perm[k]:5*perm[k]+5]=v[5*k:5*k+5]
                            out[15+perm[k]]=flips[perm[k]]*v[15+k]
                        # Then normalize each horizontal vector's orientation.
                        f=out[15+(j+2)%3]
                        b=[f*x for x in out[5*j:5*j+5]]
                        b.append(-sum(b))
                        I=tuple(k+1 for k,x in enumerate(b) if x==-1)
                        blob.append((I,f*out[15+(j+1)%3]))
                    decoded.append((j,canon[tuple(sorted(blob))]))
                images.add(tuple(k for j,k in sorted(decoded)))
        assert not seen&images
        seen|=images
        sizes.append(len(images))
    assert seen==set(choices) and sizes==[8,8,24,24]
    print("960 actual normalized seven-cliques; physical signed-block orbit sizes",sizes,flush=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fixtures',type=Path,default=None)
    args=parser.parse_args()
    check_orbits()
    for last,expected in enumerate((864,1053,864,1053)):
        ch=(0,0,last)
        D,P=entries(ch)
        assert len(P)==expected
        if args.fixtures:
            fname=args.fixtures/('00'+str(last)+'.json')
            raw=json.loads(fname.read_text())
            d=raw[1]
            assert all(v[18]*2==d and all(2*D*z%d==0 for z in v[:18]) for v in raw[-1])
            assert {tuple(2*D*z//d for z in v[:18]) for v in raw[-1]}=={tuple(x) for x in P}
        dots,branches,groups=int_signature_branches(P,D)
        results=Counter((len(v),cert(P,D,dots,v)) for v in branches)
        assert Counter(map(len,branches))==({216:4} if last%2==0 else {216:4,93:36})
        print(ch,D,len(P),'profile classes',groups,'certificates',dict(sorted(results.items())),flush=True)

if __name__=='__main__':main()
