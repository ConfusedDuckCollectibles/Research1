"""Independent small support proof. Starts at a shortest chordless cycle.

Uses fraction-free determinant and principal cofactors, independent of the
all-labelled modular full-rank screen in kernels.cpp.
"""
from itertools import combinations, permutations, product
from functools import lru_cache
from math import gcd
from collections import Counter

def det(M):
    R=[list(r) for r in M]
    size=len(R)
    if size==0:return 1
    old=1;sgn=1
    for i in range(size-1):
        k=next((j for j in range(i,size) if R[j][i]!=0),None)
        if k is None:return 0
        if k!=i:R[i],R[k]=R[k],R[i];sgn=-sgn
        a=R[i][i]
        for j in range(i+1,size):
            for k in range(i+1,size):
                t=R[j][k]*a-R[j][i]*R[i][k]
                assert t%old==0
                R[j][k]=t//old
            R[j][i]=0
        old=a
    return sgn*R[-1][-1]
def minor(M,i,j):
    return [[x for k,x in enumerate(row) if k!=j] for l,row in enumerate(M) if l!=i]
def positive_kernel(M):
    if det(M):return None
    n=len(M)
    for i in range(n):
        v=[(-1)**(i+j)*det(minor(M,i,j)) for j in range(n)]
        if any(v):
            # Rank n-1; a positive kernel exists iff the sole ray has one sign.
            if max(v)<0:v=[-x for x in v]
            return tuple(x//gcd(*v) for x in v) if min(v)>0 else None
    # Rare higher nullity: find supports of all positive extremal rays.
    # For k columns, a strictly positive extreme ray has rank k-1.
    U=set()
    for k in range(1,n):
        for J in combinations(range(n),k):
            for I in combinations(range(n),k-1):
                x=[(-1)**l*det([[M[i][j] for j in J if j!=j0] for i in I])
                   for l,j0 in enumerate(J)]
                if max(x)<0:x=[-a for a in x]
                if min(x)<=0:continue
                if all(sum(M[i][j]*a for j,a in zip(J,x))==0 for i in range(n)):
                    U.update(J)
                    if len(U)==n:raise AssertionError('unexpected positive higher-nullity case')
                    break
    return None
@lru_cache(maxsize=None)
def perms(n):return tuple(permutations(range(n)))
def canonical(n,A,p):
    E=tuple(combinations(range(n),2))
    return min((''.join(str(1+int((q[i],q[j]) in A)-int((q[j],q[i]) in A)) for i,j in E),
                tuple(p[j] for j in q)) for q in perms(n))
def positive_source_free(n):
    found=set();cases=Counter()
    for h in range(4,n+1):
        C={(i,(i+1)%h) for i in range(h)}
        E=[(i,j) for i,j in combinations(range(n),2) if j>=h]
        trips=list(combinations(range(n),3))
        for values in product((-1,0,1),repeat=len(E)):
            A=C|{(i,j) if t==1 else (j,i) for (i,j),t in zip(E,values) if t}
            if len({j for i,j in A})!=n:continue
            if any(((i,j) in A and (j,k) in A and (k,i) in A) or
                   ((j,i) in A and (k,j) in A and (i,k) in A) for i,j,k in trips):continue
            cases[h]+=1
            M=[[int((j,i) in A)-int(i!=j and (j,i) not in A and (i,j) not in A)
                for j in range(n)] for i in range(n)]
            p=positive_kernel(M)
            if p:
                assert all(sum(a*b for a,b in zip(row,p))==0 for row in M)
                found.add(canonical(n,A,p))
    return found,cases

SIX=(('000120101201010',(3,1,1,1,3,3),1,2),
     ('000010122012022',(4,1,1,1,1,4),0,1))
expect={4:{('012010',(1,1,1,1))},5:set(),6:{(s,p) for s,p,i,j in SIX}}
for n in (4,5,6):
    found,cases=positive_source_free(n)
    assert found==expect[n],(n,found)
    print('source-free',n,'cases by planted cycle',dict(cases),'supports',sorted(found))
for s,p,i,j in SIX:
    n=len(p)
    A={(a,b) if t=='2' else (b,a) for (a,b),t in zip(combinations(range(n),2),s) if t!='1'}
    assert sum(x==i for _,x in A)>sum(x==j for _,x in A)
    eligible=0
    for t in product(range(3),repeat=n):
        # 0 means outside sends to support, 1 means missing, 2 means receives.
        if sum(p[a] for a in range(n) if t[a]==1)!=sum(p[a] for a in range(n) if t[a]==2):continue
        if any(t[a]==0 and t[b]==2 for a,b in A):continue
        eligible+=1
        assert t[j]!=0 or t[i]==0
    print(s,'exact type nesting',eligible,'types; N^-_ext',j,'subset of',i)
