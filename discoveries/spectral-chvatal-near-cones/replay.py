#!/usr/bin/env python3
"""Exact controls, no floating point and no third-party packages."""
from fractions import Fraction as Q
from collections import Counter
from itertools import product, permutations

def sub(a,b):
    return not (a & ~b)

def parity(a):
    return (-1)**a.bit_count()

def matrices(D,i):
    b=1<<i
    U=[a for a in D if not a&b]
    E=[a for a in U if (a|b) in D]
    R=[a for a in U if a not in E]
    cross=[[int(not(a&t)) for t in E] for a in U]
    Vinv=[[parity(a^t) if sub(t,a) else 0 for t in U] for a in U]
    Winv=[[parity(a^t) if sub(t,a) else 0 for t in E] for a in E]
    return U,E,R,cross,Vinv,Winv

def multiply(A,B):
    V=list(zip(*B))
    return [[sum(x*y for x,y in zip(a,b)) for b in V] for a in A]

def trans(A):
    return list(map(list,zip(*A)))

def inertia(A):
    A=[list(map(Q,a)) for a in A]
    p=q=z=0
    while A:
        n=len(A)
        i=next((i for i in range(n) if A[i][i]),None)
        if i is not None:
            x=A[i][i]
            p+=x>0; q+=x<0
            idx=[j for j in range(n) if j!=i]
            A=[[A[j][k]-A[j][i]*A[i][k]/x for k in idx] for j in idx]
            continue
        pair=next(((i,j) for i in range(n) for j in range(i) if A[i][j]),None)
        if pair is None:
            z=n
            break
        i,j=pair
        x=A[i][j]
        p+=1;q+=1
        idx=[k for k in range(n) if k not in (i,j)]
        A=[[A[k][l]-(A[k][i]*A[j][l]+A[k][j]*A[i][l])/x for l in idx] for k in idx]
    return p,q,z

def verify(D,i):
    assert 1<<i in D
    U,E,R,B,V,W=matrices(D,i)
    assert all((a in E) == ((a|(1<<i)) in D) for a in U)
    bt=multiply(multiply(V,B),trans(W))
    expect=[[parity(t) if a==t else 0 for t in E] for a in U]
    assert bt==expect
    Qs=set(R)
    pair=next(((a,a|(1<<j)) for a in R for j in range(max(D).bit_length())
               if not a>>j&1 and (a|(1<<j)) in Qs),None)
    if pair is not None:
        A,C=pair
        j=A^C
        w=[V[U.index(A)][k]+V[U.index(C)][k] for k in range(len(U))]
        assert any(w)
        assert all(bool(x)==(sub(a,C) and bool(a&j)) for x,a in zip(w,U))
        assert all(sum(x*B[k][m] for k,x in enumerate(w))==0 for m in range(len(E)))
        assert all(not(w[u]*w[v]) for u,a in enumerate(U) for v,b in enumerate(U)
                   if not(a&b))
        return False
    assert all(not sub(a,t) for a in R for t in R if a!=t)
    # The complete original matrix, scaled by 2.
    C=[[-sum(parity(F)*((a==F and not(c&F))+(c==F and not(a&F))) for F in R)
        for c in U] for a in U]
    CI=multiply(multiply(V,C),trans(V))
    idx=[U.index(r) for r in R]
    assert [[CI[i][j] for j in idx] for i in idx]==[
        [-2*int(i==j) for j in idx] for i in idx]
    M=[C[k]+[2*x for x in B[k]] for k in range(len(U))]
    M+=[[2*x for x in z]+[0]*len(E) for z in trans(B)]
    idx=U+[e|(1<<i) for e in E]
    assert all(M[k][l]==M[l][k] and (not(a&b) or M[k][l]==0)
               for k,a in enumerate(idx) for l,b in enumerate(idx))
    assert inertia(M)==(len(E),len(U),0)
    assert len(E)==max(sum(bool(a>>j&1) for a in D) for j in range(max(D).bit_length()))
    return True

def all_downsets(n):
    order=sorted(range(1<<n),key=lambda k:(k.bit_count(),k))
    D=set()
    def rec(j):
        if j==len(order):yield sorted(D);return
        a=order[j]
        yield from rec(j+1)
        if all((a^(1<<i)) in D for i in range(n) if a>>i&1):
            D.add(a);yield from rec(j+1);D.remove(a)
    return rec(0)

def generated(F):
    return sorted(a for a in range(1<<max(F).bit_length()) if any(sub(a,b) for b in F))

def shifted(D,p):
    S=set(D)
    return all((a^(1<<i)^(1<<j)) in S for a in S for k,j in enumerate(p)
               for i in p[:k] if a>>j&1 and not a>>i&1)

def main():
    count=Counter()
    for D in all_downsets(4):
        count['downsets']+=1
        for i in range(4):
            if 1<<i in D:
                count['antichain certificate' if verify(D,i) else 'cover-pair stress']+=1
    assert count['downsets']==168
    # Nonpure shifted example: all triples on the first four, plus {0,4}.
    D=generated([7,11,13,14,17])
    assert shifted(D,[0,1,2,3,4]) and verify(D,0)
    # Shifted toward 0 but not fully shifted under any coordinate order:
    # star at 0 and a disjoint-edge matching on the other four vertices.
    D=generated([3,5,9,17,18,12])
    assert all(not shifted(D,p) for p in permutations(range(5)))
    assert verify(D,0)
    # Conventional cross block fails on an exact non-shifted largest star.
    D=generated([3,12])
    assert not verify(D,0)
    print(dict(count))
    print('PASS: all exact certificates and cover-pair stresses verified.')

if __name__=='__main__':
    main()
