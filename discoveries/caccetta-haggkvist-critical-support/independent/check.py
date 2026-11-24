#!/usr/bin/env python3
"""Uses exact affine simplex vertices, not determinants/cofactors or floating LP."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd, lcm
import sys

def adjacent(A, i, j):
    return (A[i] >> j) & 1

def triangle_free(A):
    n = len(A)
    for i,j,k in combinations(range(n),3):
        if ((adjacent(A,i,j) and adjacent(A,j,k) and adjacent(A,k,i))
           or (adjacent(A,i,k) and adjacent(A,k,j) and adjacent(A,j,i))):
            return False
    return True

def balance(A):
    n = len(A)
    return tuple(tuple(adjacent(A,j,i)-int(i != j and not adjacent(A,i,j)
                            and not adjacent(A,j,i)) for j in range(n))
                  for i in range(n))

def solve_unique(matrix, rhs):
    """Exact overdetermined Gaussian elimination on an augmented matrix."""
    rows = [[Q(x) for x in r]+[Q(b)] for r,b in zip(matrix,rhs)]
    k = len(rows[0])-1
    for j in range(k):
        p = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if p is None:
            return False  # rank deficient; not a vertex even if consistent
        rows[j], rows[p] = rows[p], rows[j]
        c = rows[j][j]
        rows[j] = [x/c for x in rows[j]]
        for i in range(j+1, len(rows)):
            c = rows[i][j]
            if c:
                rows[i] = [x-c*y for x,y in zip(rows[i], rows[j])]
    if any(row[-1] for row in rows[k:]):
        return None
    answer = [Q(0)]*k
    for j in range(k-1,-1,-1):
        answer[j] = rows[j][-1]-sum(rows[j][i]*answer[i] for i in range(j+1,k))
    return tuple(answer)

def simplex_vertices(A):
    M = balance(A)
    n = len(A)
    whole = solve_unique(list(M)+[[1]*n], [0]*n+[1])
    if whole is not False:
        return [whole] if whole is not None and min(whole)>=0 else []
    vertices = []
    for k in range(1,n):
        for S in combinations(range(n),k):
            rows = [[r[i] for i in S] for r in M] + [[1]*k]
            p = solve_unique(rows, [0]*n+[1])
            if p is False or p is None or min(p)<=0:
                continue
            v = [Q(0)]*n
            for i,q in zip(S,p):
                v[i] = q
            assert all(not sum(x*y for x,y in zip(r,v)) for r in M)
            vertices.append(tuple(v))
    return vertices

@lru_cache(None)
def orderings(n):
    return tuple(permutations(range(n)))

def integral(p):
    p = tuple(int(x*lcm(*(v.denominator for v in p))) for x in p)
    return tuple(x//gcd(*p) for x in p)

def canonical(A, p):
    S = list(combinations(range(len(A)),2))
    return min(("".join(str(1+adjacent(A,q[i],q[j])-adjacent(A,q[j],q[i]))
                         for i,j in S),tuple(p[j] for j in q))
               for q in orderings(len(A)))

def decode(code, n):
    A = [0]*n
    for (i,j), ch in zip(combinations(range(n),2),code):
        if ch=="2":
            A[i] |= 1 << j
        elif ch=="0":
            A[j] |= 1 << i
        else:
            assert ch=="1"
    return tuple(A)

def attachment(A,p):
    """Derive all exact external types by independently adjoining the last vertex."""
    n = len(A)
    feasible=[]
    for states in product((0,1,2), repeat=n):  # 0=external out, 1=nonedge, 2=external in
        if sum(p[i] for i,s in enumerate(states) if s==2)!=sum(p[i] for i,s in enumerate(states) if s==1):
            continue
        B = [A[i] | ((s==2)<<n) for i,s in enumerate(states)]
        B.append(sum((s==0)<<i for i,s in enumerate(states)))
        if triangle_free(B):
            feasible.append(states)
    return feasible

def degree_nesting(code,p):
    A = decode(code, len(p))
    T = attachment(A,p)
    d = tuple(sum(adjacent(A,j,i) for j in range(len(A))) for i in range(len(A)))
    certificates=[]
    for i,j in permutations(range(len(A)),2):
        if d[i]<d[j] and all(s[i]!=0 or s[j]==0 for s in T):
            certificates.append((i,j,d[i],d[j]))
    assert certificates
    print(f"{code}: {len(T)} external types; (subset, superset, indegrees) = {certificates}")
    return d

def main():
    expect = {
      4:{('012010',(1,1,1,1)):6},
      6:{('000010122012022',(4,1,1,1,1,4)):90,
         ('000120101201010',(3,1,1,1,3,3)):360},
    }
    for filename in sys.argv[1:]:
        n = int(filename.rsplit(".",1)[0][-1])
        found = Counter()
        for line in open(filename):
            A=tuple(map(int,line.split()))
            assert len(A)==n
            assert triangle_free(A)
            assert all(any(adjacent(A,j,i) for j in range(n)) for i in range(n))
            vertices = simplex_vertices(A)
            if not vertices or not all(any(v[i]>0 for v in vertices) for i in range(n)):
                continue
            assert len(vertices)==1  # verify uniqueness among positive kernels
            found[canonical(A, integral(vertices[0]))] += 1
        print(f"order {n} positive labelled kernels by canonical type: {dict(found)}")
        assert found==expect.get(n,{})
        if n==6:
            for code,p in sorted(found):
                degree_nesting(code,p)
    four=attachment(decode('012010',4),(1,1,1,1))
    profile=Counter(sum(x==2 for x in a) for a in four)
    assert profile=={0:1,1:4,2:6}
    print(f"C4 attachment types by I size: {dict(sorted(profile.items()))}")

if __name__=='__main__':
    main()
