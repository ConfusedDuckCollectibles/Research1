"""Verify each potentially sharp 7D twisted prism is linearly Hanner.

Producer: verify_twisted_prisms.py. Independently generate Cartesian and
free sums of segments; test exact signed-evaluation matrix equivalence
by colored bipartite graph isomorphism, with explicit permutation check.
"""
from ast import literal_eval
from collections import Counter, defaultdict
from functools import cache
from pathlib import Path
from verify_twisted_prisms import twisted_slacks
from verify_catalog import independent_indices


def leaf_or_dim(t):
    return 1 if not t else sum(map(leaf_or_dim,t[1:]))


@cache
def terms(d):
    if d == 1:
        return ((),)
    out = set()
    for k in range(1,d):
        for a in terms(k):
            for b in terms(d-k):
                for s in (0,1):  # product, free sum
                    aa = a[1:] if a and a[0] == s else (a,)
                    bb = b[1:] if b and b[0] == s else (b,)
                    out.add((s,)+tuple(sorted(aa+bb)))
    return tuple(sorted(out))


def combine(a,b,s):
    v1, q1 = a
    v2, q2 = b
    n1 = len(v1[0])
    n2 = len(v2[0])
    def cart(A,B):
        return [x+y for x in A for y in B]
    def free(A,B):
        return [x+(0,)*n2 for x in A]+[(0,)*n1+y for y in B]
    return (cart(v1,v2),free(q1,q2)) if s == 0 else (free(v1,v2),cart(q1,q2))


@cache
def realize(t):
    if not t:
        return [(-1,),(1,)], [(-1,),(1,)]
    r = realize(t[1])
    for x in t[2:]:
        r = combine(r,realize(x),t[0])
    return r


def slack(t):
    v, q = realize(t)
    S = [tuple((1-sum(xi*yi for xi,yi in zip(x,y)))//2 for x in v) for y in q]
    assert all(sum(xi*yi for xi,yi in zip(x,y)) in (-1,1) for x in v for y in q)
    return S


def graph(S):
    m,n=len(S),len(S[0])
    adj = [tuple(m+j for j,x in enumerate(row) if x == 0) for row in S]
    adj.extend(tuple(i for i,row in enumerate(S) if row[j] == 0) for j in range(n))
    return adj, [0]*m+[1]*n


def iso(S,T):
    if len(S)!=len(T) or len(S[0])!=len(T[0]):
        return None
    A, c0 = graph(S)
    B, d0 = graph(T)

    def rec(c,d):
        while 1:
            f = {}
            def recolor(g, color):
                out = []
                for i, row in enumerate(g):
                    sig = color[i], tuple(sorted(color[j] for j in row))
                    out.append(f.setdefault(sig,len(f)))
                return out
            c1 = recolor(A,c)
            d1 = recolor(B,d)
            hist = Counter(c1)
            if hist != Counter(d1):
                return None
            stable = len(hist)==len(set(c))
            c,d=c1,d1
            if stable:
                break
        if all(z==1 for z in hist.values()):
            j = {color:i for i,color in enumerate(d)}
            p = tuple(j[x] for x in c)
            if all({p[z] for z in A[i]}==set(B[k]) for i,k in enumerate(p)):
                return p
            return None
        c_id = min((key for key,sz in hist.items() if sz>1), key=hist.__getitem__)
        i = c.index(c_id)
        c[i] = len(hist)
        for j,n in enumerate(d):
            if n == c_id:
                e = d.copy()
                e[j] = len(hist)
                z = rec(c,e)
                if z is not None:
                    return z
        return None

    out = rec(c0,d0)
    if out is None:
        return None
    m, n = len(S),len(S[0])
    r, c = out[:m], [x-m for x in out[m:]]
    assert set(r)==set(range(m)) and set(c)==set(range(n))
    assert all(S[i][j]==T[r[i]][c[j]] for i in range(m) for j in range(n))
    return r,c


if __name__ == "__main__":
    rows = [[r for r in x.splitlines() if r and r != '-']
            for x in Path("catalog6-private.txt").read_text().strip().split("-\n")]
    rows = [r for r in rows if r]
    possible = literal_eval(Path("prism7-replay.txt").read_text().splitlines()[1])
    H = terms(7)
    assert len(H) == 180  # simplify associative/commutative terms, not linear isomorphism
    targets = defaultdict(list)
    for i,t in enumerate(H,1):
        assert leaf_or_dim(t) == 7
        S=slack(t)
        g=targets[len(S),len(S[0])]
        if not any(iso(S,T) is not None for j,T in g):
            g.append((i,S))
    assert sum(map(len,targets.values())) == 40
    match = Counter()
    for i, nf, p in possible:
        assert (nf,p)==(645120,'82575360')
        S = twisted_slacks(rows[i-1],7)
        assert len(independent_indices([[1-2*x for x in r] for r in S]))==7
        candidates = targets[len(S),len(S[0])]
        identifications = [j for j,T in candidates if iso(S,T) is not None]
        assert len(identifications)==1
        match.update(identifications)
        print(i,identifications[0])
    assert len(possible)==77 and len(match)==40
    print("77 presentations matched to exactly 40 generated Hanner types")
    print("generation multiplicities histogram",sorted(Counter(match.values()).items()))
