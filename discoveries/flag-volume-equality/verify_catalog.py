"""Exact independent Winter-bound replay on published two-level slack tables.
Usage: python verify_catalog.py DIR   (files named catalog[D]-private.txt)
Input snapshot: https://github.com/mmacchia/2L_enum, files 1d.txt to 6d.txt.
The data are not created or redistributed by this source.
"""
from fractions import Fraction
from functools import cache
from itertools import product
from math import factorial, lcm
from pathlib import Path
import hashlib
import sys

BLOBS = [
    'd9550d53bf22febd6223f0b97e7f9d22c239e623',
    'deaf417147e1f6043b03e296a84c1221d6f70f62',
    '804bf5ed2e60d865090c24823bc144ba4948c396',
    '6121c2144670fe2435197922cd47f9293ddb1178',
    '02e010ce6ea8bba5859a03e57a6093a7c36c6d95',
    '47256c326a347f16de65f44fc944d0e67bbae723',
]
TOTALS = [1, 2, 5, 19, 106, 1150]
SYMMETRIC = [1, 1, 2, 4, 13, 45]
HANNER_TYPES = [1, 1, 2, 4, 8, 18]


def det(M):
    """Fraction-free Bareiss integer determinant with pivots."""
    M = [list(row) for row in M]
    n = len(M)
    if n == 0:
        return 1
    assert all(len(row) == n for row in M)
    sgn, prev = 1, 1
    for j in range(n - 1):
        pivot = next((i for i in range(j, n) if M[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            M[pivot], M[j] = M[j], M[pivot]
            sgn = -sgn
        p = M[j][j]
        for i in range(j + 1, n):
            for k in range(j + 1, n):
                v = M[i][k] * p - M[i][j] * M[j][k]
                M[i][k], rem = divmod(v, prev)
                assert rem == 0
            M[i][j] = 0
        prev = p
    return sgn * M[-1][-1]


def independent_indices(rows):
    """Rank profile; input rows can be rectangular."""
    basis = {}
    selected = []
    for i, row in enumerate(rows):
        z = list(map(Fraction, row))
        for j, p in sorted(basis.items()):
            if z[j]:
                z = [x-z[j]*y for x, y in zip(z, p)]
        j = next((j for j, x in enumerate(z) if x), None)
        if j is not None:
            basis[j] = [x/z[j] for x in z]
            selected.append(i)
    return selected


def mask(values):
    return sum(1 << i for i, b in enumerate(values) if b)


def indices(x):
    while x:
        y = x & -x
        yield y.bit_length() - 1
        x ^= y


def triangulation_and_flags(facets, n):
    """Pull with fixed vertex order, count flags separately on face lattice."""
    full = (1 << n) - 1
    assert len(facets) == len(set(facets)) and all(0 < F < full for F in facets)

    @cache
    def subfacets(F):
        cuts = sorted({F & G for G in facets if 0 < (F & G) < F},
                      key=int.bit_count, reverse=True)
        keep = []
        for C in cuts:
            if not any(C & K == C for K in keep):
                keep.append(C)
        return keep

    @cache
    def count(F):
        if F.bit_count() == 1:
            return 1
        return sum(count(G) for G in subfacets(F))

    @cache
    def pull(F):
        u_bit = F & -F
        u = u_bit.bit_length() - 1
        if F == u_bit:
            return [(u,)]
        return [s+(u,) for G in subfacets(F) if not (u_bit & G)
                for s in pull(G)]

    simplices = [s for F in facets for s in pull(F)]
    assert len(simplices) == len(set(simplices))
    return simplices, count(full)


def integer_normalized_volume(points, facets, d):
    """Return d! times volume of full-dimensional P, 0 in its interior."""
    simplices, flags = triangulation_and_flags(facets, len(points))
    assert all(len(S) == d for S in simplices)
    terms = [abs(det([points[i] for i in S])) for S in simplices]
    assert terms and min(terms) > 0
    return sum(terms), flags


def symmetric_metrics(rows, d):
    # Replace slack 0/1 entries by signed facet evaluation +1/-1.
    matrix = list(zip(*(tuple(1-2*int(x) for x in r) for r in rows)))
    assert len(matrix) == len(set(matrix))
    pairs = set(matrix)
    if not all(tuple(-x for x in r) in pairs for r in matrix):
        return None
    cols = list(zip(*matrix))
    cj = independent_indices(cols)
    assert len(cj) == d
    points = [tuple(r[j] for j in cj) for r in matrix]
    ri = independent_indices(points)
    assert len(ri) == d
    basis = [points[i] for i in ri]
    denom = det(basis)
    assert denom
    qs = []
    for col in cols:
        rhs = [col[i] for i in ri]
        numer = [det([row[:k]+(rhs[i],)+row[k+1:] for i, row in enumerate(basis)])
                 for k in range(d)]
        q = tuple(Fraction(x, denom) for x in numer)
        assert tuple(sum(x*y for x, y in zip(q, p)) for p in points) == col
        qs.append(q)
    scale = lcm(*(x.denominator for q in qs for x in q))
    dual = [tuple(int(scale*x) for x in q) for q in qs]
    # Cross-validate exact affine dimension of each indexed facet.
    for facet_col in cols:
        in_face = [p for p, x in zip(points, facet_col) if x == 1]
        assert len(independent_indices([p+(1,) for p in in_face])) == d
    for facet_row in matrix:
        in_face = [p for p, x in zip(dual, facet_row) if x == 1]
        assert len(independent_indices([p+(1,) for p in in_face])) == d

    facets = [mask(x == 1 for x in c) for c in cols]
    dual_facets = [mask(x == 1 for x in r) for r in matrix]
    primal, fa = integer_normalized_volume(points, facets, d)
    polar, fa2 = integer_normalized_volume(dual, dual_facets, d)
    assert fa == fa2
    bound = (2**d) * fa
    prod_num = Fraction(primal*polar, scale**d)
    return fa, prod_num, bound


def replay_hansen():
    """Exact volumes not needed for the forbidden-face witnesses."""
    for name, edges, edge, f in [
        ("P4", [(0,1),(1,2),(2,3)], (1,2), (1,1,1,1)),
        ("C4", [(0,1),(1,2),(2,3),(0,3)], (0,1), (-1,0,1,-1)),
    ]:
        vs = [v for v in product((0,1), repeat=4)
              if all(v[i]+v[j] <= 1 for i,j in edges)]
        F = [v for v in vs if v[edge[0]]+v[edge[1]] == 1]
        G = [v for v in vs if v[edge[0]]+v[edge[1]] == 0]
        assert len(independent_indices([v+(1,) for v in vs])) == 5
        assert len(independent_indices([v+(1,) for v in F])) == 4
        vals = [[sum(a*b for a,b in zip(f,v)) for v in S] for S in (F,G)]
        widths = [max(a)-min(a) for a in vals]
        assert widths == [1,2]
        print(f"{name}: facet rank 4; widths {widths[0]},{widths[1]}")


def replay_catalog(path):
    for d, blob in enumerate(BLOBS, 1):
        raw = (path / f"catalog{d}-private.txt").read_bytes()
        assert hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest() == blob
        items = [x.splitlines() for x in raw.decode().strip().split("-\n")]
        items = [[r for r in rows if r and r != '-'] for rows in items]
        items = [rows for rows in items if rows]
        assert len(items) == TOTALS[d-1]
        cs, upper_equal, gaps = 0, [], []
        for i, rows in enumerate(items, 1):
            assert all(len(s) == len(rows[0]) and set(s) <= {"0","1"} for s in rows)
            out = symmetric_metrics(rows, d)
            if out is None:
                continue
            cs += 1
            flags, prod_num, bound = out
            assert bound >= prod_num
            if bound == prod_num:
                assert flags == 2**d*factorial(d)
                assert prod_num == 4**d*factorial(d)
                upper_equal.append(i)
            else:
                gaps.append(Fraction(bound)-prod_num)
        assert cs == SYMMETRIC[d-1] and len(upper_equal) == HANNER_TYPES[d-1]
        print(f"d={d}: all={len(items)}, cs={cs}, equality={len(upper_equal)}, "
              f"strict={len(gaps)}, least positive scaled gap="
              f"{min(gaps) if gaps else 'n/a'}; equality indices {upper_equal}",
              flush=True)


if __name__ == "__main__":
    replay_hansen()
    if len(sys.argv) > 1:
        replay_catalog(Path(sys.argv[1]))
