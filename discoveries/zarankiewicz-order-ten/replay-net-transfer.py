#!/usr/bin/env python3
"""Exact arithmetic and positive controls. No packages, data, or search required."""
from collections import Counter
from itertools import combinations, product


def choose2(n):
    return n * (n - 1) // 2


def dimensions(q, j):
    n = q * (q + 1)
    d = q + j
    return n, d, n * (q + 1) - d


def crossing(q, n, e, s):
    return s * (n + 1 - s) - e + (n - 1) + choose2(q + 1 - s)


def metsch_rhs(j):
    return 8*j**3 - 18*j**2 + 8*j + 4


def check_arithmetic():
    candidates = sizes = 0
    for q in range(4, 129):
        j_cube = max(j for j in range(1, q//2 + 1) if 8*j**3 <= 3*q)
        assert 3*q > metsch_rhs(j_cube)
        for j in range(1, q//2 + 1):
            n, d, e = dimensions(q, j)
            budget = choose2(n) - n*choose2(q + 1) + q*d
            assert budget == (q*q+(2*j-1)*q)//2
            left = crossing(q, n, e, q + 2)
            right = crossing(q, n, e, n)
            assert left == q*q + j - 2
            assert right-left == (n-q-2)*(n-3*q-1)//2
            assert budget < left < right
            assert 2*(n-d) > q*d
            assert (q*q-1)*(q+1)-(e-((q+1)**2-1)) == d-1
            if q <= 40:
                for s in range(q+2, n+1):
                    assert crossing(q, n, e, s) >= left
                    sizes += 1
            candidates += 1
    for q,j,upper in ((6,2,285), (10,2,1197), (30,3,28796)):
        assert 3*q > metsch_rhs(j)
        assert dimensions(q,j)[2]-1 == upper
    assert all(a*a + b*b != 30 for a,b in product(range(6), repeat=2))
    print(f"arithmetic: {candidates} parameter pairs; {sizes} individual high-degree sizes")
    print("examples: q=6 -> 285 (weaker than 283); q=10 -> 1197; q=30 -> 28796")


def truncated_plane(q):
    # Affine F_q^2 and one point at infinity per nonvertical slope.
    # Equivalently delete the vertical point at infinity and the infinity line.
    points = set(product(range(q), repeat=2)) | {(-1,m) for m in range(q)}
    graph = {(m,b): {(-1,m)} | {(x,(m*x+b) % q) for x in range(q)}
             for m,b in product(range(q), repeat=2)}
    graph.update({(-1,c): {(c,y) for y in range(q)} for c in range(q)})
    return points, graph


def extract(q, punched_flags):
    j = len(punched_flags)
    n, d, e = dimensions(q, j)
    points, lines = truncated_plane(q)
    for line, point in punched_flags:
        assert (0, 0) not in lines[line]
        lines[line].remove(point)
    assert len(lines) == len(points) == n and sum(map(len, lines.values())) == e
    assert max(Counter(p for line in lines.values() for p in line).values()) == q+1
    point_pairs = set()
    for line in lines.values():
        for pair in combinations(sorted(line), 2):
            assert pair not in point_pairs
            point_pairs.add(pair)
    assert choose2(n)-len(point_pairs)+sum(choose2(q+1-len(line)) for line in lines.values()) == (q*q+(2*j-1)*q)//2
    p = (0,0)
    groups = sorted((line-{p} for line in lines.values() if p in line), key=len)
    small, *large = groups
    assert len(small) == q-1 and list(map(len, large)) == [q]*q
    assert set().union(*groups) == points-{p}
    assert sum(map(len,groups)) == len(points)-1
    raw = []
    for line in lines.values():
        if p not in line:
            row = []
            for group in large+[small]:
                hit = line & group
                assert len(hit) <= 1
                row.append(next(iter(hit)) if hit else None)
            raw.append(tuple(row))
    assert len(raw) == q*q-1
    assert sum(row.count(None) for row in raw) == d-1
    for i,k in combinations(range(q+1),2):
        visible = [(row[i],row[k]) for row in raw if row[i] is not None and row[k] is not None]
        assert len(visible) == len(set(visible))
    full = [i for i in range(q) if all(row[i] is not None for row in raw)]
    restored = []
    for i in full:
        counts = Counter(row[i] for row in raw)
        assert set(counts) == large[i]
        assert sorted(counts.values()) == [q-1]+[q]*(q-1)
        restored.append(next(symbol for symbol in counts if counts[symbol] == q-1))
    use_small = len(full) < q+1-j
    if use_small:
        assert len(full) == q-j
        counts_small = Counter(row[-1] for row in raw if row[-1] is not None)
        assert set(counts_small) == small and set(counts_small.values()) == {q}
        new_symbol = object()
        complete = [tuple(row[i] for i in full)+(row[-1] if row[-1] is not None else new_symbol,) for row in raw]
        complete.append(tuple(restored)+(new_symbol,))
        assert sum(row[-1] is new_symbol for row in complete) == q
    else:
        complete = [tuple(row[i] for i in full) for row in raw]+[tuple(restored)]
    assert len(complete) == q*q and len(complete[0]) >= q+1-j
    for i,k in combinations(range(len(complete[0])),2):
        cols = {row[i] for row in complete}, {row[k] for row in complete}
        assert len(cols[0]) == len(cols[1]) == q
        assert {(row[i],row[k]) for row in complete} == set(product(*cols))
    return len(full), use_small, len(complete[0])


def check_controls():
    tests = [
        ("q=5, two different large groups", 5, [((-1,1),(1,0)), ((-1,2),(2,2))], (3,True,4)),
        ("q=5, two holes in one large group", 5, [((-1,1),(1,0)), ((-1,2),(2,0))], (4,False,4)),
        ("q=5, one small and one large hole", 5, [((0,1),(0,1)), ((-1,1),(1,0))], (4,False,4)),
        ("q=7, three different large groups", 7, [((-1,1),(1,0)), ((-1,2),(2,2)), ((-1,3),(3,6))], (4,True,5)),
    ]
    for title,q,flags,expected in tests:
        result = extract(q,flags)
        assert result == expected, (title,result,expected)
        before, added, after = result
        print(f"control: {title}: {before} full large columns; small completion={added}; OA({q*q},{after},{q},2)")


if __name__ == "__main__":
    check_arithmetic()
    check_controls()
    print("PASS")
