#!/usr/bin/env python3
"""Independently audit all four-mark difference and span claims, Python 3 only."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

FOLDER = Path(__file__).resolve().parent

def seed(n):
    return [tuple(map(int, s.split()))
            for s in (FOLDER / "fixtures" / f"seed-{n}.txt").read_text().splitlines()]

def certify(R, perfect=True, minimum_sum=False):
    n = len(R)
    assert n and all(len(r) == 4 and r[0] == 0 and all(x < y for x, y in zip(r, r[1:])) for r in R)
    D = Counter(b - a for r in R for a, b in combinations(r, 2))
    outer = Counter(z for r in R for z in (r[1], r[3] - r[2]))
    scope = max(r[-1] for r in R)
    spans = sum(r[-1] for r in R)
    assert D.total() == len(D) == 6 * n
    assert outer.total() == len(outer) == 2 * n
    assert sum(D) + sum(outer) == 4 * spans
    assert scope >= 6 * n and spans >= 5*n*n + n
    if perfect:
        assert D == Counter(range(1, 6*n + 1))
    if minimum_sum:
        assert outer == Counter(range(1, 2*n + 1))
        assert spans == 5*n*n + n
    return n, scope, spans

def f4mul(a, b):
    """F2[X]/(X^2+X+1); bit represents coefficient."""
    m = 0
    for i in range(2):
        if (b >> i) & 1:
            m ^= a << i
    return m ^ 7 if m & 4 else m

def lift(A, B):
    perms = [tuple(f4mul(a, i) ^ b for i in range(4))
             for a, b in product(range(1, 4), range(4))]
    assert len(set(perms)) == 12
    direct_pairs = set((a,b) for a,b in product(range(4), repeat=2) if a != b)
    for i, j in direct_pairs:
        assert {(p[i], p[j]) for p in perms} == direct_pairs
    v = 12 * len(B) + 1
    out = list(B) + [tuple(v*a for a in x) for x in A]
    for x, y, p in product(A, B, perms):
        r0 = v*x[0] + y[p[0]]
        out.append(tuple(v*x[i] + y[p[i]] - r0 for i in range(4)))
    return out

def coefficient(B):
    return (5 + Fraction(1, 6*B + 1), 1 - Fraction(B, 6*B + 1))

def main():
    seeds = {}
    for n in (6, 10, 15, 16, 17, 18, 19, 20, 21):
        seeds[n] = seed(n)
        ans = certify(seeds[n], minimum_sum=True)
        assert ans[0] == n
        print("simultaneously optimal:", *ans)
    for a, b in ((16,18), (18,16), (18,18), (15,15), (21,21)):
        ans = certify(lift(seeds[a], seeds[b]))
        v = 12*b + 1
        assert ans == (v*a+b, 6*(v*a+b), v*v*(5*a*a+a)+(5*b*b+b))
        print(f"independent exact product {a} × {b}:", *ans)
    for base in (15, 16, 18, 21):
        coef_quad, coef_lin = coefficient(base)
        v = 12*base + 1
        for j in range(1, 12):
            L = (v**j - 1) // 12
            S = (5*base*base + base) * (v**(2*j) - 1) // (v*v - 1)
            assert S == coef_quad*L**2 + coef_lin*L
        print("base", base, "=> quadratic coefficient", coef_quad, "linear coefficient", coef_lin)
    broken = seeds[18][:] + seeds[18][:1]
    try:
        certify(broken)
    except AssertionError:
        print("negative control: duplicate ruler correctly rejected")
    else:
        raise AssertionError("negative control unexpectedly accepted")

if __name__ == "__main__":
    main()
