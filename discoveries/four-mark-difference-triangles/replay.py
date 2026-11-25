"""Exact independent verifier for four-mark DTSs and one Wild lift."""
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path


def read_rulers(path):
    return [tuple(map(int, line.split())) for line in Path(path).read_text().splitlines()]


def check(rulers, require_double=False):
    n = len(rulers)
    assert all(len(r) == 4 and r[0] == 0 < r[1] < r[2] < r[3] for r in rulers)
    differences = sorted(b-a for r in rulers for a, b in combinations(r, 2))
    assert differences == list(range(1, 6*n+1))
    total = sum(r[-1] for r in rulers)
    if require_double:
        assert total == 5*n*n+n
        assert sorted(g for r in rulers for g in (r[1], r[3]-r[2])) == list(range(1,2*n+1))
    return total


def f4mul(a,b):
    x = (b if a & 1 else 0) ^ ((b << 1) if a & 2 else 0)
    return x ^ 7 if x & 4 else x  # x^2+x+1 = 0


def wild(a,b):
    v = 12*len(b)+1
    maps = [tuple(f4mul(s,x)^c for x in range(4)) for s in (1,2,3) for c in range(4)]
    assert all({(p[i],p[j]) for p in maps} == set(permutations(range(4),2))
               for i,j in combinations(range(4),2))
    out = list(b) + [tuple(v*y for y in x) for x in a]
    for x in a:
        for y in b:
            for p in maps:
                r = tuple(v*x[i] + y[p[i]] - y[p[0]] for i in range(4))
                out.append(r)
    return out


def main():
    base = Path(__file__).resolve().parent
    for n in (16,17,18,19,20,21,22):
        r = read_rulers(base/f"rulers-{n}.txt")
        assert len(r)==n
        print(f"n={n} scope={6*n} sum={check(r,require_double=True)}")
    base22 = r
    t = wild(base22,base22)
    total = check(t)
    n = len(t)
    assert n==22*265+22==5852
    assert total==2442*(265*265+1)==171491892
    assert total==(5+Fraction(1,133))*n*n+(1-Fraction(22,133))*n
    print(f"n={n} scope={6*n} sum={total}: Wild lift checked")


if __name__ == "__main__":
    main()
