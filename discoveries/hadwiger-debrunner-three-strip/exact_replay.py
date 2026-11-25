#!/usr/bin/env python3
"""Exact independent line-arrangement control. No imports from the proof generator."""
from fractions import Fraction as F
from itertools import combinations, product


def points(rects, total):
    """All pairwise non-parallel intersections of the input boundary lines."""
    endpoints = [set(v for b in rects for v in b[k]) for k in range(3)]
    for i, j in combinations(range(3), 2):
        k = 3 - i - j
        for a, b in product(endpoints[i], endpoints[j]):
            p = [None] * 3
            p[i], p[j], p[k] = a, b, total - a - b
            yield tuple(p)


def contains(box, p):
    return all(a <= x <= b for (a, b), x in zip(box, p))


def arrangement(boxes, total):
    masks = {}
    for p in points(boxes, total):
        hit = sum(1 << i for i, box in enumerate(boxes) if contains(box, p))
        if hit:
            masks.setdefault(hit, p)
    return masks


def subsets(n, k):
    for indices in combinations(range(n), k):
        yield sum(1 << i for i in indices)


def meets(s, masks):
    return any(h & s == s for h in masks)


def check43(n, masks):
    bad = []
    for indices in combinations(range(n), 4):
        if not any(meets(sum(1 << i for i in t), masks)
                   for t in combinations(indices, 3)):
            bad.append(indices)
    return bad


def tau(n, masks):
    all_hits = (1 << n) - 1
    level = {0: ()}
    for k in range(n + 1):
        if all_hits in level:
            return k, level[all_hits]
        level = {hit | cover: (*terms, cover) for hit, terms in level.items()
                 for cover in masks}
    raise ValueError('Nonempty members were not all hit')


def cover(masks, n, boxes, witness):
    return all(any(contains(b, p) for p in witness) for b in boxes)


def main():
    triangles, hexagons = [], []
    for i, (a, b) in enumerate([(0, 0), (10, 0), (0, 10), (1, 0), (9, 0),
                                (0, 1), (0, 9), (1, 9), (9, 1)]):
        e = F(1 + i % 3, 40)
        triangles.append(((a,a+10),(b,b+10),(-a-b-10,-a-b)))
        hexagons.append(((a-e,a+10+e),(b-e,b+10+e),(-a-b-10-e,-a-b+e)))
    tri, hexa = [arrangement(L, 0) for L in (triangles, hexagons)]
    nerve = [meets(s, tri) for s in range(1, 512)]
    assert nerve == [meets(s, hexa) for s in range(1, 512)]
    assert not check43(9, tri) and not check43(9, hexa)
    assert all(meets(s, hexa) for s in subsets(9, 2))
    assert tau(9, tri)[0] == tau(9, hexa)[0] == 3
    for b in hexagons:
        p = {p for p in points([b], 0) if contains(b, p)}
        assert len(p) == 6
        assert all(sum(x[i] == b[i][j] for x in p) == 2
                   for i in range(3) for j in range(2))
    print('sharp control: 511 subfamilies, identical nerves; every hexagon has 6 active sides; tau=3')

    def R(*d):
        return tuple((-20,u) for u in d)

    # Sum y=10; all R cuts have min lower <=0 and max upper >=0.
    A, B1, B2, B3, V, S = (
        ((-1,9),(0,14),(-3,-1)),
        R(0,15,15), R(15,0,15), R(15,15,0), R(9,4,10), R(6,8,10))
    hub = [A,B1,B2,B3,V,S]
    segment = ((0,9),(2,11),(-1,-1))
    cut = [segment,*hub[1:]]
    original, trimmed = (arrangement(L,10) for L in (hub,cut))
    assert [meets(s, original) for s in range(1,64)] == [
            meets(s, trimmed) for s in range(1,64)]
    assert not check43(6, original)
    assert not meets(0b101, original) and all(
        meets(s, original) for s in subsets(6,2) if s != 0b101)
    assert tau(6, original)[0] == 2
    endpoints = [(0,11,-1),(7,4,-1)]
    q = (0,0,10)
    assert all(contains(A,p) for p in endpoints)
    assert not any(contains(S,p) for p in endpoints)
    assert contains(S,q) and contains(B2,q)
    assert cover(original, 6, hub, [*endpoints, q])
    for label, i in [('U=B1',1),('V',4),('interior S',5)]:
        one = [p[0] for p in points([segment,hub[i]],10)
               if contains(segment,p) and contains(hub[i],p)]
        assert one
        print(f'negative-upper hub: {label} trace [{min(one)},{max(one)}]')
    assert not any(contains(B2,p) for p in points([segment,B2],10) if contains(segment,p))
    print('negative-upper hub: same 63-intersection nerve after trimming; leaf B2; global tau=2, but chosen endpoints miss S; explicit q=(0,0,10)')

    broken = [A,B1,B2,B3,V,R(6,8,9)]
    assert (0,1,2,5) in check43(6, arrangement(broken,10))
    print('negative control: lower interior d3 from 10 to 9 -> quartet A,B1,B2,S has no triple')

    isolated = (((-1,0),(21,22),(-12,-10)), R(12,12,12),
                R(13,12,14), R(14,15,13))
    masks = arrangement(isolated,10)
    assert not check43(4, masks) and tau(4, masks)[0] == 2
    print('isolated hub control: all three other traces empty, tau=2')


if __name__ == '__main__':
    main()
