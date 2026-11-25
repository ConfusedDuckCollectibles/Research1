#!/usr/bin/env python3
"""An independent integer incidence replay of (A), not a volume computation.
Reads locally obtained Bohn et al. slack matrices; imports no candidate code.
Every accepted inclusion uses direct vertex-difference matches; every rejection
is a strict support certificate in the ambient signed-slack coordinates.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import sys


def parts(path):
    current = []
    for line in path.read_text().splitlines() + ['-']:
        if line == '-':
            if current:
                yield current
            current = []
        elif line:
            current.append(line)


def mask(seq):
    return sum(int(s) << i for i, s in enumerate(seq))


def indices(bits):
    return [i for i in range(bits.bit_length()) if bits >> i & 1]


def check_poly(facet_strings):
    # Native catalog entries are 0/1 facet slack rows. This injectively
    # embeds the polytopes and records antipodal vertices by complements.
    m, n = len(facet_strings), len(facet_strings[0])
    full_v, full_f = (1 << n)-1, (1 << m)-1
    vv = [mask(c) for c in zip(*facet_strings)]
    if {v ^ full_f for v in vv} != set(vv):
        return None
    bitfacets = {full_v ^ mask(r) for r in facet_strings}
    assert len(vv) == len(set(vv)) and len(bitfacets) == m
    strict = 0
    good = 0

    @lru_cache(None)
    def quad_diffs(face):
        return frozenset((vv[a] & ~vv[b], vv[b] & ~vv[a])
                         for a,b in combinations(indices(face),2))

    @lru_cache(None)
    def go(face):
        nonlocal strict, good
        if face.bit_count() <= 1:
            return None
        cuts = {f & face for f in bitfacets} - {0,face}
        subfacets = []
        for f in sorted(cuts,key=int.bit_count,reverse=True):
            if not any(f & h == f for h in subfacets):
                subfacets.append(f)
        for a in subfacets:
            b = face ^ a
            da, db = quad_diffs(a), quad_diffs(b)
            da = da | frozenset((j,i) for i,j in da)
            for i,j in db:
                if (i,j) in da: # exact difference vector, sufficient
                    continue
                # Test separability by the sign vector of the target.
                # Its target value is popcnt(i|j); equality for a hull
                # generator requires all signed positions to match.
                if not any((i & x == i and j & y == j) for x,y in da):
                    strict += 1
                    return (face, a, i, j)
                raise AssertionError('test requires an exact LP fallback')
            good += 1
            failure = go(a)
            if failure:
                return failure
        return None

    failure = go(full_v)
    return failure, (good,strict)


def main(path):
    expected = (1,1,2,4,13,45)
    expected_good = (1,1,2,4,8,18)
    for d in range(1,7):
        good, bad = [], []
        count_inclusions = 0
        entries = list(parts(path / ('catalog%d-private.txt'%d)))
        for i, table in enumerate(entries,1):
            out = check_poly(table)
            if out is None: continue
            failure, (checked, strict) = out
            count_inclusions += checked
            if failure:
                assert strict == 1
                bad.append((i,failure))
            else:
                good.append(i)
        assert len(good)+len(bad)==expected[d-1]
        assert len(good)==expected_good[d-1]
        print(f'd={d} cs={len(good)+len(bad)}: exact difference inclusion {len(good)}, '
              f'exact signed-support failure {len(bad)}; {count_inclusions} successful incidences')
        print('   equality indices',good)
        for i, ws in bad:
            print(f'   reject row {i}: parent={ws[0]:x}, facet={ws[1]:x}, '
                  f'plus={ws[2]:x}, minus={ws[3]:x}')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
