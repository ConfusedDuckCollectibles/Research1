#!/usr/bin/env python3
"""Independent physical CSS controls and abstract parity-quotient exhaustion."""
from itertools import combinations, product

def dot(a, b):
    return (a & b).bit_count() & 1

def expand(rows):
    words = {0}
    for row in rows:
        if row not in words:
            words |= {row ^ x for x in words}
    return words

def basis(rows):
    piv = {}
    for row in rows:
        while row:
            i = row.bit_length() - 1
            if i not in piv:
                piv[i] = row
                break
            row ^= piv[i]
    return tuple(piv.values())

def rank(rows):
    return len(basis(rows))

def coset_rep(x, vectors):
    for v in sorted(vectors, reverse=True):
        x = min(x, x ^ v)
    return x

def least_distance(n, stab, anti, allowed=None):
    """Enumerate physical supports, not block-compressed or logical weights."""
    br = basis(stab)
    if allowed is None:
        allowed = tuple(range(n))
    for d in range(1, len(allowed) + 1):
        for support in combinations(allowed, d):
            v = sum(1 << i for i in support)
            if all(not dot(v, row) for row in anti) and coset_rep(v, br):
                return d
    return float('inf')

def physical(block_sizes, b_edges, a_support=()):
    blocks = []
    xrows, zrows = [], []
    n = 0
    for w in block_sizes:
        blocks.append(tuple(range(n, n+w)))
        xrows.extend((1<<i) | (1<<(i+1)) for i in range(n, n+w-1))
        n += w
    for strut in a_support:
        xrows.append(sum(1 << blocks[b][0] for b in strut))
    zrows = [sum(1 << i for b in support for i in blocks[b]) for support in b_edges]
    return n, blocks, xrows, zrows

def regional(d, D, triples, singles):
    B, parts = [], []
    block = 0
    for a in [(d//2, d-d//2, d-d//2)]*triples + [(d,)]*singles:
        heads = []
        for size in a:
            heads.append(block)
            B.extend((i, i+1) for i in range(block, block+size-1))
            block += size
        if len(a) > 1:
            B.append(tuple(heads))
        parts.append(a)
    return physical([D]*block, B)

def extract(n, D, xrows, zrows):
    """Infer intrinsic components directly with physical X membership tests."""
    xb, zb = basis(xrows), basis(zrows)
    classes, seen = [], set()
    for i in range(n):
        if i in seen:
            continue
        block = [j for j in range(i, n) if not coset_rep((1 << i) ^ (1 << j), xb)]
        seen.update(block)
        unit = sum(1 << j for j in block)
        if not coset_rep(1 << i, xb) or not coset_rep(unit, zb):
            continue
        classes.append(block)
    L = [block for block in classes if len(block) >= D]
    zL = [sum(((z >> block[0]) & 1) << j for j, block in enumerate(L)) for z in zrows]
    K = expand(zL)
    U = expand(z for z in K if z.bit_count() <= 2)
    scores = []
    # Construct components as equivalence classes of U, unrelated to prior union-find.
    remaining = {i for i in range(len(L)) if (1 << i) not in U}
    while remaining:
        i = min(remaining)
        comp = {j for j in remaining if ((1 << i) ^ (1 << j)) in U}
        remaining -= comp
        scores.append(sum(1 << j for j in comp))
    R = {sum(dot(z, c) << j for j, c in enumerate(scores)) for z in K}
    assert len(K) == len(U)*len(R)
    assert all(not z or z.bit_count() >= 3 for z in R)
    return L, scores, K, R

def show(name, model, D, expected, restricted_expected=None):
    n, _blocks, X, Z = model
    assert all(not dot(x, z) for x in X for z in Z)
    k = n - rank(X) - rank(Z)
    d = least_distance(n, X, Z)
    dz = least_distance(n, Z, X)
    L, C, K, R = extract(n, D, X, Z)
    allowed = tuple(i for block in L for i in block)
    mask = sum(1 << i for i in allowed)
    h = len(allowed) - rank(z & mask for z in Z) - (rank(X) - rank(x & ~mask for x in X))
    rd = least_distance(n, X, Z, allowed)
    parities = sorted(x for x in R if x)
    if len(parities) > 1:  # Test the intrinsic disjoint-triple reconstruction.
        parities = [x for x in parities if x.bit_count() == 3]
    if rank(parities) == len(R).bit_length()-1 and all(
            not (a & b) for a, b in combinations(parities, 2)):
        supports = [[i for i in range(len(C)) if (x >> i) & 1] for x in parities]
        covered = {i for supp in supports for i in supp}
        patterns = [(i,) for i in range(len(C)) if i not in covered]
        patterns.extend(pair for supp in supports for pair in combinations(supp, 2))
        xb = basis(X)
        witness_sizes = []
        for indexes in patterns:
            word = sum(C[i] for i in indexes)
            lift = sum(1 << block[0] for i, block in enumerate(L) if (word >> i) & 1)
            assert all(not dot(z, lift) for z in Z)
            if coset_rep(lift, xb):
                witness_sizes.append(lift.bit_count())
        assert rd == min(witness_sizes, default=float('inf')), name
    actual = (n, k, d, dz, max((z.bit_count() for z in Z), default=0),
              len(R).bit_length()-1, h, len(allowed))
    assert actual == expected, (name, actual, expected)
    if restricted_expected is not None:
        assert rd == restricted_expected, (name, rd)
    print(f"{name}: n,k,dx,dz,Wz,rankR,h,nL={actual}; restricted dx={rd}")
    return R

def physical_controls():
    show("one odd-distance triple", regional(3,2,1,0), 2, (10,2,3,2,6,1,2,10), 3)
    show("two triples plus chain", regional(3,2,2,1), 2, (26,5,3,2,6,2,5,26), 3)
    show("even dx and rank-one", regional(4,3,1,0), 3, (18,2,4,3,9,1,2,18), 4)
    show("degree-two projection r=0", regional(3,2,0,2), 2, (12,2,3,2,4,0,2,12), 3)
    show("d=1: triple plus singleton", physical([2]*4,[(0,1,2)]), 2, (8,3,1,2,6,1,3,8), 1)
    show("d=1 confined to S", physical([2]*3+[1]*3,[(0,1,2)],[(3,4,5)]),
         2, (9,4,1,2,6,1,2,6), 2)
    show("surplus overlapping images", physical([2]*6,[(0,1,2),(3,4,5)],[(0,1,3,4)]),
         2, (12,3,2,2,6,2,3,12), 2)
    show("one new class on a triple", physical([2]*4,[(1,2),(0,1,3)],[(0,1,2)]),
         2, (8,1,2,2,6,1,1,8), 2)
    show("zero logical image on triple", physical([2]*6,[(1,2),(3,4),(0,1,3)],[(0,1,2),(0,3,4)]),
         2, (12,1,1,2,6,1,1,12), 1)
    show("h=0 and an isolated degenerate triple",
         physical([2]*5+[1]*3,[(1,2),(3,4),(0,1,3)],[(0,1,2),(0,3,4),(5,6,7)]),
         2, (13,2,1,2,6,1,0,10), float('inf'))
    R = show("overlap violates formula without disjointness",
             physical([3]*6,[(0,1,2),(2,3,4),(0,4,5)]), 3, (18,3,3,3,9,3,3,18), 3)
    assert sorted(x.bit_count() for x in R) == [0,3,3,3,3,4,4,4]
    assert 18 < 3 * (3*3 - min(3,3//2)*(3//2))
    show("rank-one four-set at Wz=4D", physical([2]*4,[(0,1,2,3)]), 2,
         (8,3,2,2,8,1,3,8), 2)

def all_subspaces(dim):
    spaces, pending = {frozenset([0])}, [frozenset([0])]
    full = set(range(1 << dim))
    while pending:
        W = pending.pop()
        for a in full - W:
            V = W | frozenset(a ^ w for w in W)
            if V not in spaces:
                spaces.add(V)
                pending.append(V)
    return spaces

def parity_profile(t, j, d):
    if d == 1 or not j:
        return t
    a = t - 2*j
    return j*d + min(a, a*(d//2))

def plans():
    for m in range(1,6):
        yield m, ()
        for t in range(3,m+1):
            yield m, (t,)
    yield 6, (3,3)
    yield 7, (3,3)
    yield 7, (3,4)
    yield 8, (4,4)

def abstract_controls():
    spaces = {}
    tally = 0
    for m, sizes in plans():
        idx = 0
        regions = []
        for t in sizes:
            regions.append(tuple(range(idx,idx+t)))
            idx += t
        regions += [(i,) for i in range(idx,m)]
        masks = [sum(1 << i for i in region) for region in regions if len(region)>1]
        K = [word for word in range(1 << m) if all(not dot(word, mask) for mask in masks)]
        dim = len(K).bit_length()-1
        if dim not in spaces:
            spaces[dim] = all_subspaces(dim)
        # Independent coordinate labels, via enumeration of span in the given order.
        generators = basis(K)
        transfer = [0]*(1 << dim)
        for i, b in enumerate(generators):
            for x in range(1 << i):
                transfer[x+(1<<i)] = transfer[x] ^ b
        word_region = [[x for x in K if not x & ~sum(1<<i for i in region)] for region in regions]
        survivors = tuple(transfer[code] for code in range(1 << dim))
        piece_cases = [x for P in word_region for x in P if x.bit_count() <= 2 and x]
        quotient_cases = []
        for w in spaces[dim]:
            W = {transfer[x] for x in w}
            # Use the full span to compute dimension increments, not local labels.
            h = dim - (len(w).bit_length()-1)
            accumulated, changes = W, []
            for words in word_region:
                next_set = {x ^ y for x in accumulated for y in words}
                changes.append(len(next_set).bit_length() - len(accumulated).bit_length())
                accumulated = next_set
            assert sum(changes) == h and len(accumulated) == len(K)
            quotient_cases.append((h, changes, [x for x in survivors if x not in W],
                                   [x for x in piece_cases if x not in W]))
        choices = (1,2,3) if m <= 6 else (1,2)
        counts = 0
        for weights in product(choices, repeat=m):
            cost = [0]*(1<<m)
            for x in range(1,1<<m):
                low = x & -x
                cost[x] = cost[x^low] + weights[low.bit_length()-1]
            weights_of_regions = [sum(weights[i] for i in region) for region in regions]
            for h, changes, nontrivial, pieces in quotient_cases:
                d = min((cost[x] for x in nontrivial), default=float('inf'))
                restricted = min((cost[x] for x in pieces), default=float('inf'))
                assert d == restricted
                # Lower targets model what happens if d_X is set in a different region.
                for target in range(1, int(d)+1 if h else 6):
                    regional_lb = 0
                    for region, j, available in zip(regions, changes, weights_of_regions):
                        if len(region) == 1:
                            lb = j*target
                        else:
                            lb = parity_profile(len(region),j,target)
                        assert available >= lb
                        regional_lb += lb
                    assert sum(weights) >= regional_lb
                    r = len(sizes)
                    if all(t == 3 for t in sizes):
                        proposed = h*target - min(r,h//2)*(target//2)
                        if target >= 2:
                            stronger = proposed + (h%2 if h < 2*r else 0) + 3*max(0,r-(h+1)//2)
                        else:
                            stronger = max(3*r,h+r)
                        assert sum(weights) >= proposed
                        assert sum(weights) >= stronger
                    counts += 1
        print(f"abstract: m={m}, regions={sizes}, kernel dim={dim}, quotient-test cases={counts}")
        tally += counts
    print(f"ABSTRACT TOTAL {tally} weighted target/quotient cases")

if __name__ == '__main__':
    physical_controls()
    abstract_controls()
