"""Conditional exact enumeration of all sharp central seven-polytopes.

Given complete published 6D two-level facet catalogue, every sharp cs 7D
has a candidate realization conv(+- (1,F)), where F is a listed facet.
All of its facets must be 2-level. There are 2^7 candidate normals,
determined by their values +-1 on a linear basis of vertices (1,F).
The script performs only necessary tests for presumed two-level
geometry; any residual candidate is cross-matched to a known sharp
Hanner type in a separate test. It does not assert that the
tentative facet list is complete for arbitrary input F.
"""
import itertools
import hashlib
from collections import Counter
from pathlib import Path
from verify_catalog import (BLOBS, det, independent_indices, mask,
                            symmetric_metrics, triangulation_and_flags)


def twisted_slacks(rows, d):
    raw_columns = [tuple(map(int, r)) for r in rows]
    count = len(raw_columns[0])
    homogeneous = [(1,)*count] + raw_columns
    piv = independent_indices(homogeneous)
    assert len(piv) == d and piv[0] == 0
    pts = list(zip(*(homogeneous[i] for i in piv)))
    assert len(set(pts)) == count
    basis = [pts[i] for i in independent_indices(pts)]
    assert len(basis) == d
    D = det(basis)
    assert D
    br = []
    for p in pts:
        b = tuple(det(basis[:j]+[p]+basis[j+1:]) for j in range(d))
        assert tuple(sum(c*r[k] for c,r in zip(b,basis)) for k in range(d)) == tuple(D*x for x in p)
        br.append(b)
    normals = []
    for suffix in itertools.product((-1,1), repeat=d-1):
        vals = (1,)+suffix
        pattern = [sum(x*y for x,y in zip(b,vals)) for b in br]
        if all(abs(x) == abs(D) for x in pattern):
            plus_slack = tuple(int(x != D) for x in pattern)
            minus_slack = tuple(int(x == D) for x in pattern)
            normals.extend([plus_slack+minus_slack, minus_slack+plus_slack])
    assert len(set(normals)) == len(normals)
    if len(independent_indices([tuple(1-2*x for x in r) for r in normals])) != d:
        return ()
    return normals


if __name__ == "__main__":
    cat = Path("catalog6-private.txt").read_bytes()
    assert hashlib.sha1(b"blob "+str(len(cat)).encode()+b"\0"+cat).hexdigest() == BLOBS[5]
    rows = [[r for r in x.splitlines() if r and r != '-']
            for x in cat.decode("ascii").strip().split("-\n")]
    rows = [r for r in rows if r]
    assert len(rows) == 1150
    counts = Counter()
    equality = []
    for i, r in enumerate(rows, 1):
        maybe = twisted_slacks(r, 7)
        if not maybe:
            counts["rank"] += 1
            continue
        # With the complete true facet list of a full-dimensional
        # polytope, the normals active at each true vertex must span.
        # Our normals are represented injectively by their evaluations
        # on a spanning set. A rank failure is a necessary rejection;
        # an unexpected verification failure below remains a hard error.
        signed = [tuple(1-2*x for x in facet) for facet in maybe]
        if any(len(independent_indices([n for n in signed if n[j] == 1])) != 7
               for j in range(len(signed[0]))):
            counts["active-normal rank"] += 1
            continue
        # A complete primal and polar facet incidence yields only
        # full-dimensional pulling simplices on every listed facet.
        primal_facets = [mask(x == 0 for x in row) for row in maybe]
        dual_facets = [mask(x == 0 for x in row) for row in zip(*maybe)]
        primal_pulls, _ = triangulation_and_flags(primal_facets, len(maybe[0]))
        dual_pulls, _ = triangulation_and_flags(dual_facets, len(maybe))
        if any(len(simplex) != 7 for simplex in primal_pulls + dual_pulls):
            counts["pulling dimensions"] += 1
            continue
        metrics = symmetric_metrics(maybe, 7)
        assert metrics is not None
        flags, vp, bound = metrics
        if vp != bound:
            if vp < bound:
                counts["strict"] += 1
            else:
                raise AssertionError(f"Upper bound exceeded in presentation {i}")
            continue
        equality.append((i, flags, str(vp)))
        counts["equality"] += 1
    print(counts)
    print(equality)
