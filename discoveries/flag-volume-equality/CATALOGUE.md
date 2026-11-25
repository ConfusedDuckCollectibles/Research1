# Reproduction on the published two-level catalogues

Use the primary [Bohn et al. enumeration](https://arxiv.org/abs/1703.01943),
Lemma 5 (2-level combinatorial equivalence implies affine equivalence)
and §6 Table 1. Download the slack tables 1d.txt through 6d.txt from
[Macchia's accompanying repository](https://github.com/mmacchia/2L_enum)
as catalog1-private.txt through catalog6-private.txt. The original
tables are not part of any release. The Python standard-library command
    python verify_catalog.py .
checks each table's Git blob digest and independently computes flag
counts and rational volume products. Every upper equality case has
both $2^dd!$ flags and volume product $4^d/d!$. In dimensions
1 through 6, the numbers of cs cases are $1,1,2,4,13,45$; the
upper equality counts are $1,1,2,4,8,18$, exactly the Hanner type
counts (Winter's primary slides 6–7). Hanner polytopes themselves
are 2-level and sharp, so this gives the entire classification in
that range, independent of the symmetric Mahler lower theorem.
Finite completeness comes from the published enumeration, not
a new proof of its coverage.

The [ProofAtlas flag/KL research map](https://www.proofatlas.ai/collaboration/mahler-volume-product-conjecture/)
is dated August 6, 2026, before Winter's October paper, and targets
a weighted flag *lower* bound on the Mahler volume product and
small-vertex/facet resonance rays. Its stated directional shadow
rigidity concerns lower Mahler minimizers, not equality in
Winter's upper bound; no face-difference criterion or two-level
rigidity for Winter equality was located there.

## Dimension seven, from a six-dimensional catalogue only

Run, in order:
    python verify_twisted_prisms.py > prism7-replay.txt
    python verify_hanner_residual.py > hanner7-replay.txt

Completeness for the problem is independent of a published dimension-seven
flag-volume computation. A sharp centrally symmetric 7-polytope is 2-level
by PROOF.md, and hence is the twisted prism over any one of its facets, namely
$T(F)=\operatorname{conv}\{\pm(1,v):v\in\operatorname{vert}(F)\}$.
The facet is itself 2-level (the general hereditary fact for 2-level
polytopes is Bohn et al. Lemma 6) and so appears, up to affine equivalence,
among the 1,150 six-dimensional input types. An affine isomorphism of
facets lifts to a linear isomorphism of the twisted prisms.

The script chooses independent homogeneous coordinates for $F$ from
the all-ones row and normalized binary slack rows, then a linear basis
of 7 vertices $(1,v)$. Any 2-level facet of $T(F)$, normalized to
have the equations $q\cdot X=1$ on it and $q\cdot X=-1$ on its
antipode, must give a $\{\pm1\}$-vector on that basis. Enumerate
these 128 vectors (64 up to common sign) and check all vertex
evaluations exactly by integer determinants. Conversely every vector
which takes only these two values on all vertices gives a true facet:
the signed basis vertices chosen on its positive level are linearly
independent and lie in that plane. Thus this lists exactly all
two-level facets, and lists all facets whenever the candidate is
2-level. No check assumes that arbitrary twisted prisms are 2-level.

For a 2-level candidate, the signed incidence matrix has full rank,
the normals active at every vertex span, and every pulling simplex
of a primal or polar facet has seven vertices. The rank is checked
directly on evaluation vectors because the vertices span. These
necessary checks are explicit; an unexpected assertion is not treated
as rejection. If the candidate is 2-level, the computed volume products
and flag counts are the true ones. We can therefore safely exclude any
candidate that fails a necessary check or has unequal product.
Of 1,150 bases, 348 fail the pulling-dimension check, 725 give strict inequality
with the tentative list, and
77 give tentative equality. These are counts of candidate base
presentations, not of 7-polytope types; the 725 may include not
truly 2-level candidates. Every residual has $645120=2^7\,7!$
flags and scaled product $82575360=4^7\,7!$ in that replay.

For the positive geometrical check, independently construct products
and free sums of intervals, including the polar vertices by the
opposite recursion. Associative/commutative deduplication gives
180 formal expressions at dimension seven; exact row/column
permutation reduces these to 40 signed evaluation arrays.
Each residual presentation has an exactly checked bijection of its
rows and columns to one such array; all 40 occur. Signed evaluation
arrays suffice even for candidates not yet known to be 2-level:
the linear map from the vertex span to the vector of evaluations
under the listed normals is injective (these normals span), and
identifies the convex hull with the convex hull of its evaluation
vectors. Therefore matching arrays up to permutations give a
linear equivalence of the actual convex hulls, without relying on
a tentative lattice or on the recent symmetric Mahler theorem.

The complete prism verification takes a few minutes on a typical CPU; the Hanner matching takes a few seconds.
The replay proves the centrally symmetric equality classification
through dimension seven, conditional only on the published
six-dimensional 2-level catalogue's completeness. No independent
dimension-eight replay was performed. Winter's November 2025 slides
report a full central 2-level check through eight, and the [ULB
enumeration repository](https://github.com/ulb/tl) points to
the 1,790 central types in that dimension, inside a larger
compressed catalogue. Treat dimension eight as source-reported.

For an independent local check through dimension six, run
`python incidence_replay.py .` on the same files, after the digest
checking in the main volume replay. This version checks the
face-difference inclusions directly with integer vectors, not
volumes; it reproduces the equality indices and prints certificates
for all 32 strict failures.
