# An exact inertia certificate for downsets compressed toward one coordinate

**Theorem.** Let D be a finite nontrivial downset (it contains at least
one nonempty set). Suppose there is a coordinate i such that every
nonmaximal member of D omitting i can be extended by adding i.
Then Ellis–Filmus–Friedgut's **Conjecture I** holds for D: there is a
real symmetric matrix supported on disjoint pairs whose number of
nonnegative eigenvalues is exactly the size of the largest star.
The certificate has no zero eigenvalues. In particular, the theorem
holds for every shifted downset, whether pure or not. The hypothesis is
exactly the classical class compressed with respect to one coordinate.

**Fixed-block characterization.** For any specified coordinate i, if
the interblock matrix is fixed to ordinary disjointness as below, an
exact inertia certificate with the star's number of nonnegative
eigenvalues exists **if and only if** the residual family is an antichain.
The obstruction in the other direction is a rank-one positive
semidefinite stress. This converse does not concern arbitrary cross
blocks.

This is a scoped analytic result; it does not settle Conjecture I for
general downsets, nor Conjecture H.

## Primary context

The problem is Conjecture I in
[Ellis, Filmus and Friedgut, *Chvátal's conjecture: a proof from The
Book*, arXiv:2609.28404v3, §5](https://arxiv.org/html/2609.28404v3).
It concerns the Kneser graph on D with the disjointness loop at the
empty set, and the real symmetric eigenvalue bound counting positive
and zero eigenvalues. The original Chvátal conjecture (not this
spectral strengthening) was proved in
[Chang–Liu–Liu, arXiv:2609.19123](https://arxiv.org/abs/2609.19123);
the Book Proof gives another argument and proves the projection
packing analog for all D. The projection-packing theorem alone
does not give the required inertia matrix. For this same class,
[Hunter Snevily, JCTA 61 (1992), 137–141](https://doi.org/10.1016/0097-3165(92)90059-4)
already proved the ordinary maximum-intersecting-star property. Our
claim is the explicit spectral certificate and the sharp fixed-cross
criterion, not a new ordinary Chvátal class. We found no matching
standalone argument for the fixed-block characterization in these
primary sources. AI-assisted; not peer reviewed; no priority claim.

## Proof

Fix a coordinate i, and let D0 consist of the members of D that omit i.
Let E be the link of i: its members are the sets T omitting i for which
T union {i} belongs to D. Thus E is a downset contained in D0. The star
at i consists of the |E| faces T union {i}. Write R = D0 minus E,
s = |E|, and m = |D0|. Within D0 the family R is upward closed.
The hypothesis that nonmaximal faces omitting i can be extended by
i is equivalent to R being an antichain: every F in R has no strict
superset in D0, and no strict superset in D containing i.

Order D by D0 first and the star at i second (indexed by E).
Write a supported symmetric matrix as

    M = [ C    B ]
        [ B^T  0 ].

The second diagonal block must be zero, since any two star members
meet. The fixed cross block in the theorem is

    B(A,T) = 1 if A and T are disjoint, else 0
    (A in D0, T in E).

For a downset U let Z_U(A,S) be 1 if S is a subset of A, else 0.
The inverse has entries (-1)^(|A|-|S|) when S is a subset of A, else
0: every interval in U is the full Boolean interval. In particular
Z_U is invertible. Inclusion–exclusion gives

    B = Z_D0[:,E] diag((-1)^|T| : T in E) Z_E^T,

because every subset of A intersect T belongs to E. Thus simultaneous
congruence by Z_D0^(-1) and Z_E^(-1) sends the cross block to a
nonsingular signed diagonal on E and zero on R.

**If R is an antichain.** For each F in R define d_F on D0 by
d_F(A) = 1 if A is disjoint from F and 0 otherwise. Set

    C = -1/2 sum_(F in R) (-1)^|F| (e_F d_F^T + d_F e_F^T),

where e_F is the F-th coordinate vector. Each summand is supported
on disjoint pairs. In the transformed C block, restrict to rows
and columns indexed by R. Antichainness immediately gives
(Z_D0^(-1) e_F)|_R = e_F. On G in R, the other transform equals

    (Z_D0^(-1) d_F)(G)
        = sum_(S subset G minus F) (-1)^(|G|-|S|)
        = (-1)^|G| * 1[G subset F].

On the antichain this is just (-1)^|F| e_F. The transformed
principal block on R is therefore exactly -I.

For completeness, split the transformed coordinates in the order
E (deletion), E (star), R. The first two blocks form a matrix
[K J; J 0] with diagonal invertible J. This matrix has inertia
(s,s,0), as is seen by eliminating K by a triangular change of
basis; its inverse has zero E-deletion × E-deletion block.
Only E-deletion couples to R. Taking its Schur complement therefore
leaves the R block unchanged. By the law of inertia, M has exactly
s positive and m negative eigenvalues and none zero. The star at i
has size s and is independent; the inertia upper bound forces
every independent family, in particular every other star, to have
size at most s. This proves the theorem with the largest-star value.

The antichain condition is precisely classical compression with
respect to i. Compression means that whenever X belongs to D,
j belongs to X and i does not, replacing j by i still gives a
member of D. Equivalently, any one-element extension F union {j}
of a face F omitting i can be replaced by F union {i}. By downclosure,
every strict extension contains a one-element extension. These are
exactly the equivalent hypothesis and antichain condition above.
If D is shifted with respect to a full order, take i first; no
purity (equal size of maximal faces) is needed.

**If R is not an antichain.** Since R is upward closed within D0 and
D0 is a downset, some A and B = A union {j} both belong to R.
Let z = e_A + e_B in the transformed D0 coordinates and put
w = Z_D0^(-T) z in the original coordinates. The support of w is
exactly the rows S with j in S and S a subset of B. Indeed, the
two inverse coefficients cancel for S a subset of A, and only the
B coefficient appears for the remaining subsets of B. The vector
is nonzero, and B_cross^T w = 0 because z is supported on R.
All its supported sets contain j, so w^T C w = 0 for every admissible
D0 matrix C, even with arbitrary signed entries on disjoint pairs.
This excludes strict negative definiteness on the residual kernel.
The same Schur argument as above shows that such strict definiteness
is necessary to have just s nonnegative eigenvalues. Equivalently
ww^T is a nonzero rank-one dual stress supported on the kernel and
zero on every disjoint pair.

## Scope of the obstruction and variable cross blocks

The fixed-cross obstruction does not prove that Conjecture I fails.
For an arbitrary D, take a largest star at i, and allow unrelated
cross-block variables at the disjoint positions. Then **generic**
real choices avoid every rank-one stress. To see this, any independent
I within D0 has at least |I| neighbors in the star: otherwise combining
I with all star vertices not neighboring it would beat the largest
star, in contradiction with the proved Chvátal theorem. This holds
for every subfamily of I, so Hall's theorem gives a matching from I
into E. A generic cross block has rank |I| on the rows of I because
the corresponding minor has a nonzero matching monomial. There are
only finitely many I, so these generic conditions hold simultaneously.
A rank-one stress ww^T zero on all disjoint pairs must have
independent support, and hence cannot also lie in this generic
kernel. This does not rule out higher-rank stresses or prove existence
of a strict definite completion in general.

## Replay

Run Python 3.10+ with no additional dependencies:

    python3 replay.py

It checks the explicit matrix support, exact Möbius congruence and
small exact inertia for every downset on four named coordinates
and every used choice of coordinate with an antichain residual.
For comparable residuals it checks one exact rank-one obstruction
per coordinate. The code also checks a genuinely nonpure shifted
example; the proof above, not the finite controls, gives all ranks.
