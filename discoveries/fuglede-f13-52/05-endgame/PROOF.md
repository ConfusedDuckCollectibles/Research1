# Size 52 over F13³: the final eight-pencil branch

Step 5 of the size-52 proof, 25 November 2026 UTC. Combined with the [layer exclusion](../02-layers/PROOF.md) and the [both-K9 exclusion](../04-both-k9/PROOF.md), this argument excludes cardinality 52. Size 78 remains open; nothing here proves the full prime-13 conjecture.

## Inputs and notation

E and A are a hypothetical spectral pair of size 52 in dual copies of F13³. D(E) is the set of projective nonzero differences. U_A comprises normal vectors (up to scaling) with all thirteen parallel planes containing exactly four A-points. Thus D(E) is contained in U_A (cyclotomic orthogonality); similarly D(A) is contained in U_E. The audited preceding reductions give line cap 4, plane cap 10, maximum uniform-pencil size at most 8 (on both), and all mixed moments through degree 7 on both. See [step 4](../04-both-k9/PROOF.md) for the last two inputs. Only the upper bound is used here.

For each projective normal v of E let q_E(v) count unordered pairs of equal level, and define k_E(v)=(q_E(v)-78)/13. Two incidence identities are
```
sum_v k_E(v) = 330,
sum_{v perpendicular to d} k_E(v) = 18 + t_E(d),
```
where t_E(d) counts unordered differences of direction d. To verify, a fixed nonzero difference has 14 orthogonal normals among all 183, and among the 14 normals perpendicular to a fixed d it has all 14 if parallel to d and exactly one otherwise. The unordered pair total is 1326; divide the sums after subtracting 78 per normal. For d in V_A, the complement of U_A, we have t_E(d)=0. Summing on V_A and reversing the sum gives
```
18(183-m) = sum_v k_E(v)(14-n_v),
m = |U_A|, n_v = |U_A ∩ v^perp|.
```
Since n_v≤8 and k_E is nonnegative, this implies m≤73 for either side.

## 1. Every nonuniform plane profile has k at least 2

For each v, its integer profile g(s) counts members in the level v·x=s. As an F13-valued function, write this unique reduced polynomial of degree at most 12. The moments give sum_s s^j g(s)=0 for j=0,...,7. Each of these isolates minus the coefficient at 12−j, since the exponent never reaches 24 and the constant sum is 13=0. Thus g has degree at most 4. The recorded counts are its least nonnegative residues, since plane cap is 10. Every profile therefore appears in the complete enumeration over 13^5 polynomials, retaining only all counts ≤10 and integer total 52. Exact replay finds no k=1, and k=0 only for the constant 4 profile. k is integral, also directly because sum g²=0 in F13 (degree ≤8). Nonnegativity follows from minimizing equal-level pairs at four points per level. Therefore k_E(v)≥2 for every v in V_E, and |V_E|≥110.

## 2. A direction line with at most three D-points is forced

Suppose every projective line meets D(E) in at least four points. Then the m-point set U_A has line occupancies n between 4 and 8. The standard point-pair incidence identities are sum n=14m and sum n(n−1)=m(m−1); there are 183 lines. Thus
```
0 ≥ sum_lines (n−4)(n−8) = m²−155m+5856.
```
This forces integer m≥66 (the quadratic decreases through 65, where it is 6). Let L8 be the number of lines with n=8. Rearranging the earlier weighted identity gives
```
18(73−m) = sum_v k_E(v)(8−n_v) ≥ 2(110−L8).
```
Indeed at least 110 normals belong to V_E; at most L8 of them can have n_v=8, by point–line duality. On the other hand, for each n∈[4,8],
```
choose(n,2) ≥ 6+4(n−4)+6·1_{n=8}.
```
Sum this to obtain L8≤floor[(choose(m,2)−56m+1830)/6]. For m=66,...,73 this upper bound is respectively 46, 48, 50, 52, 54, 56, 59, 61 (replay recomputes the simple values). Every one of the eight possibilities contradicts the deficit inequality. Hence some projective direction line meets D(E) in at most three points.

## 3. Three orientations, covariance, contradiction

Choose that projective direction line as the horizontal direction plane. Parallel affine slices y∈F13 have average occupancy four. Any noncollinear four points determine at least four different projective directions in odd characteristic. Since line cap is four, a slice of more than four points also contains four noncollinear points. So all slices have exactly four points on one affine line; their horizontal orientations use at most three distinct directions, all in the chosen low pencil.

Two parallel four-point lines in different slices are contained in an affine plane. Fix an affine coordinate along their common direction. The cross-directions, excluding the common direction, are in bijection with the difference set of their respective four abscissas. By Cauchy–Davenport there are at least seven; the upper pencil cap eight in U_A includes the common direction, so there are exactly seven. Vosper's theorem for two 4-sets with sumset size 7≤13−2 makes both 4-term arithmetic progressions with the same step up to sign. Equivalently, enumerate all 715 4-subsets of F13 and their pairs; the [exact replay](check-endgame.js) checks the identical local fact. The step of a 4-term arithmetic progression is unique up to sign, since the ordered nonzero difference multiplicities are 3 at ±step, 2 at ±2step, and 1 at ±3step, with all six residues distinct. Thus each horizontal orientation used in at least two slices has a single fixed step vector up to sign, independent of slice offsets.

Let P∈F13² be the horizontal vector of a point, and C(y)=1/4 times the sum of the four P in slice y. The mixed-moment identities for P y^j, j=0,...,6, imply each of the two centroid coordinates as a function of y is a polynomial of degree at most 5, by the same coefficient isolation as above. It follows that sum_y C(y)C(y)^T is the zero matrix in F13, since each entry has degree at most 10. Also sum_{E} P P^T=0 by degree-two moments. Hence the sum over all slices of the unnormalized, centered second-moment matrices is zero:
```
sum_y sum_{P in slice y} (P−C(y))(P−C(y))^T
 = sum_E P P^T − 4 sum_y C(y)C(y)^T = 0.
```
Any slice contributes a scalar times d d^T for some nonzero vector d of its orientation. The rank-one symmetric matrices for at most three distinct directions in F13² are independent: with two directions choose axes; any third has both coordinates nonzero and hence nonzero off-diagonal entry. For a repeated orientation use its common step d. Each associated slice has four offsets (−3/2,−1/2,1/2,3/2)d from its centroid, so its matrix is 5 d d^T. A repeated orientation occupying 2 through 12 slices contributes nonzero coefficient 5 times its multiplicity in F13, contradiction by independence. Some orientation must repeat, so the sole remaining possibility is that one orientation occupies all thirteen slices.

In that final possibility project along its common d. The thirteen distinct projected rich points lie in an affine plane and have no three collinear (three parallel rich lines in one affine plane would give 12 selected points, violating cap 10). Their 78 unordered pairs therefore determine at least 13 distinct quotient directions, since no direction has more than six disjoint pairs on thirteen points. For each quotient direction the associated projective line through d contains at least seven non-d points in D(E), by the same cross-direction count. These sets of non-d points are disjoint for distinct quotient directions. Therefore |D(E)|≥1+13·7=92, contradicting |U_A|≤73.

This proves there cannot be a residual spectral pair of size 52, completing the preceding reductions. The size-78 case remains open.

From this directory run `node check-endgame.js | diff - check-endgame.expected` and `node independent.js | diff - independent-expected.json`. The [independent audit](AUDIT.md) also gives the proof arithmetic.
