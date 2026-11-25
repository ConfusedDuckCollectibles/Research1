# Excluding 26-point spectral sets in the 13-element prime field cubed

**Theorem.** There is no spectral subset of size 26 in F_13^3. The general Fuglede problem in F_13^3 is not decided here.

The finite computations used below have exact, bounded input domains and are included. The strategy uses the elementary reductions of Fallon–Mayeli–Villano and the Radon method and last geometric step of Ferguson–Iosevich. The new step needed at p=13 is to use both members of a spectral pair to eliminate feasible projection profiles instead of trying to prove they do not exist.

## 1. Basic reductions

Set p=13. For B in F_p^3, let D(B) be its projective directions of nonzero differences, and U(B) the projective normal directions for which every one of the p parallel affine planes has |B|/p points of B. Suppose (E,A) is a spectral pair with |E|=|A|=2p. Its character matrix, with entries zeta_p^(e·a), has both orthogonal rows and columns.

For each nonzero difference of two frequencies, the sum of the corresponding character over E vanishes. If the p plane occupancies are c_0,...,c_(p−1), then the polynomial sum c_i X^i vanishes at zeta_p. Because its degree is at most p−1, it is a constant multiple of 1+X+...+X^(p−1). Thus the occupancies all equal two. Likewise for the transposed pair. Hence

    D(A) ⊆ U(E),    D(E) ⊆ U(A).

Both sets have at most two points on a line and at most p points on a plane. These are the standard reductions [1, Corollaries 1–2 and Lemma 2; 2, Lemma 8.2]. For the line bound, fix any line direction d belonging to E. Two points of A lie in a common plane normal to d by the pigeonhole principle (there are p planes and 2p points). Their nonzero difference is perpendicular to d and gives a uniform normal for E. The original line lies in one of the planes with that normal, which has exactly two E-points. For the plane bound, if more than p points of E lie in one affine plane, their set of differences determines every direction in that plane (p parallel lines of each direction). Then U(A) would contain a full projective line. A full projective line cannot be contained in U(B): the p+1 affine planes through any chosen line of common annihilator direction would all have intersection size two. Counting multiplicities of B in the pencil gives 2p+p k, where k is the number on the common line, instead of 2(p+1). This is impossible for p=13.

In particular neither difference direction set contains a full projective line. A set with no three collinear points will be called a cap.

## 2. Radon identity and the two-layer exclusion

Given a projective line of normal directions W, project B along its one-dimensional annihilator to F_p^2 and denote its fibre count by f(x,y), so f is in {0,1,2}. For each of the p+1 projective directions of W, let its profile be the sum of f on the parallel lines of that direction. Suppose p+1−L are designated uniform, hence have constant profile 2. Choose coordinates so one of those directions is represented by the form x. The other L profiles can be written g_i(y−s_i x), with distinct finite slopes s_i. Counting the p+1 lines through (x,y) gives the exact integer identity

    sum_i g_i(y−s_i x) = 2(L−1)+p f(x,y).                   (1)

If L=1, this identity makes f constant along each line in the only bad direction. Because every plane occupancy is at most p, we have f in {0,1} in that case. Its total mass 2p now implies B lies in exactly two parallel affine planes, each containing p points.

Neither member can lie in two parallel planes (compare [1, Lemmas 3–4; 2, Lemma 9.1]). Indeed its two p-point cap slices each determine at least p tangent directions, since the partition of all binom(p,2) pairs has at most (p−1)/2 pairs in each direction. The union therefore consists of exactly p directions in the layer plane, and every one is uniform for the other member. Applying the preceding L=1 argument to that member puts it in two parallel planes as well.

Assume E and A are each in two parallel planes. For one of the p directions d determined by the two E slices, every slice comprises (p−1)/2 disjoint pairs parallel to d and one leftover point: the pair-count bound is attained in each direction. For a normal in U(E) perpendicular to d, all pairs lie wholly in planes of occupancy two. The two singletons must therefore lie in the same plane, since counts on that normal are all even. Their cross-layer difference is not parallel to d. Thus there is at most one such projective normal perpendicular to both directions. The p projective lines d^perp are part of a pencil of p+1 lines; the last line has at most p U(E) points. Hence |U(E)|≤2p. But A's two p-point slices also determine p in-plane directions. After putting the layers at heights 0,1, their cross-layer directions correspond bijectively to the planar difference set Q−P of the two slices. This has more than p elements: if it had p, all translates Q−x, x in P, would coincide, so Q would be invariant under some nonzero translation. Its size p would make it an affine line, contradicting the cap condition. Therefore |D(A)|≥2p+1, impossible.

## 3. Exact profile completeness

Work modulo p in (1). Represent each profile as the unique polynomial of degree at most p−1 with its listed arguments as residues. Since each substituted term has total degree at most p−1, the bivariate polynomial identity on F_p^2 is a coefficient identity: such a polynomial has individual degrees at most p−1 and evaluation is injective. With c_(i,j) the coefficient of degree j of g_i, the homogeneous coefficients of each degree j>0 give

    sum_i c_(i,j) s_i^e = 0        for e=0,...,j.

All binomial coefficient factors are nonzero, as j<p. Vandermonde rank implies every coefficient with j≥L−1 vanishes. For 1≤j≤L−2 the full solution space is

    c_(i,j) = w_i R_j(s_i),
    w_i = 1 / product_(k ≠ i)(s_i−s_k),
    degree(R_j) ≤ L−j−2.                                  (2)

Indeed, the given form has the required vanishing sums because the coefficient of X^(L−1) in the Lagrange interpolation formula for a polynomial of degree at most L−2 is zero. The dimensions are equal, L−j−1. The remaining identity is sum_i c_(i,0)=2(L−1).

There cannot be nine or more good directions on a projective line. Designate any nine of them, giving L=5. The polynomials in (1) have degree at most three, and their minimum masses for nonconstant polynomials, with values interpreted as least nonnegative residues, are respectively 78, 52 and 39 at degrees 1, 2 and 3. These elementary finite minima are checked by the included Python program, by evaluating all degree d polynomials. Every ordinary profile has total mass 26, so its residue polynomial is constant. As the integer values lie in [0,13], the possibilities are the constant 2 or a profile with two entries 13 and all others 0. If k of the five profiles are constant 2, reducing (1) gives 2k=8 modulo 13, forcing k=4. There is exactly one bad direction; using the L=1 argument produces the excluded two-layer set.

The exact checks for 5–8 good directions are recorded below. Here M is the ordinary cap on the plane occupancy, strictly less than p. The checks apply to *exactly* L bad directions; no constant bad profile is permitted. Their domain and coverage are as follows.

* At least one good direction exists in every case in the table. Send it to infinity (form x). Every bad mask is then a subset of F_p. Affine maps on its slopes, s -> a s+b with a≠0, are realized by invertible linear maps on the plane. Exhaust all masks, retaining each affine orbit's minimum integer bit mask. We do not assume that the omitted infinity choice is unique.
* Enumerate all M-bounded profiles of total mass 26 and polynomial degree at most L−2 by interpolating every admissible word of L−1 values at 0,1,...,L−2. The other values and the mass are checked directly. Discard constants.
* Coordinate translations (x,y)->(x+x0,y+y0) independently translate the arguments of the first two distinct slopes. For a profile of positive degree d<p, the next-to-leading coefficient after translation by h is c_(d−1)+d c_d h; so there is a unique h centering that coefficient. Restrict the first two profiles that way, including when their degrees differ.
* Descend over j=L−2,...,1. At level j choose any L−j−1 slopes, enumerate all their presently possible coefficients, and interpolate R_j from (2). The choices at the other slopes are uniquely determined. Restrict each profile catalogue to compatible coefficients and recurse. At j=0 the remaining constant choices are checked against (1) modulo p. Recompute f using the integer (1); reject nonintegral values and values outside {0,1,2}. For every surviving f the program directly recomputes all p+1 line profiles to check agreement. These checks only use necessary conditions for projection, not any assumption that a surviving profile can lift to a spectral set.

| Bad L | Cap M | Affine masks | Profiles per unconstrained column | Tuples after positive degrees | Survivors |
|---:|---:|---:|---:|---:|---|
| 6 | 8 | 14 | 117 | 3 | 3; one weight-2 projected point each, peak occupancy 6 |
| 7 | 7 | 14 | 585 | 12 | 0 |
| 8 | 6 | 10 | 3497 | 290 | 84; peak 6, and good masks in the two projective classes below |
| 8 | 5 | 10 | 1612 | 20 | 0 |
| 9 | 4 | 7 | 3562 | 72 | 0 |

For L=8, M=6, the canonical projective bit masks of the good directions are 95 and 111. Here a bit in position t denotes the finite slope t; position 13 denotes infinity; canonicalization takes the minimum bit mask over every invertible projective linear map. The column “tuples” counts all full centered profiles that pass the positive-degree tests, before the constant and exact grid test. All runs are exhaustive and uncapped.

## 4. Coupling the two members

One additional exact planar check gives the following facts about caps in F_13^2:

1. An eight-point cap cannot determine at most eight directions.
2. A five-point cap cannot determine at most five directions.
3. The direction set of a six-point cap determining exactly six directions has projective canonical mask 189, distinct from 95 and 111.
4. A cap in which parallel chords always have the same difference vector up to sign has at most seven points.

These are small affine-invariant tests. Send an ordered noncollinear triple to (0,0),(1,0),(0,1). Extend by lexicographical search through all compatible points. For the first three statements extend each direction set to the specified size and include all masks with the triangle's three directions, requiring all later chords to lie in the mask and never allowing collinear triples. For the fourth, instead record the only permitted difference vector up to sign in each used direction and reject conflicts. The program also finds a seven-point control. Thus neither classified source data nor a solver's absence of results is used.

Since there are at most eight uniform normals on a projective line, a nine-point planar slice is impossible: its 36 pairs give at least nine directions, each contributing at most four pairs. The first planar check excludes an eight-point slice as well. Thus the affine plane cap is at most seven. The L=7, M=7 gate eliminates exactly seven good directions on any projective line.

Now fix any normal v that is **nonuniform for E**. It does not belong to D(A) since D(A)⊆U(E). Hence projection of A along v is injective. Its line of dual normals v^perp cannot have eight good directions: the L=6 gate says that every such projection has a weight-two fibre. It cannot have seven either. Therefore the directions of each E slice normal to v, lying in D(E) and in v^perp, are contained in a set of at most six directions, namely U(A)∩v^perp. A seven-point cap determines at least seven directions, so these slices contain at most six points. Uniform normals have slices of size two. Apply symmetry to get the plane cap six for both members.

If some E slice normal to v has six points, their directions must number exactly six: at most six by the preceding argument, and not at most five by the five-point cap check on a subset. The six directions then equal U(A)∩v^perp and their projective mask is 189. But the L=8, M=6 gate says that this good-direction mask belongs to orbit 95 or 111, a contradiction. Hence both plane caps improve to five.

The L=6 gate requires peak six, excluding eight good directions at the improved plane cap. Seven were already excluded; the L=8, M=5 gate excludes six. Now every five-point plane slice would require a planar cap with at most five directions, which is impossible. Thus the plane cap is four. The L=9, M=4 gate excludes five good directions, and each projective line meets U(E) and U(A) in at most four points.

## 5. Final pair count

Let D=D(E). Each projective line meets D in at most four points, so a pencil of 14 lines through an element of D gives |D|≤1+14·3=43. Fix d in D. Every secant in direction d contains exactly two points; let r_d be the number of such secants.

Any two of the parallel secants have four endpoints in an affine plane. These four points form a parallelogram: their six edges have at most four directions and no direction occurs more than twice, so at least two pairs of opposite edges are parallel. The two secant displacement vectors therefore agree up to sign. Orient all secants with one displacement v and write them as {b_i,b_i+v}.

Project the b_i along v. This planar set of r_d points is a cap, since three collinear projections would put six points of E in one affine plane. Chord vectors of this projected cap in the same direction agree up to sign. To check this, take two lifts w and w'=lambda·w+mu·v of parallel projected chords, with lambda≠0. The two secant pairs for the first chord use all four possible directions in their linear plane, namely [v], [w], [w+v], [w−v]. The second chord must therefore satisfy

    {mu/lambda,(mu+1)/lambda,(mu−1)/lambda} = {0,1,−1}.

Comparing the sums and sums of squares gives mu=0 and lambda²=1. The fourth planar check therefore gives r_d≤7. Finally

    325 = binom(26,2) = sum_d r_d ≤ 43·7 = 301,

the contradiction.

## References

1. T. Fallon, A. Mayeli and D. Villano, [*The Fuglede conjecture holds in F_p^3 for p=5,7*](https://arxiv.org/abs/1902.02936), [*Proceedings of the AMS*](https://doi.org/10.1090/proc/14750).
2. S. J. Ferguson and J. Iosevich, [*Fuglede’s conjecture holds in F_p^d for p=2,d=8 and p=11,d=3*](https://arxiv.org/html/2610.30712v1), especially §§8–9. Their proof treats p=11 and precedes the present profile bootstrap. In p=13 even realistic profiles in some high-uniformity pencils do survive; the coupled-member exclusion here is needed.
3. R. D. Malikiosis, [*A linear programming approach to Fuglede's conjecture in Z_p^3*](https://link.springer.com/article/10.1007/s43670-023-00078-2), an exclusion near cardinality p² that does not cover 2p at p=13.
