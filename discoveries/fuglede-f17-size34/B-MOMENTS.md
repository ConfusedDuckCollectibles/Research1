# Appendix B. Rich sections and polynomial moments

The order of the geometry is A3 → B1 → B2 → B3 (reciprocal richness). Once both sides are rich, B4 gives degree-five moment vanishing by spectral uniformity; the separate evaluation-rank argument in B3 §2 then applies. See [Reproduce](REPRODUCE.md) for the four decisive arithmetic gates in this appendix.

- [B1. Boolean Radon equations and the zero-rich-family gate](#b1)
- [B2. At most two rich families are impossible](#b2)
- [B3. Exactly three rich families are impossible; evaluation rank](#b3)
- [B4. Marked five-plane obstruction and degree-five moments](#b4)

<a name="b1"></a>

## B1. Boolean Radon equations and the zero-rich-family gate

### Inputs and exact Boolean equations

Use the prior [A3](A-AFFINE.md#a3) and the size-34 K7 theorem: the two members have line cap2, plane cap6, and at most six uniform normals on a projective line. At least one member has a rich plane. Every rich section has the hexagonal six-direction palette, forcing the opposite mate's pencil of exactly six regular directions. That quotient is injective since its annihilator is the normal of the rich plane, nonuniform on the rich member; so the opposite mate has no secant along it.

If that opposite member is not rich, all lines of this Boolean quotient Q have at most four points (each is the projection of an affine preimage plane). Normalize the six directions to {0,1,2,9,16,∞}; slope and normal coordinates differ by projective dual identification, as in the previous notes. The twelve bad finite slope labels are S=(3,4,5,6,7,8,10,11,12,13,14,15). Each bad profile g_s(t) has mass34, values 0..4 and is nonconstant. Writing w=indicator of Q, the exact identity is
    Σ_s g_s(y−sx)=22+17w(x,y).
The real variance needed is
    Σ_s q_s=510,    q_s=Σ_t(g_s(t)−2)^2.
Indeed across all 18 direction families Σ g² =18·34+34·33=1734, since every ordered distinct pair belongs to one direction. Subtract 4·18·17=1224. The six good profiles contribute zero variance. This use of 510 is valid because w is Boolean, not for a general weighted quotient.

Over F17 each g_s is its unique reduced polynomial. Comparing homogeneous terms gives for every j>0 the conditions Σ_s c_(s,j)s^h=0, 0≤h≤j. With W_s=(∏_(t≠s)(s−t))^−1 these say c_(s,j)=W_s R_j(s), deg R_j≤10−j for j≤10, and c_(s,j)=0 for j≥11. This is the Lagrange Vandermonde kernel identity (the Σ_s W_s P(s)=0 condition for polynomials of degree≤10), so full degree is at most ten. Constants must sum to22 modulo17.

### Complete exact local inventory

Center each nonconstant profile so its coefficient at d−1 vanishes; a translate and scale of its argument preserve mass, cap and variance. Enumerate centered forms by splitting into g(t)=E(t²)+t O(t²), with E of degree k=d/2 and fixed top coefficient when even, or of degree≤k−1 with O's degree-k top fixed when d=2k+1. Centering forces the other parity's highest term to vanish. Enumerate every choice of the remaining O coefficients and work in leading coefficient representatives modulo dth powers; scale the input to obtain each other lead.

At each j=1,…,8 the unique possible signed integer difference of g(j),g(−j) must be within [−4,4] and have residue 2j O(j²). If its absolute value is a_j, the pair mass is a_j+2b for 0≤b≤4−a_j and the field midpoint is a_j/2+b. At j=0 the domain is exactly 0..4. After subtracting the fixed leading even term if present, any k midpoint locations determine the still-unknown even polynomial of degree≤k−1. Choose the k smallest domains (distinct squares), interpolate, then check every remaining point and exact mass34. If k=0 there is no unknown even polynomial. This is full enumeration of all possible capped centered profiles, not a random or numerical feasibility test.

Counts across all leading coefficients are 0 for d≤6, then 16,20,112,2096 for d=7,8,9,10 respectively. In degree7 all q_s=34. In degree8 the only physical leads are ±1 and ±2; in degree9 all 16 leads have seven profiles. In degree10 there are 128 profiles for each square lead and 134 for each nonsquare lead.

### Complete one-palette coupling, all possible top degrees

Let d be the highest local degree. If d=7 all rows are degree seven, and energy is 12·34=408, contradicting 510. If d=8, the weighted quadratic lead must be nonzero and have each nonzero actual value W_s R_8(s) in {±1,±2}. Enumerate all 17³ quadratics, excluding zero; no such polynomial exists. The code quotients only by common scaling (8th powers {±1}, so arrange first nonzero actual lead in {1,…,8}); a lower-degree slot is not disallowed merely because its octic lead is zero.

For d=9, the weighted leading coefficient is a nonzero binary linear form. Enumerate its root in all 18 projective positions, including infinity (=constant in affine coordinates), with scalar normalized to 1 by common dilation, since ninth powers are a bijection. The weighted octic coefficient is quadratic. Translate the quotient plane to make it vanish at two selected nonroot slopes, leaving γ(s−s_i)(s−s_j) with all 17 values of γ. Each actual nonzero lead fixes the unique shift of any centered shape with that lead. A zero top lead is allowed and, if its octic lead is admissible, every centered octic and every shift is included. If both leads vanish the exceptional row is septic, and every centered septic of all 16 leads and all of its translates is included. No lower degree is possible by the cap-four inventory.

For d=10 every bad slot has nonzero decic lead since R_10 is constant; normalize it to one of the square classes {1,3}. Simultaneous translation kills the weighted affine degree-nine coefficient at every slot (because R_10 is constant, a quotient translation contributes an arbitrary affine slope function). All local rows are thus centered.

Up to d−2, enumerate coefficient-value buckets with q_s min/max. For each j, select 11−j pivot labels (the smallest current domains) and try all their possible normalized coefficient values. The remaining values follow by exact interpolation and must have nonempty buckets. Prefix and bucket variance intervals discard only choices that cannot sum to 510. The constant sum and full integer Radon equation would be checked on any surviving tuple, but the relevant stages vanish earlier. For d=10, 9,826 degree-eight prefixes give 598 patterns; 24,114 degree-seven prefixes give 40, and all 40 degree-six prefixes fail. For d=9, 306 root/shift gauges were tested, 162 had local options, and all 9,804 degree-seven prefixes fail.

<a name="b2"></a>

## B2. At most two rich families are impossible

### Complete exact Boolean atlas with up to two elevated directions

If a forced quotient has at most two rich families, pad to a marked two-set of its 12 exceptional (bad) slopes. All other ten have cap4, and the marked two have cap6. A projectivity stabilizing the six good directions preserves Boolean mass and transports the marked two-set. The exact 12-element stabilizer acts on all C(12,2)=66 pairs and gives 9 orbits, of sizes 6,6,6,6,12,6,6,12,6. Their representatives in this normalized slope chart are {3,4}, {3,5}, {3,6}, {3,7}, {3,8}, {3,10}, {3,11}, {3,13}, {3,15}, in that order. A padded slot is not required to be rich. This argument covers 0,1,2 elevated directions, without treating the orbits as independent coordinate relabels without geometry.

Use the exact cap4 and cap6 centered local catalogues from the even/odd interpolation method in [B1](B-MOMENTS.md#b1). All nonconstant profiles have d≥6, and all ten low cap4 slots have d≥7. The inventories for d=6,7,8,9,10 are
    cap4: 0,16,20,112,2096;
    cap6: 16,48,272,3008,53264.
The finite choices include all shifts from centered forms; for each slot its cap comes from membership in the selected pair. The Boolean total variance is 510. In degrees ≤7, these two caps have per-profile variance upper bounds 34 (ten low slots) and 68 (two elevated slots). Thus all top degrees ≤7 are impossible by energy≤476. For top8 all 17³ quadratic weighted leading coefficients are tested up to simultaneous scaling; the two elevated labels use their cap6 allowable lead set, while the remaining ten use {±1,±2} or admit a zero leading coefficient. None is feasible, even before checking the next coefficient.

The top9 and top10 normalizations are exactly those in [B1](B-MOMENTS.md#b1) and choose any two nonzero top labels for the translation. Top9 includes every legal degree≤7 row at a slot when both its nonic and octic coefficients vanish, including centered sextics and all their shifts at elevated slots; when only the nonic vanishes, it includes every legal octic lead and prescribed/all shifts. Top10 is nonzero at all twelve, with all rows centered simultaneously. All positive coefficient stages test the exact normalized-degree interpolation conditions and prune only when attainable variance intervals exclude 510; no approximated variances, Boolean solver, or presumption about profiles' rootlessness is used.

Aggregate across 9 geometric orbits:
* Top9: 2,754 root/shift gauges, 1,638 have all local option sets nonempty. Exactly 1,007,154 degree-seven interpolation prefixes yield 176 patterns. All 198 resulting degree-six prefixes fail.
* Top10: 18 scale gauges; 88,434 degree-eight prefixes yield 88,124 patterns; 5,538,422 degree-seven prefixes yield 99,787 patterns; all 117,762 degree-six prefixes fail.
* No complete profile tuple survives, even before constants and exact Radon reconstruction.

The geometric implication follows by fixing a rich normal of E and taking its forced quotient of A along that normal. Every projected line has as preimage a genuine affine plane parallel to that normal, with the same number of points by injectivity. Its projective normal is therefore orthogonal to E's marked normal. Distinct line families in the quotient have distinct preimage projective normals. The construction is symmetric, giving the asserted minimum degree on both partitions.

<a name="b3"></a>

## B3. Exactly three rich families are impossible; evaluation rank

### 1. Exact three-rich closure

#### Prior input and notation

Use [A3](A-AFFINE.md#a3), [B1](B-MOMENTS.md#b1), [B4](B-MOMENTS.md#b4) and [B2](B-MOMENTS.md#b2). Both members have line cap 2, plane cap 6, and each uniform pencil has ≤6 directions. Each rich normal (a normal of an affine plane with at least five points) induces on the opposite mate the single fixed G6 orbit of uniform pencils. Projection along the axis injects into a 34-site Boolean planar Q: the axis is nonuniform on the rich member and cannot be a difference direction of the mate. Every projected line has its genuine parallel-plane preimage, with exactly the projected count; cap6 is valid. The previous computation excludes 0, 1 or 2 rich line families in this quotient. Thus to strengthen it, all three marked families in an exactly-three quotient really do have maximum occupancy ≥5; this is not a padding argument using cap6 rows that might remain low.

Use good slopes {0,1,2,9,16,∞}. The 12 exceptional finite slopes S are (3,4,5,6,7,8,10,11,12,13,14,15), and every exceptional line-count profile g_s is nonconstant. Its weights are in [0,4] at unmarked slopes and [0,6] at marked slopes, with max≥5 at each marked slope. For 34 Boolean sites the exact Radon identity is Σ_s g_s(y−sx)=22+17w(x,y); the variance Σ_s q_s=510 where q_s=Σ_t(g_s(t)−2)^2. Indeed the total squared occupancy across all 18 directions is 18·34+34·33 (each distinct ordered site pair meets in one direction); subtract 18·17·4. Good profiles contribute zero variance. Use of 510 depends on Boolean, not general weighted charts.

Let W_s=(∏_(t≠s)(s−t))^−1 in F17. The positive-degree coefficients of the reduced profiles are c_(s,j)=W_s R_j(s) for j≤10 with deg R_j≤10−j; no coefficients above degree ten. To see this, the Radon identity mod17 gives Σ_s c_(s,j) s^h=0 for 0≤h≤j by comparing homogeneous terms (j≤16 has invertible binomial coefficients). The Lagrange kernel gives the claimed degree bound, with trivial kernel when j≥11. Constants sum to 22 mod17.

The 12-element physical PGL stabilizer of the six good slopes has 19 orbits on the 220 unordered marked triples, sizes 12 except one orbit of size 4. The script enumerates all matrices with first nonzero entry 1, verifies coverage, and transports only genuine projective slope relabelings (with actual GL2 representatives on the plot). It is safe to rotate to a representative because line masses and maximum occupancies are unchanged by a linear plane map.

#### Complete local inventories and higher stages

The accompanying verifier reconstructs the previously validated centered local catalogues by the even/odd interpolation method explained in [B1](B-MOMENTS.md#b1). Briefly center the degree-d profile by killing its coefficient of degree d−1, split as E(t²)+t O(t²), and for one leading representative per dth-power class enumerate all free odd coefficients. For each ±t pair the unique signed difference representative in [−8,8] must have magnitude ≤cap. The legal even-midpoint domains then have size cap+1−|difference| (t=0 gives cap+1); choose the smallest number needed to recover the even polynomial, interpolate after subtracting any fixed even top term, and check the rest and exact mass34. Dilating t produces every actual leading coefficient in that class. The code checks the exact degree, center, cap, variance and sum for each retained profile. A genuine degree6,7,8,9,10 has respectively 0,16,20,112,2096 centered cap4 profiles (over all actual leads). The corresponding cap6 counts are 16,48,272,3008,53264, of which 16,32,252,2896,51168 have peak≥5. All actual translates of in-catalogue profiles are considered.

Top degree ≤6 is impossible with nine low slots. For top degree 8, the weighted quadratic R_8 is tested over all 17³ possibilities except zero; all nonzero leads must belong to their exact degree8 lead set. A zero at a slot is allowed and treated only as a lower-degree possibility. There are zero surviving leading quadratics across all 19 orbits.

For degree 9, R_9 is either a nonzero linear polynomial or a constant (root at infinity). Common argument dilation normalizes its scale uniquely because gcd(9,16)=1. Common translation of the chart replaces t by t+L(s) with L affine, so the degree-eight coefficient can be killed at any two nonroot slopes. Enumerate the remaining λ(s−s_i)(s−s_j), λ in F17. At a zero top slot cover all legal lower degrees and translates; there is no implicit rootless assumption. For degree10, R_10 is constant nonzero at every exceptional label, hence one of two square classes after argument dilation. Common translation kills all of R_9 because R_10 is constant and R_9 affine. For each lower coefficient stage solve R_j from 11−j pivot labels (fewest current possible residues) and filter the remaining rows exactly. Prune only if the attainable real variance interval cannot contain 510. No maximal-profile assumption is needed other than the prescribed peak thresholds. If all positive coefficient stages passed, the code would check the constant sum and exact integer Radon reconstruction. Here no such tuple reaches the final gate.

High-degree counts, totals over the 19 normalized marked cases:
- Top9: 5,814 gauges, 3,534 with all local option sets nonempty; 2,208,750 degree-seven prefixes give 7,160, and all 14,591 degree-six prefixes fail.
- Top10: 38 gauges; 186,694 degree-eight prefixes give 186,694; 16,235,478 degree-seven prefixes give 1,357,324; 1,483,634 degree-six prefixes give ten; all ten degree-five prefixes fail.
- No final integer reconstruction survives. This is a complete test of an actual three-rich chart, not a sample of high-stage gauges or a timeout exclusion.

#### The separately compressed degree-seven equality case

At a low slot the inventory is one centered septic per nonzero leading coefficient, q_s=34, with centered fifth coefficient κ c_7³ for κ=7. At a rich slot there are two centered septics per lead, q_s=68, with κ=6 or12; or one centered sextic per lead, also q_s=68. The total is exactly 9·34+3·68=510. At a nonzero septic slot, translation gives c_6=7 c_7 u and adds binom(7,2)c_7 u² to its centered c_5. Eliminating u over F17 gives
    T(s)=R_7(s) R_5(s)−15 R_6(s)^2=κ W_s² R_7(s)^4.
Here T has degree≤8. The nine low slopes determine it by interpolation for each R_7. At a rich zero of R_7 the sextic leading coefficient must be nonzero, so T(s)=−15 R_6(s)^2 implies T(s)/(−15) is a nonzero square. At a rich nonzero slot T must satisfy the displayed equation for one of κ=6,12. This is a necessary prefilter without enumerating R_6.

Normalize R_7 monic by common argument dilation (gcd(7,16)=1), allowing every degree ≤3. Exactly 57,399 normalized choices across all marked orbits are nonzero at the low slots, and 143 pass the T prefilter (73 with a marked zero). Kill R_6 at the first two low slopes by common chart translation. Enumerate all 17³ residual quadratics after dividing those two factors out, for 702,559 checks; 377,704 satisfy the root-square equalities. Low slots determine R_5, of degree≤5, from its first six low values. With known c_7 and c_6 at a nonzero slot its shift is unique; at a zero c_7 slot c_6 and c_5 similarly determine the sextic shift. Check the remaining c_5 coefficients exactly. There are 264 resulting full high-prefix profile assignments. Every one fails the degree-four condition. This finishes the low branch, with no integer candidate.

#### Geometric consequence

Each distinct rich projected line family gives an actual rich normal orthogonal to the original rich normal of the other member. The exact exclusion gives degree≥4 in both parts of the reciprocal rich-normal graph. There are no 4-cycles because two distinct projective normals have only one common orthogonal normal. Starting at one vertex, four neighbours and their ≥3 additional, pairwise disjoint neighbours imply at least 13 vertices in its part, and likewise in the other part.

If one part in a component has exactly 13, the other has at most 13 by counting distinct pairs being covered by vertex neighbourhoods of size≥4 (78≥6 times its size). Equality forces 13 on both sides, all degrees 4, and every pair covered. The chosen points and reciprocal chosen lines then form a projective subplane of order 3 in PG(2,17). For an elementary contradiction, normalize a quadrangle of its points to the affine square (0,0),(1,0),(0,1),(1,1). The standard horizon is the join of intersections of the opposite sides, hence belongs to the subplane, as do the horizontal and vertical infinity points and the infinity point of slope −1 from the square diagonal. For any subplane point (x,0), go vertically to (x,1), then in slope −1 to (x+1,0). Incidence closure of joins and intersections holds, forcing all 17 points (x,0) on a subplane line, impossible. Thus each part has at least 14.

### 2. Coupled cubic rank obstruction for the sixth moments

This is analytic and does not need the newly increased minimum degree. Its inputs are size 34=0 in F17, both members' mixed moments through total degree five ([B4](B-MOMENTS.md#b4) and [B1](B-MOMENTS.md#b1)), and spectral direction–uniformity duality. Use the natural coordinates on the paired spaces; notation M_6^X(t)=Σ_(x in X)(t·x)^6.

Let V_j(X) be the image of affine polynomial evaluation through degree j in F17^X, and K_j(X) the kernel. Put r_X=dim K_3 and q_X=dim K_2. The dimension of V_3 is 20−r_X. For the ten multi-indices |α|=3, put
    C_X(α,β)=Σ_(x in X) x^(α+β).
This is a middle catalecticant of M_6^X up to invertible normalizations. Mixed moments through degree five, including degree zero, say that in the Gram matrix for all degree≤3 evaluations the blocks using any row or column of degree≤2 are zero. Its rank is therefore rank C_X, and is also the rank of the dot product restricted to V_3. Since the ambient dot product on 34 coordinates is nondegenerate, the radical of the restriction has dimension at most 34−dim V_3. Thus
    rank C_X >= 2(20−r_X)−34 = 6−2r_X.

For the upper bound on rank C_E, let L_A be the vector space of cubic moments of all **weighted** functions on A whose affine moments through degree two are zero: taking weighted sums α→Σ a_x x^α for |α|=3. Its dimension is rank eval through degree3 minus rank eval through degree2, i.e. 10+q_A−r_A. Since M_6^E vanishes on D(A) (every actual difference normal makes E uniform), M_6^E(a−b)=0 for all a,b∈A, also for a=b. Multiply these equations by two independent such weight vectors, one for a and one for b, and sum. Every bidegree term vanishes except (3,3), because the remaining terms have one variable group of degree at most two. The (3,3) coefficient is
    −binom(6,3) Σ_(e in E) (e·a)^3 (e·b)^3,
so scaled L_A (diagonal scaling by the nonzero multinomial coefficients of degree three) is totally isotropic for C_E. A symmetric form on a 10-dimensional space admitting a totally isotropic subspace of codimension d has rank at most 2d by writing its top-left block as zero in adapted coordinates. Consequently
    rank C_E <= 2(r_A−q_A).
Swap A and E for the dual. This gives the asserted sum inequality. In the boundary sum =3, both q's are zero and both matrices reach their respective bounds. Each V_3(X) then contains its dot-product orthogonal complement (the radical reaches ambient orthogonal-complement dimension).

The same middle-block argument without using the paired lower bound applies to any homogeneous sextic F vanishing on D(A): its middle bidegree coefficient matrix has rank at most 2(r_A−q_A). In particular, if the cubic evaluation is injective then that matrix is zero, forcing F=0 because every degree-six monomial splits into two of degree three and all relevant binomial coefficients are nonzero in F17. Hence cubic-unisolvency of either mate would force the other's full degree-six moments to vanish. This fully includes the hexagonal binary sextic contribution at infinity; it does not set that contribution to zero by assumption. If M_6^E actually vanishes then r_E≥3, since V_3(E) is entirely isotropic. The same coupled derivation in F_p^3 with size 2p and degree-five moment assumptions gives sum ≥20−p+max(q_A,q_E) when p>6, but no degree-five premise for another prime is claimed here.

<a name="b4"></a>

## B4. Marked five-plane obstruction and degree-five moments

### 1. Exact marked-direction rigidity in three dimensions

Let E be a 34-point cap in F17^3, with plane cap six, with a planar five-point cap that extends to the centered hexagon and for which all differences in E parallel to that plane use only its six directions. This includes every rich-plane member of a putative size-34 spectral pair. The target **proved below**: no nonzero homogeneous polynomial of degree five vanishes on D(E). Consequently none of smaller degrees vanishes either, since multiplication by a nonzero homogeneous linear form preserves nonzero polynomial and vanishing.

Take the marked plane z=0 and normalize its five points to -Q for Q={(1,0),(-1,0),(0,1),(0,-1),(1,1)}. Indeed the earlier cap atlas says the five extend to a centered hexagon {±u,±v,±(u+v)} and all deleted vertices are equivalent by linear symmetries of this hexagon. A six-point section just supplies such a choice of five and is otherwise not needed. Let G be their six distinct planar directions. If F is homogeneous degree five vanishing on D(E), its binary restriction F(x,y,0) must vanish at six distinct projective points and hence is zero. Thus F=z C, and f(x,y)=C(x,y,1) is a nonzero polynomial of total degree≤4. Evaluating at affine z=1 retains all homogeneous coefficients of C, so nonzero is automatic.

At least 28 points have z≠0 by the plane cap six (and at least 29 for an exact five-point plane). With 16 nonzero levels there is a pair at one common nonzero level. After z dilation and the shear adding a multiple of z to x and y, normalize that pair to e0=(0,0,1) and e1=(dx,dy,1), without changing the five marks. Its nonzero horizontal difference has a direction in G. Exchanging which endpoint is e0 changes the sign of the displacement (and again shears by its full displacement), and swapping x,y preserves the marks. There are 6·16=96 raw displacement vectors and exactly 32 representatives under those operations. Nothing else is taken to be equivalent.

For every exterior point e=(x,y,t), t≠0, the five points ((x+u)/t,(y+v)/t) for (u,v)∈Q must all be roots of f, using the difference e−(-Q,0). For t≠1 the two extra roots x/(t−1),y/(t−1) and (x−dx)/(t−1),(y−dy)/(t−1) come from e0 and e1. Hence:
1. Form the vector space of degree≤4 polynomials satisfying the ten (possibly duplicate) root equations of e0,e1. Its dimension is 5,6 or7.
2. For every one of 15·289 possible e with t≠0,1, restrict that space by the seven new equations and enumerate its nonzero projective kernel. Count each resulting polynomial once for that point. This avoids enumerating all raw polynomials.
3. At least 28−6=22 points occur off both heights 0 and 1. Accordingly retain every (displacement,projective polynomial) whose possible count is ≥22, and reconstruct its complete exterior pool including height1. The exact attainable counts jump from ≤20 to 30.
4. On the pool, discard any pair violating either the common G when its height is equal, or f((x−x')/(t−t'),(y−y')/(t−t'))=0 when its height differs. Also discard pairs collinear with a fixed anchor or a marked point. Search exact subsets with at most six per height and no further collinear triple, using greedy independent coloring as a clique upper bound. Pairs on the same level do not introduce a spurious root of f; equality of heights is handled separately. This is only a necessary feasibility problem; disallowing any 28 exteriors is enough.

The verifier checks all 138,720 third-point tests over the 32 displacements and counts 434,232 projective occurrences. There are 70 surviving (displacement,form) cases; their exact subsets require 21,719 search nodes and have global maximum 6 outside points, counting the two anchors. No spectral assumption beyond the direction and cap hypotheses is used in this rigidity proof. Greedy coloring classes are independent sets in the compatibility graph, so any remaining clique has at most the largest available color index. The recursion still explicitly rejects forbidden triples and overfull heights.

For a spectral pair, the opposite homogeneous degree-$j$ moment form vanishes on $D(E)$ for $1\le j\le5$, because every field value occurs twice along each such direction. The statement above excludes every such nonzero form, by multiplication by a homogeneous linear form to the required power. Polarization gives every mixed moment through degree five on both now-rich mates, since the pair can be interchanged. The degree-zero moment is 34, zero in the field.
