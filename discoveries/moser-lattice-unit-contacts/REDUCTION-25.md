# The 25-point frontier and a near-Moser reduction

Write u(n) for the largest number of unit pairs in a set of n distinct points in the real plane.

**Status audited on 25 November 2026.** Alexeev–Mixon–Parshall, [The Erdős unit distance problem for small point sets, arXiv:2412.11914v2](https://arxiv.org/html/2412.11914v2), Theorem 1, established exact values through 21, and reported 72≤u(25)≤78. Their original graph6 list is in the arXiv ancillary file; their paper gives detailed algorithmic procedures in §§2–4. The first-hand Artificial Wasteland project published certificate-based u(22)=60, u(23)=64, u(24)=68 in September–October 2026, not yet externally refereed: [u(22)](https://doi.org/10.5281/zenodo.22468977), [u(23)](https://artwaste.land/strata/the-sixty-fifth-distance/), [u(24)](https://doi.org/10.5281/zenodo.23272696) with corresponding [original Python checker and evidence](https://github.com/artwasteland/Maths/tree/main/unit-distance-24). The latter records 1,920 accepted slices and separately checked nonembeddability of the sole 24/69 survivor; it is not a list of all planar 24/68 extremizers. Adopting their 68, a random vertex deletion shows u(25)≤73 because any 74-edge graph has a 24-vertex induced graph with at least ceil(74·23/25)=69 edges. I found no later public exact 25 result. Our work does not independently rerun their enumeration, so the following is conditional on their 22–24 exact inputs.

The established lower 72 occurs in Engel et al., [arXiv:2406.15317v3, Table 2](https://arxiv.org/html/2406.15317v3). As a simple independent exact coordinate check, put A={(1,0),(2,0),(0,1),(1,1),(0,2)} and take all 25 basis tuples (a,b,c,d) with (a,b) and (c,d) in A. The distance formula in MOSER-CONTACTS.md gives precisely 72 edges (70 within fibers of the two projections and two between distinct fibers). This re-establishes the known lower bound, not a new construction record.

## Exact conditional conclusions

Suppose P is a 25-point set with 73 unit pairs. For any isometric copy of M containing at least 22 members of P, all 25 members are in that copy. If some copy of M contains at least 20 members, there is a short list of possible ambient arithmetic groups: the entire set is in that M or one of its two explicit index-three supergroups
    M + Z(1+eta)/3,       M + Z w(1+eta)/3.
This is a reduction, not a proof of u(25)=72 or 73, and does not imply the presence of a 20-point M core in an arbitrary extremal configuration.

## First gate: at most three outside

Put S=P∖M after normalizing the copy. By u(24)≤68, every point has degree at least five. Our exterior contact theorem says each point of S has at most three neighbors in M, so G[S] has minimum degree at least two. S cannot have size one or two. If it has size three, it is a triangle, and all three vertices have exactly three neighbors in M; deletion of S removes 3+9=12 edges, leaving 61 on 22 points, contradicting u(22)=60.

## Second gate: four or five outside

Here is the complete finite structural gate, run by verify-reduction.js. For r=4,5 enumerate every labelled graph F on the outside set with minimum degree at least two and at most u(r) edges. Reject K_4 (four planar equidistant points do not exist) and K_{2,3} (two centers have at most two common unit neighbors). For each outside vertex i enumerate the contact count t_i with the actual finite core, where max(0,5−d_F(i))≤t_i≤3. For every nonempty T⊆S test the exact necessary deletion inequality
    73 − sum_{i in T}(d_F(i)+t_i) + e_F(T) ≤ u(25−|T|).
Each t_i=3 proves that point is in (1/3)M (indeed it has three contacts with the infinite M). Propagate this membership over triangles (the last point is obtained from either endpoint by rotating their difference by ±pi/3, multiplication by eta or its inverse, which preserves the module) and over four-cycles (the last point is obtained by the opposite-vertex sum identity). That identity holds even without presumed noncrossing drawings: viewing two opposite vertices as centers, the other two distinct opposite points lie at the two common circle intersections and are reflections across the center line.

The exact count of labelled graphs/contact profiles after the deletion test and propagation, by (outside edges, still unassigned), is:
    r=4: (4,0) 3; (5,0) 18.
    r=5: (5,0) 12; (6,0) 210; (7,0) 480; (7,3) 60.
All 60 not yet propagated are the same configuration up to relabeling: the gem, a four-vertex path attached to a universal central vertex, with the two rim endpoints having t=3 and the other three vertices t=2. The checked labels include the shape of the path. The outer path lies on the unit circle about the center and steps successively by ±60°. A change of sign would retrace a vertex, impossible. Its endpoints A,E are therefore antipodal, with length |A−E|=2. The half-difference lies in (1/6)M and is a unit vector. In fact it belongs to (1/3)M by the denominator lemma below. Hence the midpoint and all the rim vertices belong to (1/3)M, since the rim positions differ from the midpoint by eta-powers of the half-difference.

**Denominator lemma.** For an integer tuple (a,b,c,d) with |a+b eta+c w+d eta w|/(2m)=1, the irrational norm coefficient vanishes, so ad=bc. The nonzero rank-one matrix with rows (a,b),(c,d) can be written with rows u(s,t),v(s,t), where (s,t) is a primitive integer vector and u,v are integers. Substituting the norm form gives
    (s²+st+t²)(3u²+5uv+3v²) = 12m².
Either quadratic form vanishes mod 2 exactly when its two entries are both even. The first factor is odd by primitivity, so both u,v are even. This proves a unit vector in (1/(2m))M belongs to (1/m)M; use m=3.

## The two possible index-three supergroups

All unit vectors of (1/3)M have a numerator tuple with P=54 and Q=0 and coordinates in [−10,10], since p=9 and p≥(a²+b²+c²+d²)/12. The script enumerates exactly 42. By residue class modulo 3, 18 are in M (class 0000), each special class 0011,0022,1100,2200 has three, and there are 12 further classes with one each. No unit residue has both of its two row blocks nonzero. Subtracting residue classes gives this count matrix for possible unit differences, in the class order 0000,0011,0022,1100,2200:
    18 3 3 3 3
     3 18 3 0 0
     3 3 18 0 0
     3 0 0 18 3
     3 0 0 3 18

Any off point whose t_i is at least two must be in one of the four special classes. The script also confirms that t_i≤1 can occur only for t_i=1 at the center of a seven-edge gem, whose rim endpoints already have three contacts. The center's class is the half-sum of the endpoint classes. If those are of different types (0011/0022 versus 1100/2200), the half-sum has two nonzero row blocks and cannot have even one M-contact. If they are opposite within one type, the half-sum is zero, incompatible with being off M. Thus that point too is in a special class. For r≤5 the outside graph of minimum degree two is connected (each component needs at least three vertices). The zero entries in the matrix now imply all off points lie in only one of the two pairs of classes. This establishes the stated two-group reduction.

A separate exact check in the same script confirms that any four-point outside diamond whose rim tips and one shared-edge vertex have three contacts forces its fourth vertex either into M or into the same special class as that shared-edge vertex (72 off, 12 on across four modulo-translation anchors and unordered tip pairs). This is an optional refinement; no four- or five-point outside configuration is globally excluded by the ambient-group result.

**Reproducibility.** Run `node verify-moser.cjs` and `node verify-reduction.cjs`. These original scripts use only exact integer arithmetic; no third-party data or code are required. The derived partial upper bound is only the previously known 73 from deletion, conditional on the October reports. No priority claim or general asymptotic inference is made. The reflection z ↦ w conj(z) swaps the two index-three supergroups. Indeed, writing M = A + wA for A = Z[η], conjugation preserves A and this reflection swaps the summands. It sends (1+η)/3 to w(2−η)/3 = w−w(1+η)/3, so the two listed groups are isometric. The original M itself is a subgroup of either.
