# Excluding 104 and 117 from the spectrum cardinalities in F₁₃³

**Theorem.** No spectral subset of the additive group F₁₃³ has size 104 or 117.

We follow the general reductions and the moment method of Ferguson and Iosevich,
[arXiv:2610.30712v1, §§8 and 13](https://arxiv.org/html/2610.30712v1), which
prove the field-level result at p=11. The new finite input operates in the
affine plane of order 13. Two original replay implementations accompany this
note. The rest of the proof is elementary.

## 1. Spectral reductions and affine pencils

Let (E,A) be a spectral pair of common size 13m, for m=8 or m=9. In other
words, their character matrix with entries ζ^(e·a), ζ a primitive thirteenth
root of unity, has mutually orthogonal columns. The matrix is square, so its
rows are also mutually orthogonal. For either member B write U(B) for the
projective normal directions v for which the thirteen planes v·x=t each have
exactly m points of B, and D(B) for its projective difference directions.
Cyclotomic irreducibility gives D(A) ⊆ U(E) and D(E) ⊆ U(A): the polynomial
whose coefficients are the thirteen slice counts, if zero at ζ, is a multiple
of 1+X+...+X^12.

If an affine line l contains a points of B, the fourteen planes through it
have total occupancy 13m+13a, since each point not on l appears once and each
point on it appears fourteen times. Their normals are the projective line
given by the orthogonal complement of the direction of l. Therefore no
projective line can be contained in U(B): that would give 13m+13a=14m,
contrary to 0<m<13. Any affine plane with at least fourteen points would
determine every one of its fourteen tangent directions, by partition into
thirteen parallel lines in any one direction. This would give a forbidden
complete projective line in the partner's U. Thus every affine plane has
at most thirteen points.

For any affine line l, at least one plane containing it is uniform for B:
partition the partner's >13 points by the thirteen cosets of the orthogonal
complement of the direction of l; a pair in the same coset supplies the normal.
Combining the pencil count with one uniform plane and the thirteen-point
plane bound gives 13m+13a ≤ m+13², hence a ≤ 13−m. More generally, if
every projective line contains at least k points of U(B), then
13m+13a ≤ km+(14−k)13.

## 2. Three exact planar gates

The following are statements about subsets of F₁₃². A line cap c means that
no affine line contains more than c points.

| Size | Line cap | Lower bound on projective difference directions |
|---:|---:|---:|
| 8 | 5 | 6 |
| 8 | 3 | 8 |
| 9 | 4 | 8 |

The two verifiers cover the same gates via different normalizations.
In the Python verifier, a noncollinear triple exists since size exceeds the
cap. Affine-normalize it to (0,0),(0,1),(1,0), which forces directions
0,12,∞. If at most k directions occur, extend their set to a k-element mask
containing all three. For each such mask the algorithm traces all unordered
extensions of the triangle, requiring pairwise slope compatibility. The cap
test at insertion counts how many existing points share each slope from the
new point; a cap violation on any line must appear at the moment its final
point is inserted. The branch cutoff merely tests whether enough candidate
points remain to reach the target. No target set is reached for k=5,7,7.

In the C++ verifier, any hypothetical set in these gates has at least two
missing directions. Send them to the axes; distinct set points now have
distinct x-coordinates and distinct y-coordinates. Translate and scale the
coordinates separately so that two points are (0,0),(1,1). Its directions
then lie in the nonzero finite slopes and include 1. Extend them to a k-element
mask and enumerate all choices of occupied x's and their y's, with the
corresponding slope and line-cap tests. Again no target set is reached. The
different counts are normalized search states, not numbers of geometric
orbits.

## 3. The eight-point uniform-normal threshold

In every plane orientation, some parallel slice of A contains at least m
points, by averaging over thirteen planes. A subset of m of these gives
difference directions inside the corresponding projective line, all belonging
to U(E).

For m=9, the initial line cap is four. The third planar gate immediately
gives at least eight points of U(E) on every projective line.

For m=8, the initial line cap is five, and the first planar gate gives at least
six points of U(E) on every projective line. The same holds for U(A). For any
affine line in either member the refined pencil count now gives
104+13a ≤ 6·8+8·13 = 152, so a ≤ 3. Applying the second planar gate to the
partner's eight-point slices gives at least eight points of U(E) on every
projective line. Thus in both cases U(E), and symmetrically U(A), meets every
projective line in between eight and thirteen points.

## 4. Moments contradict both cardinalities

Put U=U(E), N=|U|, and n_L=|U∩L| for a projective line L. Suppose a nonzero
homogeneous degree-ten polynomial F over F₁₃ vanishes on U. There are 183
projective lines, fourteen through each point, and exactly one through each
pair. Thus Σ n_L=14N and Σ n_L²=N²+13N, and
$$\sum_L(n_L-8)(n_L-10)=N^2-239N+14640\ge 360$$
for integer N. On a projective line not dividing F there are at most ten
roots: its restriction is a nonzero degree-ten binary form. The summand
there is nonpositive. There are at most ten distinct projective line factors,
and each summand on one of them is at most (13−8)(13−10)=15. The sum is
at most 150, a contradiction. Multiplication by a power of a nonzero linear
form gives the same rigidity in every lower degree.

For 1≤j≤10, the moment polynomial M_j(v)=Σ_{e∈E}(v·e)^j vanishes for every
v∈U, by uniformity and Σ_{t∈F₁₃}t^j=0. Hence M_j vanishes identically. Also
the degree-zero moment is zero modulo thirteen because |E|=13m. For any
nonzero normal v, write n_t=|{e∈E:v·e=t}|. Its count vector, reduced modulo
thirteen, is orthogonal to t^j for 0≤j≤10. Their evaluations have rank eleven,
and the evaluations of 1 and t lie in their orthogonal complement and span
it. Consequently n_t=a_v t+b_v modulo thirteen.

A nonconstant affine residue profile uses residues 0 through 12 once each.
Because 0≤n_t≤13, its integer mass is 78 or 91, the latter only if the lone
zero residue is represented by 13. Since our masses are 104 and 117, every
profile is constant modulo thirteen. Such an integer profile is either
uniform m or consists of exactly m thirteens and 13−m zeroes.

Take an affine secant with a≥2 points of E, arbitrarily choosing two distinct
points to define it. The initial line cap gives a≤13−m. All fourteen planes
through the secant are nonempty, so each has occupancy m or 13. If b are
uniform, their total gives
$$(13-m)b=13(14-m-a).$$
Since 13 is prime, b is a multiple of 13, and 0≤b≤14 leaves b=0 or 13.
These give a=14−m and a=1 respectively, both forbidden. This last profile
argument is the p=13 instance of Ferguson–Iosevich Proposition 13.1. ∎
