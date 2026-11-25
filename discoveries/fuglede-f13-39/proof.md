# Excluding cardinality 39 in finite Fuglede for F₁₃³

**Theorem.** There is no 39-point spectral subset of the additive group
\(\mathbb F_{13}^3\).

The decisive additional ingredient is a set of concurrent factors for
degree-six and degree-seven moment forms. The short finishing computation
uses a single prescribed 7-point planar fibre. This proof does not use
the rank-one low-pencil census, the 364-mask rank-two search, or any earlier
partially scanned branch.

## 1. Reductions to a planar base

Write \(p=13\) and let \(E,A\) be a spectral pair of size \(3p=39\), so the
square character matrix \((\zeta^{e\cdot a})\) has orthogonal rows and
columns, where \(\zeta\) is a primitive \(p\)-th root of unity.
Let \(D(E)\) be the set of projective directions of differences of \(E\),
and \(U(A)\) the projective normals in which \(A\) has exactly three
points in each of the thirteen parallel planes. The minimal polynomial
\(1+T+\cdots+T^{12}\) of \(\zeta\) gives \(D(E)\subseteq U(A)\) and
\(D(A)\subseteq U(E)\): an integer occupancy polynomial of degree at
most 12 vanishes at \(\zeta\) exactly when all coefficients are equal.

Every affine line of either member has at most three points. Indeed some
two points of the partner lie in the same one of the thirteen cosets of
the plane perpendicular to its direction, so its direction lies in a
plane with a uniform normal. An affine plane has at most thirteen
points. Otherwise it determines all fourteen directions in its direction
plane and gives the partner a complete line of uniform normals. Such a
complete line is impossible: the fourteen planes through any fixed
affine line containing \(a\) points have total occupancy \(39+13a\), not
\(14\cdot3\). These preliminary reductions are the general finite-field
arguments of [Ferguson–Iosevich §8](https://arxiv.org/html/2610.30712v1).

### Projection identity and high-normal bounds

For a projective line of normals, project along its one-dimensional
annihilator, writing \(f\) for the resulting integer multiplicities on
\(\mathbb F_{13}^2\). Choose one good direction to be the vertical family.
If \(K\) good directions are removed, write \(L=14-K\) and \(g_i(t)\) for the
remaining line sums, parametrized as \(y-s_i x=t\) at distinct finite
slopes. Summing all lines through a point in the projection gives
\[
13 f(x,y)=\sum_{i=1}^L g_i(y-s_i x)-3(L-1).
\]
Each \(g_i\) has mass 39. Modulo 13 interpolate its values by a polynomial
of degree at most 12. The vanishing of the homogeneous terms of degree
\(j\ge1\) gives \(\sum_i c_{ij}s_i^k=0\) for \(k=0,\ldots,j\), where
\(c_{ij}\) is its coefficient of \(t^j\). The Vandermonde matrix has rank
\(\min(j+1,L)\), so every \(g_i\) has degree at most \(L-2\) and the possible
coefficient vectors in degree \(j\le L-2\) are precisely
\(c_{ij}=w_i Q_j(s_i)\), with \(w_i=(\prod_{h\ne i}(s_i-s_h))^{-1}\) and
\(\deg Q_j\le L-j-2\). This last assertion follows from Lagrange
interpolation and dimension.

If \(K=11\), the other three profiles are affine modulo 13, and if
\(K=10\) the other four are at most quadratic. A nonconstant affine or
quadratic has least-residue mass at least 78 or 52 respectively
(exact replay), so their residues must be constant. With actual entries
between 0 and 13, a profile is then either the constant 3 or has three
13s and ten zeros. The residue of the projection identity forces exactly
one profile of the latter kind. It follows from the integer identity that
the projection is Boolean and supported on three parallel lines. The
member is therefore supported on three full 13-point parallel planes,
with normal in the original projective line.

Neither member can be so supported. A 12-point planar set of line cap
three has at least eleven directions (exact gate M12), so a 13-point
plane section supplies eleven uniform directions for its partner.
Apply the \(K=11\) criterion to the partner and back to the first member:
one obtains dual coordinates in which the cross projections are injective.
More explicitly, view the first layer normal as a linear form \(n_E\) on
the first member's space and the second layer normal as a vector \(n_A\)
in that space (pairing with the partner). The criteria give \(n_E(n_A)=0\).
Take the second basis vector to be \(n_A\), take the first independent
of it in \(\ker n_E\), and take the third with value one under \(n_E\).
The dual third basis vector equals \(n_E\). Thus the first
projection is along its second coordinate axis and the partner's along
its third. The two layer forms are
\((x,f_s(x),s)\) and \((u,t,g_t(u))\), where each variable \(x,u\)
ranges over the entire field and \(s,t\) each run over three layer labels.
Their \(13\times13\) character blocks differ from the standard Fourier
matrix by diagonal unitaries. For two distinct block rows, orthogonality
would give three unitary \(13\times13\) matrices over \(\mathbb Q(\zeta)\)
whose sum is zero (divide each Fourier-product block by 13).
Left multiplication by the inverse of the first makes the equation
\(I+X+Y=0\). Unitarity implies \(X^2+X+I=0\).
But \(\mathbb Q(\zeta)\) does not contain a primitive third root, since
\([\mathbb Q(\zeta_{39}):\mathbb Q(\zeta_{13})]=2\), and an operator satisfying this irreducible quadratic would make the
13-dimensional vector space a module over a degree-two field extension,
forcing its dimension to be even. This is the same block obstruction as in
[Ferguson–Iosevich Lemma 10.2](https://arxiv.org/html/2610.30712v1).

If \(K=9\), five profiles remain, of degree at most three. If one is
cubic, all leading coefficients are proportional to \(w_i\), and the
quadratic coefficients have the form \(w_i(a+b s_i)\). A joint
translation in \(x,y\) centers the five polynomial arguments so all
quadratic coefficients vanish. The exact degree-three profile calculation
finds minimum least-residue mass 39 and precisely eight centered
polynomials of that mass: all \(a t^3+c\), with \(c=3\) or 9.
Nothing can be added to a zero residue without exceeding the mass 39.
Their five constants sum to \(15+6k\) with \(0\le k\le5\), never
\(12\bmod13\), as the identity requires. If no profile is cubic, the
same constant-profile argument forces the already impossible layering.
Thus \(U(E)\) and \(U(A)\) meet each projective line in at most eight
directions. A 10-point plane slice of line cap three gives at least
ten directions (M10), so each member has plane cap nine.

### Both members are graphs over one specified set

The exact profile census for \(K=8,L=6\), plane cap nine, yields 60 valid
normalized tuples, all with Boolean \(f\) and attained peak profile eight.
They lie in the same affine orbit. A canonical representative \(B\)
is encoded by \(13x+y\) as follows:

    0 1 2 3 10 11 12 13 14 16 17 20 28 30 31 56 58 59 73 79
    86 88 89 93 94 96 103 109 123 124 126 151 152 154 156 162 165 166 168

The exact census mechanism is recorded in the replay notes below.
If one uniform set has an eight-point line, its member has an actual
8-point affine slice from the integer profile peak. Every 8-point slice
of line cap three determines at least eight directions (M8), so the
partner also has an eight-point line of uniforms.

Suppose neither has such a line. Gate M8 gives plane cap seven. The
\(K=7,L=7\) profile census at this cap has zero survivors; thus all
projective lines meet either difference-direction set in at most six
points. For \(D=D(E)\), write \(N=|D|\), \(n_P=|D\cap P|\) for each
projective line \(P\). There are 183 such lines, fourteen through every
point, and one through every pair, so
\[
\sum_P(n_P-3)(n_P-6)=N^2-113N+3294\ge102.
\]
The left would be nonpositive if all \(n_P\in[3,6]\). Hence some tangent
plane orientation has at most two directions. Its thirteen affine
parallel planes each have at most three points, which must be collinear:
a noncollinear triple determines three directions, and the global
line cap is three. The total 39 implies thirteen collinear triples.
At least seven triples share one direction \(d\). Their projected
centres modulo \(d\) form a planar cap, for three collinear centres
would give nine points in an affine plane of cap seven.

For any pair of these \(d\)-parallel triples, the cross directions number
at least five by the three-element difference calculation; there can be
at most five because their projective line also contains \(d\).
Its equality case forces the triples to be arithmetic progressions
with common step (up to sign). Choose midpoint centres and scale the
step to 1. For chords of parallel projected centres with lifted vectors
\(w\) and \(\lambda w+\mu d\), the five cross directions must coincide.
Thus the affine map taking \([-2,2]\) to
\((\mu+[-2,2])/\lambda\) fixes that quintet. Its only affine stabilizers
are identity and negation; hence \(\lambda=\pm1,\mu=0\). Both projected
and lifted parallel centre chords equal up to sign. The exact seven
triple replay finds 210 normalized planar cap skeletons and 2730 eligible
lifted centre assignments after affine shear; each has a projective line
with 12 directions, contradicting the bound six. Therefore the graph
representation over \(B\) holds for both members.

## 2. Concurrent factors give all moments through seven

Work in an affine frame of \(E\) as a graph \(\{(x,y,h(x,y)):(x,y)\in B\}\)
and let \(u=[0:0:1]\) be the projective vertical direction. As the
projection is injective, \(u\notin D(E)\). Each base direction \(q\), including
infinity, gives a distinct projective tangent line \(\ell_q\) through
\(u\). The directions in \(D(E)\cap\ell_q\) are exactly the directions
of lifted pairs whose planar differences have slope \(q\). In particular
all \(\ell_q\) have at least one such point, since every planar direction
occurs in a subset of \(\mathbb F_{13}^2\) of size greater than 13.

Direct counting on \(B\) gives largest affine line-slice sizes for
slopes \(q=0,\ldots,13\) (13 means infinity)
\[
(3,7,8,3,7,3,3,3,7,3,7,3,3,7).
\]
Thus \(\ell_2\) has at least eight points of \(D(E)\) by M8. Each
\(\ell_q\) for \(q=1,4,8,10,13\) has at least six by the seven-point gate M7.
We also need \(\ell_8\) to have **at least seven**. For the affine line
\(y-8x=1\) its \(x\)-coordinates are \(0,2,4,6,7,11,12\).
A lifted slice is a seven-point height graph over exactly those abscissae
with at most three collinear. Subtract the unique affine function agreeing
at the first two abscissae, preserving both the line cap and the number
of slopes. Exhausting the five remaining field heights, pruning only
once seven slopes or a fourth collinear point appears, leaves zero
height graphs with at most six slopes. This is a self-contained exact
check and does not use scaling or orbit assumptions.

Consider any nonzero homogeneous degree-6 form vanishing on \(D(E)\).
The eight points on \(\ell_2\) force \(\ell_2\) to divide the form.
Restricting to any other line through \(u\), this factor gives a zero
at \(u\). For each of the remaining five large-slice lines this adds
one to at least six distinct zeros, so each is a line factor too.
Now a line through \(u\) outside these six has multiplicity at least
six at \(u\), plus a point in \(D(E)\) different from \(u\). Thus it also
divides the form, impossible in degree six. Since multiplication by a
nonzero linear power raises any lower-degree homogeneous form to
degree six without removing zeros, all lower-degree forms are excluded.

For degree seven, the first factor is again \(\ell_2\). The seven
points on \(\ell_8\), plus \(u\), force \(\ell_8\) as a distinct factor.
The two factors give multiplicity at least two at \(u\) on every other
line, so all remaining four large-slice lines (each with at least six
points) are factors. Once six distinct factors through \(u\) are known,
a line through \(u\) with at least two points of \(D(E)\) must itself be
a factor. Among the other eight lines at most one can have two such
points. At least two of \(\ell_0,\ell_3,\ell_5\) consequently have exactly
one such point each. For each of these two slopes, the lifted differences
on all base lines of that slope have one fixed height increment per
base direction vector. Subtract a linear function on the base taking
those prescribed increments on the two independent selected direction
vectors. The remaining height is constant on every edge of the graph
on \(B\) joining two points if their base direction is one of the
selected two slopes. This graph is connected for each pair from
\(\{0,3,5\}\), by direct finite graph calculation. Thus the height
would be affine and the whole 39-point member would lie in one affine
plane, a contradiction.

The use of multiplicity is elementary: restricting a nonzero homogeneous
form to any projective line which is not one of its factors yields a
nonzero binary form of the same degree; all distinct already known
concurrent linear factors pull back to nonzero multiples of the same
linear form whose root is \(u\). Hence their contributions in the
restriction add as stated and do not count any point of \(D(E)\) twice.

In the corresponding **dual** frame of \(A\), every form
\(M_j(v)=\sum_{a\in A}(v\cdot a)^j\) for \(1\le j\le7\) vanishes on \(D(E)\):
the scalar products there are uniform, and the sums of these powers over
the field are zero. The rigidity just proved forces each \(M_j\) to be
identically zero. Its multinomial coefficients are all invertible
because \(j<13\); thus every mixed moment through degree seven vanishes.
For degree zero use \(39=0\) in the field. This property is invariant
under affine changes because it is the vanishing of the sum of every
polynomial of total degree at most seven. We may now pass to \(A\)'s
own graph frame over \(B\).

## 3. Symmetry and the determinant

The base \(B=-B\) consists of the origin and nineteen antipodal pairs.
Subtract a constant from the height to make \(h(0)=0\), and subtract
a linear function to make its odd part zero at representatives
\((0,1),(1,0)\). Write \(e,o\) for the even and odd parts on the nineteen
nonzero pairs. The linear-in-height mixed moments of degree at most
seven say that \(e\) is orthogonal to the even base monomials of degree at
most six, and \(o\) to the odd such monomials. The normalization adds
the two specified odd constraints. The dimensions are 3 and 5 exactly.
The quadratic-in-height moments with odd base monomials \(w\) of degrees
1,3,5 give \(\sum w e o=0\), because
\((e+o)^2-(e-o)^2=4eo\). At every one of the 2196 nonzero \(e\)'s in
its three-dimensional space, their coefficient matrix on the five
odd variables has rank 5. Thus \(o=0\) if \(e\ne0\). Then the
quadratic-in-height moments with even base monomials \(w\) of degrees
0,2,4 require \(\sum w e^2=0\), and no such nonzero \(e\) passes.
All these modular ranks and equations are checked directly by the
standalone exact finish script. It follows that the normalized graph
has odd height, so \(A\) is centrally symmetric up to translation.

For completeness the [Ferguson–Iosevich Corollary 10.11](https://arxiv.org/html/2610.30712v1)
determinant obstruction in this case has a short verification. Translate
the centre to zero and let \(\Delta\in \mathbb Q(\zeta)\) be the determinant
of the \(39\times39\) character matrix. Negation of the 39 vertices
has sign \((-1)^{19}=-1\), hence \(\overline\Delta=-\Delta\).
Orthogonality gives \(\Delta\overline\Delta=39^{39}\).
Therefore \((\Delta/39^{19})^2=-39\). But conjugation on
\(\mathbb Q(\zeta)\) is the square of the Galois automorphism
\(\zeta\mapsto\zeta^5\), since \(5^2\equiv-1\bmod13\); it fixes
every quadratic element of this field, as any automorphism acts on the
two conjugate roots by a permutation whose square is trivial. This is
impossible for \(\sqrt{-39}\). ∎

## Replay, priority and residual cases

The exact planar gates M12, M10, M8, and M7 respectively exclude
\((12,10,3),(10,9,3),(8,7,3),(7,5,3)\) (size, maximum directions,
line cap). Normalize two **missing** directions to the two axes, and
then two points to \((0,0),(1,1)\); every remaining \(x\) and \(y\)
must be used at most once, and every secant is tested against a
direction mask. These cover all hypothetical sets. The small-profile
and common-step-AP script checks every listed centred polynomial
and every 3-set pair. For the higher profile gates, the C++ script
enumerates all degree-\(\le L-2\) positive-mass nonconstant profiles
by their values at \(0,\ldots,L-2\); the coefficients uniquely determine
the polynomials. Because the profile cap is \(<13\), their integers
are exactly the least residues. Normalize a good direction at
infinity, enumerate the finite bad-slope masks up to affine
permutation, use the Vandermonde identities degree by degree, and
use the two independent translations of the arguments of the first two
profiles to center them. Test the exact integer projection identity
and directly replay all fourteen line sums. At \((L,\mathrm{cap})=(6,9)\),
14 masks, 1040 profiles, 1032 completed positive-degree tuples and 60
survivors occur, all Boolean with peak 8, all in the orbit of \(B\)
under normalized antipodal centre and ordered affine frames. At
\((L,\mathrm{cap})=(7,7)\), 14 masks, 3328 profiles, 736 positive-degree
tuples and no survivors. The seven-triple program reports 210
normalized cap skeletons and 2730 exact lifted configurations,
with pencil-size histogram \(\{12:2730\}\). No solver timeout or
incomplete calculation is used.

The 2026 primary frontier used here is [Ferguson–Iosevich](https://arxiv.org/abs/2610.30712),
which treats \(\mathbb F_{11}^3\), including size 33 by their determinant
method and other uses of line-factor multiplicities in Section 11. The [Malikiosis upper-cardinality theorem](https://doi.org/10.1007/s43670-023-00078-2)
does not exclude 39. A bounded current primary search found no previous
size-39 exclusion in \(\mathbb F_{13}^3\), and no priority claim is made.
After the separate [size-26](../fuglede-f13-26/) and
[size-104 and size-117](../fuglede-f13-104-117/) exclusions and the primary general reductions, the unresolved cardinalities
in this \(\mathbb F_{13}^3\) ledger are **52, 65, 78 and 91**. This is
not full Fuglede for \(p=13\).
