# Size 52 over \(\mathbb F_{13}^{3}\): no nine-uniform pencil

Step 4 of the size-52 proof. **If a spectral pair of size 52 exists, neither member has a nine-uniform pencil.** Both sides then have uniform pencil cap eight, at least one attains eight, and all mixed moments through degree seven still vanish. The remaining pencil-at-most-eight case is excluded in [step 5](../05-endgame/PROOF.md).

The proof uses the independently audited prior nonlayer reduction, recorded in [step 3](../03-nonlayer/PROOF.md). Its inputs here are: (i) plane cap 10, line cap 4 and every projective line meets each uniform-normal set in at most nine; (ii) nine in either member forces nine in both, and each such member is an injective graph over the same 52-point shadow, up to an affine change of coordinates; (iii) the both-nine branch has vanishing mixed moments through degree seven. The normalization, Radon atlas, and line-factor argument establishing these are not redone by the new verifier. The main mathematical antecedent for the moment and pencil method is [Ferguson–Iosevich §§8, 11](https://arxiv.org/html/2610.30712v1); the 10-pencil gate follows the earlier [Wang planar profile argument](https://arxiv.org/html/2609.20495v1) with its explicitly checked projection scope. No priority or publication claim.

## 1. Two intersecting ten-point sections

Let \(E=\{(x,y,h(x,y)):(x,y)\in B\}\), with the frozen shadow \(B\) below (encoded \(13x+y\)):
~~~
0 1 3 9 15 17 21 22 26 33 35 37 40 45 50 51 52 56 60 63 65 68 73 76 78 79
85 86 91 92 99 101 104 107 108 115 120 122 124 127 130 137 138 140 143 150
152 155 156 163 167 168
~~~
The affine base coordinate changes and vertical shears used here preserve the vanishing mixed moments (all monomials of degree at most seven span an affine-invariant polynomial space). The sections \(y=0\) and \(y=2x\) both contain ten points; their intersection in \(B\) is exactly \((0,0)\). The sections' respective \(x\)-supports are
\[
 T_0=(0,2,4,5,6,7,8,10,11,12),\qquad T_2=(0,1,3,4,7,8,9,10,11,12).
\]
Each section is a planar graph in coordinates \((x,h)\), of line cap four. It determines at most nine directions since its directions lie in a projective line of the partner's uniform-normal set.

For each fixed labelled support, all necessary height shapes are obtained by the following finite enumeration. At most nine directions are determined and vertical is missing, so also some finite slope \(m\) is missing. Subtract \(mx\) from the height; the resulting function on the fixed support is injective. Translate and scale its heights to make the first two values \(0,1\) in the indicated increasing \(x\)-ordering. Enumerate every other height injectively in \(\{2,\dots,12\}\), checking all new directions and line caps. The first implementation also searches every size-nine direction palette containing the initial slope and excluding slope zero; direction sets with fewer than nine extend to one such palette. The second counts directions directly and does not use outer palettes. Each has 25 distinct accepted normalized injective arrays for each section. After subtracting the affine height through the first two entries and making the first nonzero residual height one, both arrays give exactly six height shapes, one of eight directions and five of nine. These operations do not move the \(x\)-supports.

Now use a single global affine shear to set the heights at the origin, the second point of the first section, and the second point of the second section to zero. Those three base points are noncollinear. The two sections are nonaffine by cap four, so one can globally scale the height to make the first section's first nonzero height one. The second section's scale is an arbitrary nonzero \(t\in\mathbb F_{13}\). Thus exactly \(6\cdot6\cdot12=432\) labelled anchors must be tested, assigning the heights on the 19-point union. No independently modifying both section scales beyond this enumeration is assumed. It is enough to exclude these anchors; the remaining graph constraints can be dropped rather than exhaustively tested.

## 2. Linear and quadratic moments

Write \(x^a y^b\) as \(P_{ab}\). The already proved moment identities imply, over \(\mathbb F_{13}\),
\[
 \sum_B P_{ab} h=0 \quad(a+b\le6),\qquad
 \sum_B P_{ab} h^2=0 \quad(a+b\le5).
\]
There are 28 equations of the first type and 21 of the second type. The linear moment matrix has rank 28. With the 19 anchor assignments the coefficient matrix on all 52 heights has rank 46 (equivalently, on the 33 unpinned coordinates the rank is 27). It is consistent for all 432 anchors. For each anchor the possible linear solutions take the form
\[
 h=h_0+\sum_{i=1}^{6} t_i b_i ,
\]
with the same basis \(b_i\) for the unpinned homogeneous kernel and an anchor-dependent particular solution \(h_0\). All data are derived by exact modular elimination from the printed shadow, coordinates and shapes.

Quadratic expansion of the second type gives, for each of the 21 weights,
\[
 \sum_B P_{ab} h_0^2
 + 2\sum_i t_i\sum_B P_{ab} h_0 b_i
 + \sum_i t_i^2\sum_B P_{ab} b_i^2
 + 2\sum_{i<j}t_i t_j\sum_B P_{ab} b_i b_j=0 .
\]
Relax this by replacing the 21 degree-two monomials \(t_i t_j\) by 21 independent field variables. Any genuine height gives a solution of this linear relaxation, so inconsistency suffices. The rank of the anchor-independent block of 21 quadratic columns is 16. The paired rank results are:

| coefficient rank | augmented rank | anchors |
|---:|---:|---:|
| 16 | 17 | 12 |
| 18 | 19 | 120 |
| 19 | 20 | 300 |

Thus no relaxation is consistent. In the primary replay each inconsistent anchor has a direct certificate: a field vector \(\lambda\) of length 21 annihilates every coefficient column, while its dot product with the required right-hand side is 1. All 432 certificates are checked explicitly against their matrices; 48 distinct certificate vectors occur. The optional --cert argument prints the 432 normalized witnesses in anchor order as hexadecimal residue rows. The second replay independently enumerates height shapes by accumulating actual directions rather than slope palettes, eliminates only the unpinned coordinates and obtains lifted coefficients by finite differences of the square-moment function. It checks exact augmented ranks with a reversed-column elimination without using the annihilator certificates. It matches every rank total. No 9-point sections, geometric fit beyond the two anchors, time-truncated search, or solver tolerance enters the exclusion.

## 3. Consequence for the no-nine branch (additional analytic lemma)

The following argument supplied independently in this research pass keeps the degree-seven moments available even after excluding both-nine. Let \(B_*\) be a member of a hypothetical remaining spectral pair, and \(D=D(B_*)\) its projective directions. All tangent projective line sections have at most eight D-points because the partner has uniform pencil cap eight. Suppose a nonzero homogeneous degree-seven polynomial \(F\) vanishes on \(D\). A lower degree may be multiplied by a linear power first.

If every projective line has at least four D-points, set \(N=|D|\), \(n_L=|L\cap D|\). Nonfactor lines have \(n_L\le7\); at most seven line factors have \(n_L\le8\). Projective-plane counting gives
\[
 N^2-141N+5124=\sum_L(n_L-4)(n_L-7)\le7\cdot4=28.
\]
The left-hand minimum \(615/4\) contradicts this.

If a projective line has at most three D-points, its thirteen parallel affine slices each have occupancy four and are collinear. Indeed they average four, and any four noncollinear points would give at least four directions; a larger slice under line cap four contains such a four-subset. At least five of these distinct 4-lines have the same direction \(d\); their five images in the quotient by \(d\) are distinct and form a planar cap (three collinear quotient points would give twelve coplanar original points). These five quotient points determine at least five directions, since each of their ten pairs belongs to a slope with at most two pairs. Each such quotient direction arises from two parallel original 4-lines and gives at least seven cross-directions by Cauchy–Davenport, as well as \(d\). Therefore at least five distinct rational lines through \(d\), each with eight D-points, must be factors of \(F\). Fix one 4-line \(\ell\) of direction \(d\). Any selected point outside \(\ell\) gives four different non-\(d\) directions on the projective line associated to its affine plane with \(\ell\). If that projective line is not a factor, the five concurrent line factors already give multiplicity at least five at \(d\); adding the four new roots exceeds degree seven. Hence every selected point is on at most seven factor planes through \(\ell\), containing at most \(4+7(10-4)=46\) points, contradiction.

Thus no such \(F\) exists. The directional moments of the partner vanish on \(D\) by spectrality, so they vanish identically through degree seven. Polarization is valid in characteristic 13. Apply symmetrically. This analytic input is not needed for the both-K9 exclusion in Sections 1–2; it is a separate result narrowing the residual.

## Replay

Run from this directory:
~~~sh
node check-both-k9.js | diff - check-both-k9.expected
~~~
The default source uses no packages, only exact integer arithmetic with field reduction, and both replays ran together in V8. The [independent audit](AUDIT.md) includes a third, separate verification. The [endgame](../05-endgame/PROOF.md) uses the universal degree-seven rigidity above, with no K8 profile-atlas assumption.
