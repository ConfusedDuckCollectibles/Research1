# Independent audit: off-module Moser unit circles

**Pass.** For any translate, rotation or reflection of the standard rank-four Moser module \(\mathcal M\), any point of the real plane outside that copy has at most three module neighbors at distance 1. Three is attained, and when attained the three neighbors form an equilateral triangle of side \(\sqrt3\). For the standard copy the points attaining three are exactly four additive cosets of \(\mathcal M\). This is not a result for the larger Moser ring or the full number-field plane. The consequence for a 25-point 73-distance configuration is a conditional split, not a global bound at 25.

## Independent complete geometric reduction

Write
\[
\eta=(1+i\sqrt3)/2,\qquad \omega=(5+i\sqrt{11})/6,\qquad
\mathcal M=\mathbb Z\langle 1,\eta,\omega,\eta\omega\rangle,
\]
and \(s=\sqrt{33}\). For integer coordinates \((a,b,c,d)\), write the physical point as \((X/36,\sqrt3Y/36)\), with
\[
X=36a+18b+30c+15d-3sd,\qquad Y=18b+15d+s(2c+d).
\]
The two rational components of each of \(X,Y\) determine \((a,b,c,d)\) uniquely: first use the coefficient of \(s\) in \(X\) for \(d\), then the coefficient of \(s\) in \(Y\) for \(c\), then the rational part of \(Y\) for \(b\), and finally the rational part of \(X\) for \(a\). This directly certifies a faithful coordinate system, with no assumption that abstract integer vectors have distinct images.

Suppose a real point has at least three unit neighbors in \(\mathcal M\). Translate one such neighbor to 0; this preserves \(\mathcal M\). Any other two neighbors \(A,B\) are noncollinear with 0 since a circle meets a line at most twice. In scaled coordinates the center is the unique solution of the two linear equations \(2(X_A X+3Y_A Y)=X_A^2+3Y_A^2\), and the same for \(B\); thus both center coordinates are in \(\mathbb Q(s)\). Every anchor then lies on the unit circle both in the physical embedding and after conjugation \(s\mapsto-s\), since its exact equation is over this field. The conjugate scaled metric still comes from the positive real quadratic form \(X^2+3Y^2\). The anchor 0 is fixed by conjugation. Consequently each difference vector from 0 must have length at most 2 in *both* embeddings. This use of conjugation does not require the arbitrary starting real center to belong to a prescribed field: the three anchors imply it.

To obtain a small complete coefficient bound, put \(u=a+b\eta\) and \(v=c+d\eta\). Conjugation in scaled coordinates keeps \(\eta\) and exchanges \(\omega\) with \(\bar\omega\). The mean of the two squared lengths is
\[
q=\frac{|u+\omega v|^2+|u+\bar\omega v|^2}{2}
=|u|^2+|v|^2+\frac53\Re(u\bar v)
\ge \frac{|u|^2+|v|^2}{6}
\ge \frac{a^2+b^2+c^2+d^2}{12}.
\]
Thus \(q\le4\) forces integer \(a,b,c,d\) within \([-6,6]\) and squared coefficient sum at most 48. We do not bound points by small physical length alone: \(\mathcal M\) is not a discrete planar lattice. The Euclidean isometries in the statement are handled simply by pulling back the configuration to this copy, no fundamental-region or symmetry-orbit completeness assumption needed.

## Independent exact enumeration

Let \(S(A)=X_A^2+3Y_A^2\). The script independently generates all 204 nonzero differences surviving both length bounds by using the coordinate map above and comparing signs of \(u+vs\) via squared integer magnitudes when necessary. For each unordered pair \(A,B\) with \(\Delta=X_A Y_B-X_B Y_A\ne0\), use exact Cramer coordinates (independently implemented before examining the author's replay):
\[
D=6\Delta,\qquad
N_X=3(S(A)Y_B-S(B)Y_A),\qquad N_Y=X_A S(B)-X_B S(A).
\]
Keep the center iff \(N_X^2+3N_Y^2=1296D^2\); normalize the two field quotients \(N_X/D,N_Y/D\) by rational conjugation and positive reduced denominators, then check every difference from the full list by the same perpendicular-bisector equality. Recover integer membership by the triangular coefficient inversion above. No calculations or signs use floating arithmetic. Since any circle with at least three module points supplies one enumerated noncollinear pair and all its normalized anchors pass the finite difference bounds, the resulting census is complete. There is no need to enumerate centers with exactly two anchors.

| Test | Exact count |
|---|---:|
| Nonzero vectors of length at most 2 in both embeddings | 204 |
| Noncollinear unordered candidate pairs | 20,196 |
| Pairs producing radius exactly 1 | 2,460 |
| Distinct normalized centers with at least 3 anchors | 30 |
| On-module centers, each with 18 anchors | 18 |
| Off-module centers, each with 3 anchors | 12 |

For each off-module center, all three pairwise squared lengths of its anchors are 3. After multiplying its center coefficient vector by 3 and reducing coordinatewise mod 3, the four residues are \(0011,0022,1100,2200\), each occurring three times. Every point with three module neighbors is a lattice translate of a normalized center, and every lattice translate retains those neighbors, so each residue indeed describes a full coset. An explicit sharp fixture is the center with coefficient vector \((-2,1,4/3,-2/3)\), with three neighbors \( (0,0,0,0),(-3,0,2,0),(-3,3,2,-2)\). Its nonintegral coefficients together with faithfulness prove it off-module.

## Exact implication for 25 points, and sources

Assume \(u(24)\le68\) and \(u(22)\le60\). In a 25-point set with at least 73 unit pairs, deletion of one vertex proves every vertex has degree at least 5. For a fixed isometric module copy meeting that set in at least 22 points but not containing it, write \(k\in\{1,2,3\}\) for the number outside. Any outside vertex would have degree at most \(3+k-1\), ruling out \(k=1,2\). For \(k=3\), there are at most \(9\) cross edges and \(3\) inside the outside set, so the remaining 22 points carry at least \(73-12=61\) edges, also impossible. Thus a putative 73-edge set either is wholly contained in some copy or meets each copy in at most 21 points. This does not exclude the all-module or genuinely sparse-in-every-copy alternatives. Upper bounds on a heuristic beam search do not replace an exact all-module branch proof.

- [Engel, Hammond-Lee, Su, Varga and Zsámboki, arXiv:2406.15317v3](https://arxiv.org/html/2406.15317v3#S2.SS1), Definition 2.2 and Theorem 2.5 give this module, its exact squared norm and 18 *on*-module unit directions. Their Table 2 yields the 72-edge lower construction at 25 and does not certify a 25-point upper bound. The off-module classification above is not stated there. [Alexeev, Mixon and Parshall, arXiv:2412.11914v2](https://arxiv.org/html/2412.11914v2#S1), Theorem 1, establishes through 21, lists older bounds \(u(22)\le61,u(24)\le72,u(25)\le78\), and explicitly notes an extremal 17-vertex graph outside the Moser ring; no assumption that all dense real configurations are on this module is justified.
- The primary computational notes [Artificial Wasteland, The Sixty-First Distance](https://artwaste.land/strata/the-sixty-first-distance/) and [The Sixty-Ninth Distance](https://artwaste.land/strata/the-sixty-ninth-distance/) report \(u(22)=60\), \(u(24)=68\) (the latter also depends on their \(u(23)=64\)). These are released, computer-assisted and described by their authors as not externally refereed. I checked the native statements and their reported trust boundary, not those large enumerations. Their updated [small-value table](https://artwaste.land/strata/every-edge-but-one/) records \(72\le u(25)\le73\); its upper bound also follows directly by summing deletions, \((25-2)e\le25\cdot68\). No full later 25-point result or same off-module three-neighbor classification surfaced in the targeted search. This is only a source-relative priority check.

The independent replay is `node independent-moser.cjs`; the exact expected transcript is `independent-moser-expected.json`. It does not import source datasets or the author's checker.