# Proof of the size-34 exclusion

**Theorem.** There are no spectral mates of cardinality $34$ in paired copies of $\mathbb F_{17}^3$.

Write $X\subset V$ and $Y\subset V^*$. For a set $E$, $D(E)$ is its set of projective secant directions, and $U(E)$ is the set of projective normals for which all seventeen affine levels contain two points. The transpose of a unitary character matrix is unitary. Orthogonality and the prime cyclotomic identity therefore give $D(X)\subset U(Y)$ and $D(Y)\subset U(X)$. Coordinates in opposite spaces are always paired dually.

The four steps below are in dependency order. The appendices contain the proofs of the intermediate lemmas; the finite gates and their exact inputs/outputs are indexed in [Reproduce](REPRODUCE.md).

## 1. Affine caps and the forced six-direction quotient

The odd-prime geometric reductions and the two-layer argument of [Ferguson–Iosevich, §§8–9](https://arxiv.org/html/2610.30712v1) provide the starting constraints; in particular each mate has at most two points on an affine line. The explicit $p=17$ profile gates successively exclude pencils with at least eight uniform normals and then pencils with seven. By the planar seven-arc direction bound, this also forces at most six points in every affine plane. The planar fact follows from [Beutelspacher–Wettl, Theorem 2.3(d)](https://doi.org/10.1007/BF01388411); an exact finite check is supplied with the seven-normal gate. See [Appendix A](A-AFFINE.md).

Both $U$ and $D$ have at most 86 directions. Indeed every projective line contains at least 12 points of the complement of $U$. [Ball, Theorem 4.2(ii)](https://web.mat.upc.edu/simeon.michael.ball/blocking.pdf) gives at least $17\cdot13=221$ points in the complement, out of the 307 of $\mathrm{PG}(2,17)$. The bound 86 rather than 91 is used later.

Call a normal **rich** if some affine plane with this normal contains five or six points. If both plane caps were four, perpendicular secant pairs define a bipartite graph of degree at most two on the two direction sets. The local weights come from two-point intersections in the 17 parallel planes. Its possible nontrivial cycles force alternating matchings whose polygons are coplanar, hence of length four, incompatible with $34\not\equiv0\pmod4$; the matching weights in the remaining components do not reach $\binom{34}{2}$. Details are in [the rich-plane lemma](A-AFFINE.md#a3). Thus some mate has a rich plane. Five or six points in one plane and at most six determined directions force the unique six-direction hexagon palette, via its five- and six-point planar gate.

A rich normal is nonuniform for its own mate, so it cannot be a secant direction of the other. Projection of that other mate along the rich normal is therefore injective, and its image is a Boolean set of 34 sites with precisely the prescribed six uniform slope families. No weighted or many-to-one projection is used in the rest of the proof.

## 2. Reciprocal richness and vanishing through degree five

In a forced Boolean chart the twelve nonuniform profiles have values in $\{0,\ldots,6\}$, mass 34 each, and total squared deviation from 2 equal to 510. Radon inversion and the Vandermonde kernel give the exact low-degree polynomial couplings for the profiles. The three complete atlas gates exclude zero rich families, at most two rich families, and exactly three rich families, respectively. These are gates on a single actual projective palette and its stabilizer; they do not relabel the twelve profiles independently. See [Appendix B](B-MOMENTS.md).

It follows that both mates are rich, and every rich normal is perpendicular to at least four opposite rich normals. The resulting incidence graph has no four-cycles. Each side of any component has at least 13 vertices. Equality would give an order-three subplane, impossible in a projective plane over the prime field $\mathbb F_{17}$, so there are at least 14 on each side.

The exact marked-quintic gate now shows that no nonzero homogeneous quintic can vanish on either whole difference-direction set. Its input is a real five-point affine section and the accompanying normalized exterior pair; after completing the linear kernels it checks exact pair compatibility in the remaining candidate pools. The same argument applies to both mates. Since any opposite degree-five moment form vanishes at the actual difference directions, the degree-five form is zero. Multiplying a lower homogeneous form by a suitable monomial gives the same directional obstruction at lower degrees. Polarization, with factorials below 17, consequently gives vanishing of every mixed affine moment through total degree five. The associated evaluation-space argument yields $r_X+r_Y\ge3$, where $r_E$ is the dimension of the affine degree-at-most-three vanishing space.

## 3. Reduction to two integral septic curves

The quadric argument eliminates global quadratic relations. Its companion moment-isolation lemma proves more: a quadratic vanishing at at least 24 marks vanishes at all of them. It is this stronger lemma that is used to control factors of the cubic relations. See [Appendix C](C-CURVES.md).

If $r_E\ge2$, two geometrically coprime cubic relations cut a complete intersection in projective three-space. The reduced fixed curve has degree at least four. If $r_E\ge3$, accounting for residual components, including a doubled line, gives degree at most seven. The exact linkage analysis in characteristic 17 excludes degrees four and five; the remaining fixed curves in this range are an integral genus-three sextic carrying at least 33 marks, an integral genus-five septic carrying all 34, or an integral genus-two/three sextic carrying at least 30 marks together with a residual line. The finite quartic-plus-twisted-cubic check belongs to this classification. The genus-two point-count gate gives at most 32 rational smooth points on the normalization in its application.

A second argument treats small cubic relation spaces without assuming this degree-seven cap. Any nonzero sextic form vanishing on the actual difference directions has middle catalecticant rank at least five. The proof uses the restricted binary hexagon and the rich-normal incidence to rule out the small apolar supports. The sixth-moment pairing bounds the rank of the opposite sextic moment form by $2r_E$. Thus $r_E\le2$ makes the opposite sixth moments zero. The resulting zero-moment curve classification and translate-cone lemma kill the sixth moments on the first side as well; its cubic evaluation space would then be isotropic of dimension at least 18 in a nondegenerate space of dimension 34, impossible.

When both $r$'s are at least three, apply the translate-cone lemma to the actual septic or the 30+ marks on the integral sextic component. Its filtration calculation rules out a sextic form on the secants: a product of a quintic difference and a linear difference vanishes to order two along the geometric diagonal, not order six. This gives sixth-moment vanishing on both mates. The residual genus-three sextic has only birational secant projections: an off-curve double projection would induce a reflection with six transverse fixed points, contrary to Riemann–Hurwitz. The consequent weight cap seven contradicts the necessary two-direction pencil weight 17. Thus the only residual is two integral, arithmetically Cohen–Macaulay, degree-seven genus-five curves, each containing all of its 34 marks and with moments through six zero.

## 4. Exact horizon identities and the final obstruction

Let $H_X,H_Y$ be the two **actual** length-seven infinity schemes, not generic hyperplane sections. Linkage and residue duality put every selected mark in the smooth locus and yield a distinguished smooth rational point $P_*$ on each curve, not assumed to be at infinity. Its degree-five relation has zero divisor consisting of the selected marks and $P_*$. The degree-seven sum functional on the coordinate algebra of each horizon is a perfect multiplication pairing. Quadratic restrictions inject, cubic restrictions surject, and the degree-five section $g$ annihilates the quadratic restrictions under that pairing. Isolating bidegrees $(2,5)$ and $(3,4)$ gives
$$
\ell^5=g_X a_Y+a_X g_Y,\qquad \lambda(K_3K_4)=0,\qquad \operatorname{rank}K_4\le3 ,
$$
where $\ell$ is the ordinary pairing and $K_j$ is a one-side image of $\ell^j$. These are identities of the actual finite algebras, including their local lengths. A line meeting a horizon in length at least four would force a nonzero quadratic relation or violate the one-dimensional upper bound for $\operatorname{Ann}(g)$. See [Appendix D](D-HORIZONS.md).

The mixed reduced/nonreduced cases are impossible by the local jet and fifth-power argument. For two nonreduced horizons, the actual septics give secant weights at most ten and $|D|\ge59$. The directional root bound and the long-row lemma eliminate local length four or more and the collinear three-support branches. The high-support pair lemmas reduce the remaining possibilities to $(3,3,1)$ and $(3,2,2)$ on at least one side. An active fat point has a globally unit value row; the quartic rank then forbids active triples and two active points. Meanwhile two inactive points would have scalar-zero groups carrying mass at least four each, whose forced overlap contradicts the separated fifth identity. Thus neither partition survives.

It remains to treat geometrically reduced horizons. If one $g$ is zero-free, the conic
$$
x\longmapsto\lambda_Y\!\left(a_Y\ell(x,\cdot)^2\right)
$$
vanishes identically, showing $a_Y$ proportional to $g_Y$. Applying the reverse argument either contradicts spanning immediately or reduces to $\ell^5=Cg_Xg_Y$ with both $g$'s nonzero. Short Frobenius orbits and the rigidity of affine maps on triples of fifth roots then place a reduced horizon on a conic, forbidden by quadratic injectivity. Consequently in the only possible reduced case both $P_*$'s lie in their horizons.

The last step has no orbit classification. The two distinguished representatives pair to zero. The separated fifth identity lets one scale eligible points to have ordinary pairing 1 with the opposite mark. In dual affine coordinates the six eligible sites on each side have forms $X=(1,c_i,a_i)$, $Y=(b_j,d_j,1)$. No more than two sites can share a value of $c$ or $d$, since the corresponding line includes its own mark. The functional $L_Y$ obtained from $\lambda_Y$ by dropping its marked coefficient is nonzero and kills all affine quadratics. The cubic
$$
T(C)=L_Y((b+Cd)^3)
$$
is nonzero because affine cubic evaluations span the functions on the eligible sites. Cubic–quartic orthogonality against the quartic row of the opposite mark gives $T(c_i)=0$. Therefore there are exactly three $c$ values. By symmetry there are exactly three $d$ values. The total weights on their three $d$-groups vanish by the Vandermonde system for $1,d,d^2$, so $L_Y(d^3)=0$ as well. This is the coefficient of $C^3$ in $T$; a nonzero quadratic cannot have three roots. This exhausts the last case and proves the theorem. $\square$

### Source context

The starting odd-prime strategy is from [Ferguson–Iosevich, *Fuglede’s conjecture holds in ...*, §§8–9](https://arxiv.org/html/2610.30712v1), and the broader prime-dimensional Fourier setting is described in [Malikiosis, *A linear programming approach to Fuglede’s conjecture*](https://link.springer.com/article/10.1007/s43670-023-00078-2). The affine pencil uses [Ball](https://web.mat.upc.edu/simeon.michael.ball/blocking.pdf) and [Beutelspacher–Wettl](https://doi.org/10.1007/BF01388411) as above. The characteristic-sensitive curve arguments use [Ballico–Miyazaki, §2](https://arxiv.org/pdf/math/9809042), [Peskine–Szpiro, *Liaison des variétés algébriques I*](https://doi.org/10.1007/BF01425554), and [Gruson–Lazarsfeld–Peskine, *On a theorem of Castelnuovo, and the equations defining space curves*](https://commack.math.stonybrook.edu/robert.lazarsfeld/Reprints/Gruson.Laz.Peskine.Castelnouvu.Bound.pdf). These links identify mathematical inputs, not independent review of this proof.
