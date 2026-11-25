# Equality in the flag–volume inequality for simple polytopes

Let $P$ be a full-dimensional simple polytope. Then equality in [Winter's flag–volume inequality, Theorem 1](https://arxiv.org/html/2610.09990v1) holds **if and only if $P$ is affinely equivalent to a Cartesian product of simplices**. In particular these sharp examples are faces of Hanner polytopes. This is a classification in the simple class, with no symmetry assumption. It uses the face-difference equality criterion in the companion proof, not a result on Mahler minimizers.

**From equality to two-face geometry.** Suppose $P$ is sharp. By the criterion, for every face $F$ of $P$,
$$
(P-P)\cap \operatorname{span}(F-F)=F-F.
$$
For any edge $E=[v,w]$ in a two-face $F$, it follows that
$(F-F)\cap\operatorname{span}(w-v)=E-E$: the inclusion from right to left is automatic, and the reverse follows by restricting the edge identity for $P$. Thus the edge is a longest chord of $F$ parallel to itself.

To analyze the adjacent angles of an $n$-gonal face, take a side from $(0,0)$ to $(\ell,0)$ with the polygon above it. Let $0<\alpha,\beta<\pi$ be its interior angles at the left and right endpoints. For sufficiently small $t>0$, the slice of the polygon at height $t$ has $x$-endpoints $t\cot\alpha$ and $\ell-t\cot\beta$, since the boundary intersections are then on the neighboring sides. Its length cannot exceed $\ell$. Therefore
$\cot\alpha+\cot\beta=\sin(\alpha+\beta)/(\sin\alpha\sin\beta)\geq0$,
so $\alpha+\beta\leq\pi$. Summing this for all $n$ sides gives $2(n-2)\pi\leq n\pi$, hence $n\leq4$. For $n=4$, all neighboring angle sums equal $\pi$. The two consecutive turn angles along the oriented boundary then sum to $\pi$ at every pair, making opposite edges parallel. The quadrilateral is a parallelogram. Thus every two-face is a triangle or a parallelogram.

**From the combinatorial product to the affine product.** For simple polytopes of dimension $d\ge3$, [Wiemeler, Proposition 4.5](https://arxiv.org/abs/1401.0403), as stated for simple convex polytopes and proved directly in [Yu–Masuda, Theorem 2.1(b) and Appendix](https://arxiv.org/pdf/1609.05761), shows that triangular or quadrilateral two-faces imply the face lattice of a product of simplices. In dimension $2$, the preceding calculation already leaves a triangle or parallelogram, and dimension $1$ is a segment (in dimension $0$, the point may be viewed as a zero-dimensional simplex).

For a face-lattice isomorphism to $\Delta^{n_1}\times\cdots\times\Delta^{n_k}$, label the actual geometric vertices by
$$
 x(a_1,\ldots,a_k),\qquad a_j\in\{0,\ldots,n_j\}.
$$
An edge changes one coordinate, since the graph of each simplex is complete. A face obtained by changing two different coordinates between two fixed choices each is a two-face with four vertices, and so is a parallelogram. Consequently the vector of an edge changing coordinate $j$ from $0$ to $b$ does not change if any one of the other coordinates changes: use the opposite, consistently oriented edges of that two-face. Iterate over the other coordinates. Write $x_0=x(0,\ldots,0)$ and $w_{j,b}=x(0,\ldots,b,\ldots,0)-x_0$, with $w_{j,0}=0$. Changing the coordinates successively gives
$$
 x(a_1,\ldots,a_k)=x_0+\sum_{j=1}^{k} w_{j,a_j}.
$$
The nonzero $w_{j,b}$ are the $d=\sum_j n_j$ incident edge directions at $x_0$; they span the full tangent space and so are independent. The affine map from the product of the simplices with these factor vertices to $P$, given by the displayed sum, is therefore an isomorphism.

**Converse and Hanner scope.** On a simplex, an affine functional on a face extends preserving its minimum and maximum: assign values in that interval arbitrarily to the remaining vertices, then interpolate affinely on the entire simplex. For products every face is a product of faces and every affine functional is a sum of factorwise affine functionals. Their extrema over a product add; independently extend each summand within its range. The companion criterion therefore gives equality on every product of simplices. Finally, an $n$-simplex is affinely equivalent to a facet of the $(n+1)$-crosspolytope, so the product is a face of a product of crosspolytopes, which is Hanner. No implication that all sharp polytopes are simple is asserted.
