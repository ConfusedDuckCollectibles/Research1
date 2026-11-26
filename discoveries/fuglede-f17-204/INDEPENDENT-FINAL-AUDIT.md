# Independent atlas and final contradiction at size 204

**Status:** Full nonexistence of 204-point spectral pairs in \(\mathbb F_{17}^3\) passes this independent audit. This uses the linked-layer reduction in INDEPENDENT-AUDIT.md, an independent finite atlas here, and then a rank contradiction. It does not address any other unresolved cardinality.

## Complete 12-graph atlas

For every 12-element field subset \(T\) and every set \(I\) of at most 11 finite slopes, define
\[
\mathcal H(T,I)=\{h:T\to\mathbb F_{17}\ :\ h(t_0)=0,\ D(\{(t,h(t)):t\in T\})\subseteq I,\quad\text{line cap }5\},
\]
for any fixed \(t_0\in T\). Then \(|\mathcal H(T,I)|\le4\), and \(\operatorname{affdim}\mathcal H(T,I)\le2\) for nonempty cells. It suffices to extend \(I\) arbitrarily to 11 slopes, since restriction to a subset only decreases cardinality and affine dimension. Vertical is automatically absent from every graph. Independently normalize \(T\) and the **complement** of the allowed 11 slopes: an affine \(x\)-coordinate change normalizes \(T\), which temporarily rescales slopes. Then \(y\mapsto a y+b x+c\) acts by the full affine group on slopes without further changing \(T\); the constant \(c\) restores the choice of normalized point. For fixed transformations, the induced map on the height vectors is a coordinate permutation, multiplication by \(a\ne0\), and addition of a fixed vector, with a possible change of zero-reference via subtracting one coordinate. That last map is invertible between the two normalized height spaces. Hence the cardinality and affine dimension are unchanged.

My independent verifier enumerates all 272 maps \(x\mapsto a x+b\), not someone else's representative list. This partitions all 6188 domain masks into 25 orbits, and all 12376 six-missing-slope masks into 49 orbits. For each of the 1225 cells it fixes \(h(t_0)=0\), uses 17-bit candidate-value domains for all other variables, and chooses the smallest remaining domain. After selecting a new height it removes/disallows each candidate whose pair slope is not allowed, and removes the candidate lying on any affine line that just reached five. Thus both filters are necessary, and all choices remaining at each step are explored. At a completed tuple it independently rechecks all 66 pair slopes; each line's pair count at most ten additionally verifies line cap five. It computes affine rank by exact elimination of the height-vector differences from the first accepted tuple. This is separate from the author's census (different branching and orbit generation).

Exact output: 558992 visited nodes, 207 nonempty cells, 353 normalized rows, maximum 4 rows and maximum affine dimension 2. Cell-cardinality histogram, including empties: \(1018,84,108,7,8\). Among nonempty cells the affine-rank histogram is \(84,108,15\) in ranks \(0,1,2\). The exact script and output are independent-atlas.cjs and independent-atlas-expected.json; remove fences and run with Node.js, no imports or data files. The exact V8 run took 8.2 seconds.

## How the atlas excludes the linked pair

Use the perpendicular layer coordinates
\[
E=(x,f_s(x),s),\quad A=(u,t,g_t(u)),\quad |S|=|T|=12
\]
from the independently verified preliminary reduction. There are at least fourteen active \(r\) for the \(E\)-layer Fourier coefficients, and for each active \(r\), all cross slopes \(I_r\) formed between the \(A\)-layers at \(u-u'=r(t-t')\) satisfy \(|I_r|\le11\), as proved in INDEPENDENT-AUDIT.md. For fixed such \(r\) and varying \(x\in\mathbb F_{17}\), the cross sections through \(u=x+rt\) give the normalized height rows
\[
H_r(x)_t=g_t(x+rt)-g_{t_0}(x+rt_0).
\]
Their pair slopes are in \(I_r\), and they inherit line cap five from a linear parametrization of the affine plane \(u-rt=x\) in the original 3-space. Thus all rows lie in one atlas cell, with at most four rows and affine rank at most two.

Take the reduced polynomial representatives of the \(g_t\), all of degree at most 16. For each active \(r\), interpolation in \(x\) shows that the nonconstant coefficient-row vectors of \(H_r(x)\) all belong to the two-dimensional-or-smaller linear direction space of its affine hull. To see this explicitly, subtract the constant vector \(H_r(0)\) from all rows and then apply the linear inverse Vandermonde map to their values, leaving every nonconstant coefficient unchanged. Restoring the baseline polynomial \(g_{t_0}(x+rt_0)\) to every column adds a matrix of rank at most one. The nonconstant coefficient matrix of the twelve polynomials \(g_t(x+rt)\) therefore has rank at most three at every active \(r\).

For any pair of layer indices the difference \(g_t(x+rt)-g_{t'}(x+rt')\) is determined by a linear function of \(H_r(x)\) and takes at most four distinct values. If both polynomials have degree at most four, their difference is necessarily constant for every active \(r\): any nonconstant polynomial of degree at most four has fiber size at most four and therefore image at least \(\lceil17/4\rceil=5\). Pick two distinct active \(r_1,r_2\). The constant identities express \(g_t(y)\) as \(g_{t'}(y+r_i(t'-t))\) plus a constant for each \(i\). Subtraction gives a constant finite difference of \(g_{t'}\) along the **nonzero** step \((r_1-r_2)(t'-t)\). Iterating a nonzero step traverses all 17 field values and shows \(g_{t'}\) and \(g_t\) both affine, contradicting the original line cap (each affine layer would be a 17-point line). Hence at least eleven polynomials have degree five or larger.

At any fixed active \(r\), shifted polynomials of different nonzero degrees have independent nonconstant coefficient columns by leading degree. Thus there are at most three distinct degrees among these eleven and at least four share one degree \(d\ge5\). Let their distinct labels be \(t_1,\dots,t_4\) and their nonzero leading coefficients \(a_i\). Use the \(4\times4\) minor of their shifted nonconstant coefficient matrix on coefficient rows \(d,d-1,d-2,d-3\). Its row \(j=0,1,2,3\) entry in column \(i\) has degree at most \(j\) in \(r\), with leading term \(a_i\binom dj t_i^j r^j\). The determinant is of degree at most \(6\) and its coefficient at degree 6 is
\[
\left(\prod_i a_i\right)\left(\prod_{j=0}^{3}\binom dj\right)\det(t_i^j)_{0\le j\le3,1\le i\le4}\ne0.
\]
The binomial factors are nonzero as \(5\le d\le16<17\), and the \(t_i\) are distinct. But the determinant vanishes at every active \(r\) by rank at most three, at least fourteen roots. This is impossible. The full size-204 exclusion follows without the optional 12-point \(\ge10\)-directions lemma; the atlas supplies precisely what is needed directly.

## Positive controls

The independently generated 17-point witness
\[
(0,0),(2,1),(1,3),(3,9),(9,13),(8,4),(10,8),(14,6),(4,7),
(16,11),(5,10),(15,12),(11,5),(12,2),(13,14),(7,15),(6,16)
\]
has affine line cap five and exactly 14 directions (missing \(0,1,4,\infty\)), by direct checks of all 18 line-family profiles. The separate 12-point witness in the preliminary audit has cap five and exactly ten directions. These confirm sharpness of the two planar lower bounds, not possible spectrality.

## Prior scope

The sources establish the conditional layer cutoff (Malikiosis), the low-special-direction bound (Ghidelli), and classification of exactly four special directions for Boolean sets (Wang), not this field/cardinality contradiction. The affine atlas and rank bound supply the missing linkage. No externally peer-reviewed size-204 conclusion was located in the examined prior sources; that is a literature-check limit, not a definitive priority claim or a classification of all \(\mathbb F_{17}^3\) Fuglede cases.