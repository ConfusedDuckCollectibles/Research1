# Mixed alphabets and a quantitative gap

The [theta certificate](proof.md) also classifies **all** maximum intersecting families when alphabet sizes \(s_i\ge2\) vary. On positive cells, each maximum is an original largest star, or is given by an antipodal orientation on one class of binary coordinates. That class consists of binary coordinates with identical original positive-support incidences and maximum original point-star weight. For every support containing the class, choose all assignments whose restriction to it belongs to that orientation. An orientation selects exactly one of each pair of opposite binary assignments. Thus nonstar maxima exist precisely when such a class contains at least three binary coordinates, even when other coordinates have larger alphabets.

## Equality

Complete singleton weights exactly as in the [proof](proof.md), so the cap \(M\) is unchanged. A maximum family containing a singleton is already its full original star. Otherwise its weighted indicator \(y_U=\sqrt{w'_U}\,\mathbf1_{\mathcal A}(U)\) has the same mass under completion and lies in the kernel of \(P=H-f'f'^{\mathsf T}/M\), as in [the common-alphabet proof](equality.md).

The sector \(J=\varnothing\) has matrix \(I-\Pi\); singleton sectors have the \(z_i\) kernel on cell indices, tensored with the point-label contrasts. If \(|J|\ge2\) and \(J\) contains a nonbinary coordinate, its sector has eigenvalues at least \(1/2\). For even \(|J|\), use \(W_J\succeq-I\) and the fact that the divisor \(\prod_{i\in J}(s_i-1)\) is at least two. For odd \(|J|\) one can use a binary \(i\in J\) and \(\lambda_{\max}(W_J)\le1\); if all coordinates of \(J\) are nonbinary, the product of the other factors is at least four and the gap is at least \(3/4\).

After dividing out \(\sqrt{w'_U}\), the indicator on a cell is therefore a function of its binary labels plus a sum of first-order functions of its nonbinary labels. Each latter function is shared across cells containing its coordinate, since the normalized singleton kernel carries \(z_i\), and undoing the constant-cell normalizations yields the same function up to a fixed factor \(\sqrt{s_i}\) independent of the cell. If a nonbinary function is active, Boolean additivity implies all other summands are constant in every cell containing that coordinate: independent nonconstant summands would have integer-valued differences and total range at least two. The accepted label set is identical across these cells because the shared varying function must be shifted to have range \(\{0,1\}\). It is exactly one label, since two would permit disjoint accepted faces in one cell. Those cells already provide a completed star of weight \(M\). Positive weights force the optimizer to be that star, hence also an original largest star.

Otherwise the indicator depends only on binary labels. A cell with no binary labels cannot be accepted, since it would contain disjoint selections. In any other cell at most one of two antipodal binary assignments is accepted: if both were, the nonbinary labels could be chosen different on every nonbinary coordinate to give disjoint faces. Thus the odd function \(g(x)=h(x)-h(-x)\) still encodes \(h=\mathbf1_{g=1}\). The [binary block argument](equality.md) applies without change to the binary coordinates: odd-sector equality forces identical completed positive-support incidences; the half-integer case cannot attain \(M\), using the containing-cell mass \(2M\) for every binary coordinate; and the single integer class must choose a full antipodal orientation. For a nonstar class, a positive added singleton at any center would split its incidence from the other members, so the centers are originally maximal and their incidences coincide in the original family. Conversely a binary class of size at least three supports a nonstar optimum by majority on any three of its coordinates. Other labels have no effect. This proves the statement.

## Energy estimates for near-maxima

Let \(\mathcal A\) be intersecting and contain no singleton, with weight \(a\le M\). Its completed weighted indicator \(y\) then satisfies \(\lVert y\rVert^2=a\) and the exact identity
\[
y^{\mathsf T}P y=a-a^2/M=a(M-a)/M=:E.
\]
Indeed \(H\) is identity on the intersecting support of \(y\). With all alphabets at least three, every sector \(|J|\ge2\) has gap at least \(3/4\) (even sectors have divisor at least four, odd sectors have at least two other factors at least two). In the constant and singleton sectors the nonkernel eigenvalues are at least one. The kernel is exactly the weighted image of the global additive incidence space, that is, functions that sum independent single-coordinate functions shared across cells. Hence, in the completed cell-weighted squared \(\ell^2\) norm, the distance of the indicator to this space is at most \(4E/3\) **when squared**. This does not by itself bound distance to a star.

For mixed alphabets the sum of squared weighted norms of the components in higher-order sectors involving a nonbinary coordinate is at most \(2E\). For an odd all-binary sector \(J\), put
\[
\rho_J=\frac{\sum_{S\supseteq J}N_Sw'_S}{2M}.
\]
On the supports containing \(i\in J\), \(W_i\) has eigenvalue one on \(z_i\) and nonpositive quadratic form on its orthogonal complement. Thus \(W_i\preceq z_i z_i^{\mathsf T}/\|z_i\|^2\); its restriction to the supports containing \(J\) gives \(W_J\preceq\) a rank-one matrix of norm \(\rho_J\). Thus the odd binary sector contributes at least \((1-\rho_J)\|y_J\|^2\) to \(E\). These contributions can be summed over all odd all-binary \(J\), since the sectors are orthogonal.

## Exact checks

Run the [standalone replay](mixed-replay.cjs) and compare with [expected output](mixed-expected.json). It enumerates maximal intersecting families across 254 weighted mixed-alphabet support systems, plus two four-coordinate examples with exactly ten nonstar maxima each. These finite checks illustrate the analytic classification; they do not establish the general theorem.

AI-assisted research; independent deductive and exact checks, not externally peer reviewed.
