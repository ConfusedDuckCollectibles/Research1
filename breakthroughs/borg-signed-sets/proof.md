# Exact weighted theta for complete signed cells

## 1. Statement

Let \(\mathcal S\) be any finite family of nonempty subsets of \([r]\). Coordinate \(i\) has an alphabet of size \(s_i\ge 2\). A face of support \(S\in\mathcal S\) is an assignment of one label at each coordinate in \(S\), viewed as a set of labelled points. Include every face in the cell of \(S\), each with weight \(w_S\ge 0\). Two faces intersect if they share a coordinate and have the same label there. Omit zero-weight vertices. Let \(M\) be the maximum weight of the star consisting of all faces containing a given labelled point. Then the weighted Lovász theta of the disjointness graph (edges join distinct disjoint faces) is exactly \(M\). If there are no positive-weight faces both sides are zero.

Use the convention
\[
\vartheta(G;w)=\max\{\langle ff^{\mathsf T},Y\rangle:Y\succeq0,\ \mathrm{tr}\,Y=1,\ Y_{UV}=0 \text{ on edges}\},
\qquad f_U=\sqrt{w_U}.
\]
The stars have support on independent sets, giving \(\vartheta\ge M\). In fact \(Y=(f\mathbf1_{\rm star})(f\mathbf1_{\rm star})^{\mathsf T}/M\) for a largest star attains it. For the upper bound it suffices to find a symmetric \(H\) with diagonal 1, zero entries on distinct intersecting faces, and \(H-ff^{\mathsf T}/M\succeq0\). Such an \(H\) satisfies \(\langle H,Y\rangle=1\) on primal feasibility.

## 2. Balance point stars, not cell sizes

The weight of a star at any one label of coordinate \(i\) depends only on \(i\); call it \(M_i\). Add \(M-M_i\) to the weight of each singleton face of coordinate \(i\), inserting the complete singleton cell if necessary. Each added singleton intersects exactly its own labelled-point star, so after completion all point stars weigh \(M\). In Sections 2–4 use the completed weights \(w'_S\) and keep only positive-weight cells. Write
\[
N_S=\prod_{i\in S}s_i,\qquad z_S=\sqrt{N_Sw'_S},\qquad V_{Si}=z_S\mathbf1_{i\in S},\qquad
W=zz^{\mathsf T}/M-\Pi,
\]
where \(\Pi\) is the orthogonal projection onto the column space of \(V\), with no independence assumption. For each coordinate \(i\), using \(z^{\mathsf T}V_i=s_iM\), we have
\[
W\succeq-I,\qquad W V_i=s_i z-V_i. \tag{1}
\]
If \(W_i\) denotes the principal block on supports containing \(i\), and \(z_i\) the corresponding strictly positive restriction, then \(W_i z_i=(s_i-1)z_i\). On \(z_i^\perp\) the rank-one term vanishes and the principal block of \(-\Pi\) is negative semidefinite. Hence \(s_i-1>0\) is the largest eigenvalue of \(W_i\), and its eigenspace is precisely \(\mathbb R z_i\). This conclusion does not require entrywise nonnegative \(W\).

## 3. Lift to labelled faces

For two faces \(U,U'\) with supports \(S,T\), let \(A_{UU'}=0\) if they intersect. If disjoint, set
\[
A_{UU'}=\frac{W_{ST}}{\sqrt{N_SN_T}\prod_{i\in S\cap T}(1-1/s_i)} .
\]
Then \(H=I+A\) has the required diagonal and intersecting entries. Fix a labelled point \(a=(i,\ell)\) and let \(v_a(U)=\sqrt{w'_S}\mathbf1_{a\in U}\), \(f'_U=\sqrt{w'_S}\). We show \(Hv_a=f'\). If \(a\in U\), the adjacency term vanishes and the identity term is \(f'_U\). Otherwise, in a target cell \(T\ni i\), the fraction of all faces in that cell that are disjoint and have the prescribed label at \(i\) is \(\prod_{j\in S\cap T}(1-1/s_j)/(s_i-\mathbf1_{i\in S})\). Summing against \(v_a\) gives
\[
(Av_a)(U)=\frac{(WV_i)_S}{\sqrt{N_S}(s_i-\mathbf1_{i\in S})}=\frac{z_S}{\sqrt{N_S}}=f'_U
\]
by (1).

## 4. Positivity in every tensor sector

The Euclidean space of each cell \(S\) is \(\bigotimes_{i\in S}\mathbb R^{s_i}\). Split each \(\mathbb R^{s_i}\) into the normalized constant direction and its orthogonal mean-zero space. In the block from cell \(T\) to cell \(S\), the normalized transfer for a common coordinate \(i\) is \((J_{s_i}-I)/(s_i-1)\), with \(J_{s_i}\) the all-ones matrix. Its eigenvalues on the two factors are 1 and \(-1/(s_i-1)\). For a coordinate belonging to only one support the transfer is the normalized constant vector or its transpose. Indeed the product of those transfers with multiplier \(W_{ST}\) gives exactly the displayed formula for \(A\) (zero unless every shared coordinate has distinct labels).

Thus the sector indexed by a set \(J\) of mean-zero coordinates is preserved and consists of cells \(S\supseteq J\). On its support index space, tensored with the identity on the common mean-zero label factors, \(A\) acts by
\[
(-1)^{|J|}W_J\big/\prod_{i\in J}(s_i-1).
\]
Here \(W_J\) is a principal block of \(W\) on supports containing \(J\). If \(|J|\) is even, including \(0\), the bound \(W\succeq-I\) suffices for \(A\succeq-I\) on this sector. If odd, choose any \(i\in J\); \(W_J\) is a principal block of \(W_i\), so its largest eigenvalue is at most \(s_i-1\le\prod_{j\in J}(s_j-1)\). This also bounds \(A\) below by \(-I\) on that sector. All sectors together give \(H\succeq0\).

As \(Hv_a=f'\) and \(v_a^{\mathsf T}Hv_a=M\), Gram Cauchy–Schwarz gives \(H-f'f'^{\mathsf T}/M\succeq0\). On original positive-weight faces restrict this matrix and put \(D_{UU}=\sqrt{w_U/w'_U}\). The scaled witness
\[
H_0=D H D + I-D^2
\]
has diagonal 1, zero at distinct intersecting pairs, and satisfies \(H_0-ff^{\mathsf T}/M\succeq0\), since the latter is the congruence of the completed certificate plus the nonnegative diagonal \(I-D^2\). This gives the asserted exact theta with the original weights.

## 5. Scope and relation to earlier work

The [classification of equality](equality.md) proves the strictness part of [Borg's original 2007 signed-set conjecture](https://doi.org/10.37236/959) and describes every maximum family when there are two labels per coordinate. The ordinary maximum-star bound is not a priority claim: [Chang–Liu–Liu, 2026, Proposition 5.3](https://arxiv.org/html/2609.19123v1#S5) proves the maximum-star bound for nonnegative decreasing weights on the Boolean cube. That applies through [Borg's 2011 Lemma 4.3](https://staff.um.edu.mt/peter.borg/publications/conjectures%20gen.pdf#page=8) to the ordinary signed-set statement, which his 2011 paper calls Conjecture 3.1. The same reduction with cell-constant weights and varying alphabet sizes counts compressed signed faces by the subset \(T\) of coordinates carrying a distinguished label. Its induced weight is \(\sum_{S\supseteq T} w_S\prod_{j\in S\setminus T}(s_j-1)\), a decreasing function because each factor \(s_j-1\ge1\). The new contribution here is a different, explicit sharp semidefinite certificate and its equality information. Earlier primary signed-set results include [Borg's 2007 compressed and uniform EKR cases](https://doi.org/10.37236/959) and [his large-alphabet theorem](https://doi.org/10.1016/j.disc.2008.09.039). The [general spectral Chvátal conjecture](https://arxiv.org/html/2609.28404v3#S5) is not resolved here.

AI-assisted mathematical research. The tensor transfer, signed eigenvalue bound and arithmetic controls were independently checked. No external peer review or priority claim.
