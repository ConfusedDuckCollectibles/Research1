# Unequal deficits in weighted transversal products

This companion proves a restricted transversal-family theorem for weighted Lovász theta. It does not settle spectral Chvátal H in general or establish a perfect-matching factorization.

## 1. Statement and normalization

Let \(r\ge3\) and \(s\ge400r^3\) be integers. There are \(r\) disjoint parts of size \(s\), and the full top product \(\mathcal T\) consists of the \(r\)-tuples with one label in each part. Write \(K=s^{r-1}\). The initial weight at each top tuple is any \(p_U\in[0,1]\), including zero. For each point \(i\) let
\[
d_i=\sum_{U\ni i}(1-p_U),\qquad D=\max_i d_i,\qquad m=\min_i d_i.
\]
If \(D\le K/(40r^3)\), the weighted Lovász theta of the graph on positive top weights with edges joining disjoint tuples is exactly
\[
M=K-m .
\]
If every proper nonempty transversal is also included at unit weight, the value is exactly
\[
M+\sum_{k=1}^{r-1}\binom{r-1}{k-1}s^{k-1}.
\]
There is no balance, deficit-spread, or individual balancing-addition assumption. Hypotheses are near-completeness within an equal-part product, not near-completeness of all \(r\)-sets on the ground set. No newcomer markers or independently weighted lower faces are included. Binary weights with all lower faces yield an ordinary downset on adding the empty face, which is a universal weight-one vertex and leaves the formula unchanged.

In the weighted theta primal, \(Y\succeq0\), \(\operatorname{tr}Y=1\), and \(Y_{UV}=0\) at disjoint pairs; maximize \(\langle \sqrt p\sqrt p^\top,Y\rangle\). If there is a PSD \(P\) whose diagonal and intersecting entries equal those of \(I-\sqrt p\sqrt p^\top/M\), then testing \(\langle P,Y\rangle\ge0\) gives the upper bound \(M\). The star of any point with deficit \(m\) attains it by the rank-one feasible matrix made from \(\sqrt p\) restricted to that star. This uses no unverified theta monotonicity convention.

## 2. Balance by adding mass, with no bound on a single output weight

Put \(b_i=d_i-m\). Every part has the same \(\sum d_i\) because each original deficit \(1-p_U\) contributes once per part. Therefore \(T=\sum b_i\) on any one part is independent of the part. If \(T=0\) let \(c_U=0\); otherwise take
\[
c_U=\frac{\prod_{i\in U}b_i}{T^{r-1}} .
\]
Summing over all but one coordinate gives \(\sum_{U\ni i} c_U=b_i\). This is just the independent product coupling of these prescribed added marginals. The proof would work with any nonnegative coupling with the same marginals.

Set \(w_U=p_U+c_U\). Its point marginals are all \(M\), its total mass is \(sM\), and \(\sum_{U\ni i}c_U\le D\). Define \(D_a=\max_i\sum_{U\ni i}|1-w_U|\). The triangle inequality gives \(D_a\le\max_i(d_i+b_i)\le2D\). The completed weights are allowed to exceed one. The signed loss \(\Delta=I-\operatorname{diag}(w)\) must **not** be treated as PSD.

It suffices to construct the PSD certificate for \(w\) with \(\sqrt w=f\). In fact, restrict that certificate to indices with original weight \(p_U>0\), let \(S=\operatorname{diag}\sqrt{p_U/w_U}\) there, and use \(SPS+(I-S^2)\). Since \(0\le S^2\le I\), this is PSD and has the exact original locks. This is a direct diagonal congruence-plus-slack proof of the needed monotonicity, including newly restored zero rows and arbitrarily large \(w_U\).

## 3. Signed overlap supplies an ordinary spectral gap

Let \(X\) be full tuple-by-point incidence, with row \(x_U^\top\), and \(q=(s-1)^{r-1}\). Let \(L=J/K-A/q\), where \(A\) is the adjacency of disjoint full tuples. The kernel entries are \(1/K-\mathbf1[U\cap V=\varnothing]/q\), and \(LX=X\): there are \(q\) disjoint tuples through a point outside \(U\) and none through a point in \(U\). The projector \(H\) onto \(\operatorname{im}X\) has entries
\[
H_{UV}=\frac{h}{K}-\frac{r-1}{sK},\qquad h=|U\cap V|.
\]
Indeed, coordinate incidence spans the constants and the \(r\) mutually orthogonal single-coordinate mean-zero spaces. The constant projector has entry \(1/(sK)\), and the mean-zero projector from one part has entry \((\mathbf1[U_i=V_i]-1/s)/K\). Consequently
\[
B_{UV}:=(L-H)_{UV}=\frac{1-h}{K}+\frac{r-1}{sK}-\frac{\mathbf1[h=0]}q .
\]

Denote \(F=\operatorname{diag}\sqrt w\), \(\widehat X=FX\), and \(R=I-FLF\). Work on the completed positive support when a reciprocal or relative norm is used; pad by zeros to evaluate the original incidence and \(L\) formulas. For \(y\in\ker\widehat X^\top\), \(Fy\perp\operatorname{im}X\) and \(y^\top R y=\|y\|^2-y^\top FBFy\). The Schur bound with positive vector \(f=\sqrt w\) on this support gives
\[
\|FBF\|_2\le\max_U\sum_V |B_{UV}|w_V .
\]
It follows from the standard symmetric weighted row test: \(|y_U y_V|\le\frac12(y_U^2 f_V/f_U+y_V^2 f_U/f_V)\), so it applies without a positive-semidefiniteness condition on \(B\).

For \(h\ge1\), bound \(|B_{UV}|\le[(h-1)+(r-1)/s]/K\). For original weights \(p\le1\) the sum of \((h-1)\) over intersecting \(V\) is at most \(\binom r2 K/s\), by \(h-1\le\binom h2\) and counting \(V\) agreeing on a chosen pair of parts. For added weights use \(h-1\le h\) and \(\sum_V c_V h=\sum_{i\in U} b_i\le rD\). The total intersecting mass is at most \(rM\).

For \(h=0\), Taylor's formula for \((1-t)^{-(r-1)}\) at \(t=1/s\) gives
\[
0\le K/q-1-(r-1)/s\le r(r-1)/s^2.
\]
Its second derivative multiplier beyond \(r(r-1)\) is at most \((1-1/s)^{-(r+1)}<2\), by Bernoulli and the size cutoff. The disjoint weights sum to at most \(sM\). Combining,
\[
\|FBF\|_2\le \frac{(5/2)r(r-1)}s+\frac{rD}K
\le \frac{r+3}{160r^2}\le\frac1{240}.
\]
In particular \(R\) has quadratic-form gap at least \(239/240\) on \(\ker\widehat X^\top\). This is a bound on the compression only; \(R\) can be indefinite on the incidence space. The signed cancellations are taken before absolute values. The estimate is polynomial in \(r\) and independent of \(\max w_U\); the old sector-by-sector \(2^r\) bound is unnecessary.

## 4. Pseudoinverse and relative entrywise estimates

Let \(G_f=X^\top X\), \(G=\widehat X^\top\widehat X=G_f-X^\top\Delta X\), and \(Z\) be the part-constant column vectors whose part constants sum to zero. The full Gram has diagonal blocks \(KI_s\) and cross blocks \((K/s)J_s\); its kernel is \(Z\), and its nonzero eigenvalues are \(K,rK\). In particular
\[
G_f^+=\frac{I-\operatorname{blockdiag}(J_s/s)}K+\frac{J}{r^2sK},\qquad
\|G_f^+\|_\infty\le3/K,
\]
where \(\|\cdot\|_\infty\) here is maximum absolute row sum. Every row of \(X\) kills \(Z\). Even with signed \(\Delta\), \(X^\top\Delta X\) has maximum absolute row and column sum at most \(rD_a\). Thus \(G\) has exactly the same kernel, is positive on \(Z^\perp\) with eigenvalues at least \(K-rD_a\), and the Neumann expression \((I-G_f^+X^\top\Delta X)^{-1}G_f^+\) equals \(G^+\). Its absolute row-sum norm is at most \((3/K)/(1-3rD_a/K)\le4/K\). This explicitly inverts the rank-deficient incidence Gram on its support; it does not assume the columns of \(X\) are independent.

Let \(\Pi=\widehat XG^+\widehat X^\top\) and \(E=I-\Pi\) (ordinary orthogonal projection). On positive rows put \(\|Q\|_*=\max_{U,V}|Q_{UV}|/(f_U f_V)\). Bernoulli gives \(q\ge.99K\). The two exact identities
\[
R\widehat X=FL\Delta X,\qquad
\widehat X^\top R\widehat X=X^\top\Delta X-(\Delta X)^\top L(\Delta X)
\]
follow by \(LX=X\), \(X^\top L=X^\top\) and \(F^2=I-\Delta\). They hold with signed losses and zero rows. For the first one, \(\lvert(L\Delta X)_{U,i}\rvert\le D_a/q\) since \(|L_{UV}|\le1/q\). In \(\Pi R\), dividing by \(f_U f_V\) leaves one row \(x_U^\top G^+\) and the bounded \(L\Delta X\) entries. Therefore \(\|\Pi R\|_*\le4rD_a/(Kq)\le5rD_a/K^2\).

For the projected double product, the relative expression sandwiches the second identity between \(x_U^\top G^+\) and \(G^+x_V\). These rows have 2-norm at most \(\sqrt r/(K-rD_a)\) and 1-norm at most \(4r/K\). The first term has operator norm \(\le rD_a\), so contributes at most \(2r^2D_a/K^2\). The second has maximum absolute entry \(\le D_a^2/q\) by bounding the two signed column \(\ell_1\) norms before summation, so contributes at most
\[
(4r/K)^2 D_a^2/q\le r^2D_a/K^2
\]
since \(16D_a/q\le1\) at the stated cutoff. This step requires neither sign positivity nor a bound on \(\Delta^2\). The independent replay exhibits a negative direction for the full matrix \(\widehat X^\top R\widehat X\), so importing its PSD from the balanced \([0,1]\)-weight special case would be invalid.

## 5. Compatible matrix inversion and recovery of the original certificate

Let \(\mathcal M\) retain the diagonal and distinct intersecting entries. A symmetric matrix \(Q\) supported on this mask is compatible if \((Q\widehat X)_{U,i}=0\) whenever \(i\in U\). Define
\[
Q_0=\mathcal M(I-ff^\top/M-ERE).
\]
It is compatible because masking leaves the relevant incidence sums unchanged, \(ERE\widehat X=0\), and \(f^\top\widehat X_i=M\). The diagonal \(I\) cancels before taking any relative norm: \(\mathcal M(R)=\mathcal M(I-ff^\top/K)\). Expanding \(ERE=R-\Pi R-R\Pi+\Pi R\Pi\), using \(|1/K-1/M|\le2D/K^2\), and substituting \(D_a\le2D\) gives
\[
\|Q_0\|_*\le(20r+6r^2+2)D/K^2\le13r^2D/K^2 .
\]

For compatible \(Q\), \((\widehat X^\top Q)_{i,V}/f_V\) is zero at \(i\in V\). Outside \(V\), its absolute value is at most \(\|Q\|_*[(r-1)K/s+D]\): among original \(p\le1\) there are at most \(r-1\) choices of the second, intersecting coordinate (its part must differ from that of \(i\)), and through two prescribed points there are \(K/s\) tuples. All added mass through \(i\) is at most \(D\). It follows that
\[
\|\Pi Q\|_*\le4r[(r-1)/s+D/K]\|Q\|_* .
\]
Also \(\sum_U|\Pi_{UV}|f_U/f_V\le M\|G^+x_V\|_1\le4r\); this is a **weighted** absolute column sum. It controls the second factor in \(\Pi Q\Pi\) in the same relative norm. The map \(T(Q)=\mathcal M(EQE)\) preserves symmetry and compatibility, because \(EQE\widehat X=0\) before masking. Thus on the compatible space,
\[
\|T-I\|_{*\to *}\le (2+4r)4r[(r-1)/s+D/K]
\le \frac{(2r+1)(r+9)}{50r^2}\le\frac{14}{75}<1/4.
\]
In the penultimate inequality the numerator after division by \(r^2\) is \(2+19/r+9/r^2\), decreasing for \(r\ge3\). This gives the exact inverse via a convergent Neumann series on the finite-dimensional compatible space, not an assertion of invertibility on arbitrary matrices.

With \(Q=T^{-1}Q_0\), set \(P=ERE+EQE\). For this symmetric masked \(Q\), the same weighted Schur test gives \(\|Q\|_2\le rM\|Q\|_*\) since each intersecting row mass is \(\le rM\). Compression may be dense; it is bounded by \(\|EQE\|_2\le\|Q\|_2\), not from a claim that \(EQE\) is itself masked. Hence
\[
\|EQE\|_2\le \frac{13r^3D}{K(1-14/75)}
\le\frac{195}{488}<\frac{239}{240}.
\]
Thus \(P\succeq0\) and \(P\widehat X=0\). By definition of \(Q_0\) it has exactly the locks \(I-ff^\top/M\). Rows where \(w=0\) are an isolated identity block if retained temporarily; projecting to positive support loses no contribution. Coming back to \(p\) uses the diagonal slack in §2. At \(D=0\), all original weights are one and no correction is needed; \(T=0\) simply means the original deficits were already equal. All numerical cutoffs permit equality in the allowed deficit bound.

## 6. Complete proper layers and attribution

For a proper \(k\)-transversal let \(S\) be its \(k\)-part support. Set \(b=\binom{r-1}{k-1}\), \(c=\binom{r-2}{k-2}\) (zero for \(k=1\)), \(\epsilon=b/[s(b-c)+kc]\), and \(\gamma=\epsilon(b-c)\). For disjoint \(U,V\) with supports \(S,T\), put
\[
A'_{UV}=\frac{1-\epsilon|S\cap T|}{(s-1)^{|S\cap T|}s^{k-|S\cap T|}},
\]
and zero otherwise. The unweighted \(H_k=I+A'/\gamma\) maps every point-star indicator \(z_i\) to \(\mathbf1\). Indeed, outside a given \(U\), the unnormalized adjacency sum through a specified point is either \((b-\epsilon kc)/s=\gamma\) or \((b-\epsilon[b+(k-1)c])/(s-1)=\gamma\), according to whether that point's part is absent or present in \(U\). If the point is in \(U\), the identity matrix provides the required one.

For PSD decompose into tensor label modes. The support matrix is \(W_{ST}=1-\epsilon|S\cap T|\ge0\) since \(s\ge k\). In the constant-label sector \(W=J-\epsilon B_0B_0^\top\), where \(B_0^\top B_0=(b-c)I+cJ\) for the \(k\)-support incidence; its least eigenvalue is \(-\gamma\). A fixed choice of \(j\ge1\) mean-zero label coordinates gives the corresponding support principal submatrix of \(W\) scaled by \((-1)^j/(s-1)^j\), since disjoint label averaging annihilates any mean-zero mode in an unshared part. The restricted \(W\) has row sums at most \((s-1)\gamma\), by restricting only to a single fixed shared mode. Hence all of \(A'\) has spectrum bounded below by \(-\gamma\). Thus \(H_k\succeq0\), including \(H_1=J\). Gram Cauchy–Schwarz with \(z_i\) gives \(H_k-J/M_k\succeq0\) for \(M_k=b s^{k-1}\). No divisibility assumption appears.

For the top layer, conjugating its normalized certificate by \(\operatorname{diag}(1/\sqrt{p})\) and adding \(J/M\) produces its physical \(H_{\rm top}\). For the block sum of the top and proper-rank \(H\) matrices, Cauchy–Schwarz gives \(\sum_j(\mathbf1^\top y_j)^2/M_j\ge(\sum_j\mathbf1^\top y_j)^2/\sum_jM_j\). Cross-layer intersections have zero Gram entries as required, and the same old point with deficit \(m\) attains all layers. The lower-rank certificate and balanced predecessor are credited to the [pinned public balanced companion, commit 29145a](https://github.com/ConfusedDuckCollectibles/Research1/blob/29145a4d3ed0f154e051fa6071a2e54c01744b60/discoveries/spectral-chvatal-growing-deletions/balanced-transversals.md), fetched directly at that pin. The latter has better parameters for equal deficits (\(s\ge100r^3,\ d\le K/(10r^3)\)); this extension does not replace that stronger special-case bound.

## 7. Exact controls, concentration and literature scope

Run `node unequal-transversal-replay.cjs` and compare its output to [the expected JSON](unequal-transversal-expected.json). The fresh independent fractional \(r=3,s=3\) fixture deletes the diagonal matching and changes original tuple weights \(001,010,100\) to \(1/4,1/2,0\). It has four original zeros, one zero after completion, \(D=11/4,m=1,M=8\) and maximum \(w=34/27\). Direct rational lifted Gauss-Jordan reconstructs both pseudoinverses and 729 unsigned overlap entries, and checks 243 signed-annihilation entries, 676 weighted-projector entries, both signed identities, and 156 compatible lock positions. A second \(r=3,s=3\) fixture deleting only \(000,111\) has \(m=0,w_{\max}=5/4\); the incidence test vector \(e_1-e_0\) gives \(-3/8\) on \(\widehat X^\top R\widehat X\), a negative control for the tempting but invalid PSD shortcut.

The tiny \(r=3,s=2\) end-to-end control solves the exact linear compatible matrix equations after congruencing away all square roots. Its 32 allowed symmetric entries have a combined system of rank 31 (outside the theorem's cutoff, no uniqueness is claimed). The completed and original-support certificates each have exact PSD rank 3 by rational Schur elimination. One original zero is restored in completion, then discarded when recovering the original certificate. A separate balanced fixture exercises \(T=0\). These are precise miniature checks, not proofs of theorem-domain positivity.

For a theorem-domain family, delete binary tuples with **exactly one zero**, where the remaining labels (viewed bijectively modulo \(s-1\)) sum to one of \(\ell\) selected residues. A zero point loses \(D=\ell(s-1)^{r-2}\), a nonzero point loses \(m=\ell(r-1)(s-1)^{r-3}\). Hence all additional mass goes to the all-zero tuple, with final weight \(1+\ell(s-r)(s-1)^{r-3}\). Example \(r=3,s=10801,\ell=5\) has \(D=54{,}000,m=10,M=116{,}661{,}591\) and a completed weight \(53{,}991\). A second example, \(r=13,s=900000,\ell=9\), exercises the same polynomial budget at higher rank. Large examples verify integer formulas and cutoffs without materializing massive Grams. Even zero partite codegree can occur under unequal deficits: remove the \(s\) diagonal-translation matching classes that meet one fixed admissible \((r-1)\)-prefix, then delete one remaining tuple; all coordinate deficits are \(s\) or \(s+1\), within the budget at the stated size.

The exact complete signed-set independence number is classical, as in [Borg and the cited Deza–Frankl/Bollobás–Leader results](https://arxiv.org/abs/1912.10324). The general downset spectral conjecture is [Ellis–Filmus–Friedgut §5](https://arxiv.org/html/2609.28404v3#S5); this is not its resolution. A graph fractional clique decomposition [arXiv:2604.25206v2](https://arxiv.org/abs/2604.25206) concerns graph edges, not an SDP on weighted transversal hyperedges. [Keevash–Mycroft](https://people.maths.ox.ac.uk/keevash/papers/geom-match-journal.pdf) and [Counting Perfect Matchings in Dirac Hypergraphs](https://link.springer.com/article/10.1007/s00493-025-00194-8) distinguish fractional vertex matching (satisfying point totals) from a decomposition into perfect matchings with prescribed tuple marginals. [Cheng–Glock–Lane Theorem 1 and Conjecture 49](https://arxiv.org/html/2611.27781v1#S10) concern nonpartite integral factorization and a separate high-codegree multipartite analogue, respectively. These do not supply this weighted theta assertion at possible partite codegree zero. No identical primary theorem was found in the bounded search, not a priority claim. Classical matching arguments prove some special controls and are not challenged by this theorem.

AI-assisted mathematical research with an independent deductive audit and exact arithmetic checks. Not externally peer reviewed. No priority claim.
