# Weighted theta for independently weighted triples and multiple markers

**Theorem.** For $n\ge1000$, old singleton/pair weights one, arbitrary old triple weights $p_U\in[0,1]$, and arbitrary finitely many new markers with independently weighted singletons and old–new pairs in $[0,1]$, define
$$
 K=\binom{n-1}{2},\quad d_i=\sum_{U\ni i}(1-p_U),\quad D=\max d_i,\quad
 w_i=\sum_{e\text{ bridge at }i}w_e,\quad W=\max w_i,\quad
 B_*=D+W,\quad\alpha=\sum_i w_i .
$$
If $B_*\le n^2/400$ and $\alpha\le n^2/8$, weighted theta equals $n+\max_i(K-d_i+w_i)$. There are no new–new pairs or new triples. Zero weights are omitted. This is not arbitrary weights on all faces, an unrestricted residual with arbitrarily heavy bridges, arbitrary downsets, or general Conjecture H.

## 1. The heterogeneous seed preserves the full gap

Work temporarily on the complete triple layer, $X$ its old-point incidence. Put $F=\mathrm{diag}(\sqrt p)$, $f=F\mathbf1$, $\Delta=I-F^2$, $\widehat X=FX$. The complete Kneser residual is $P_*=I+A_{\rm disj}/q-J/K$, $q=\binom{n-4}{2}$. It kills $X$; its exact eigenvalues off incidence are $1+(n-5)/q$ and $\mu=1-1/q$. Set $L_* = I-P_*$ and
$$
 R=F P_*F+\Delta=I-F L_* F.
$$
This is PSD on the whole space, and its locked entries are $I-ff^\top/K$ since $P_*=I-J/K$ there. For $y\in\ker\widehat X^\top$, $Fy\in\ker X^\top$, so $y^\top R y\ge\mu\|Fy\|_2^2+y^\top\Delta y\ge\mu\|y\|_2^2$. There is no weighted/unweighted subspace angle loss. A zero-weight triple has $R_{UU}=1$, its other entries zero, and a zero incidence row.

Because $L_* X=X$, multiplication gives
$$
 R\widehat X=F L_*\Delta X,\qquad
 \widehat X^\top R\widehat X=X^\top\Delta X-X^\top\Delta L_*\Delta X .
$$
The first identity cancels the entire bare diagonal contribution. The second comes from replacing $F^2$ by $I-\Delta$ on the remaining side, not from bounding the two terms independently without signs. All entries of $L_*$ have magnitude at most $1/q$. Thus $\|F L_*\Delta X\|_{\max}\le D/q$. For $G=\widehat X^\top\widehat X$, relative to $aI+(n-2)J,\ a=\binom{n-2}{2}$, the Gram loss is $X^\top\Delta X$ with absolute row-sum norm at most $3D$. The inverse row norm is $\le2/(a-6D)$, and least eigenvalue of $G$ is $\ge a-3D$. Therefore for $\Pi=\widehat XG^{-1}\widehat X^\top$, $\|\Pi R\|_{\max}\le6D/[(a-6D)q]$. For the quadratic term the middle matrix above is PSD; $\|\Delta X\|_2^2\le\|X^\top\Delta X\|_2\le3D$ because $0\preceq\Delta^2\preceq\Delta$, and $\lambda_{\min}(L_*)=-(n-5)/q$. Hence its norm is $\le3D[1+(n-5)/q]$. A row $\hat x_U^\top G^{-1}$ has Euclidean norm $\le\sqrt3/(a-3D)$; it follows that $\|\Pi R\Pi\|_{\max}\le9D[1+(n-5)/q]/(a-3D)^2$. Using $q,a\ge(99/200)n^2,\ a-6D\ge(12/25)n^2,\ a-3D\ge(39/80)n^2$, the exact coefficients are $2500/99<26$ and $634880/16731<39$ at $D/n^4$.

## 2. Rank fits, bridge multiplicities, and signed shift

Let $t_i=K-d_i,\ M=\max_i(t_i+w_i)$, $\lambda_i=t_i/(M-w_i)$, and $L=M-\alpha$. Then $0<\lambda_i\le1$ and $|M-K|,\,M-w_i-t_i,\,w_i\le B_*$. Each of $M,t_i,M-w_i$ is at least $(12/25)n^2$; $L\ge n^2/3$. The nonnegative weighted losses in avoiding-row Gram sums are bounded by $3D$ as before. At a positive aggregate bridge fit a physical $C_i$, zero on incident triples, to $X^\top\operatorname{diag}(p)C_i$ equal $\lambda_j$ at $j\ne i$, and 0 at $i$. Starting at $1/a$, the adjustment demand is a difference between two nonnegative quantities, a weighted avoidance deficit divided by $a$ and $1-\lambda_j$; both are bounded by $B_*/[(12/25)n^2]$. Fit $V_i$ at every positive slack, even if it has a positive bridge, to $t_i e_i$ with value one on incidence. Its off-incidence bound is $13/n$. In the $C$ estimate use $|M-a|\le n+B_*$, which permits $M>K$ as well as $M<K$; the bound off incidence is $c=5/n^3+32B_*/n^4$.

Use hats to multiply fitted columns by $F$. Set $g=f-\widehat C w$, $F_{\rm slack}=\sum_i[(1-\lambda_i)/t_i]\widehat V_i\widehat V_i^\top$, and $B_{\rm rank}=\widehat C\operatorname{diag}(w)\widehat C^\top+gg^\top/L+F_{\rm slack}$. Exact multiplication gives $g^\top\hat x_i=L\lambda_i$ and $B_{\rm rank}\hat x_i=\lambda_i(f-w_i\widehat C_i)+(1-\lambda_i)\widehat V_i$, which is $f_U$ on each incidence. For actual bridge $e=\{i,b\}$ use the same column $\widehat C_i$ and reciprocal diagonal $1/w_e$. Put **all** distinct bridge-block entries zero. For intersecting pairs this is required; for disjoint pairs it is an allowed free choice. Its Schur subtraction is $\sum_e w_e\widehat C_{i(e)}\widehat C_{i(e)}^\top=\widehat C\operatorname{diag}(w)\widehat C^\top$, not weights squared. For an old tight star, all its bridge contributions together cancel $-w_i\widehat C_i$ on the old rows. Each bridge at the tight coordinate gets exactly one from its own diagonal and zero from the fit; all other bridges get one from the fit. This holds with two or more tight bridges at the same old coordinate.

Exact first-order rank cancellation centers $B_{\rm rank}$ at $ff^\top/M$. The three entrywise rank bounds are $3B_*/M^2+\alpha c^2,\ (3B_*/M+\alpha c)^2/L,\ [3B_*+nB_*(13/n)^2]/[(12/25)^2n^4]$. Use both $\alpha\le nB_*$ and $\alpha\le n^2/8$ to obtain $\|B_{\rm rank}-ff^\top/M\|_{\max}\le30B_*/n^4$, by the previously audited rank calculation. With $E=I-\Pi,\ P_0=ERE$, the signed compatible defect on the mask is
$$
 Q_0=\mathcal M\big[-(B_{\rm rank}-ff^\top/M)+(R-P_0)+(1/K-1/M)ff^\top\big].
$$
The sign changes with $M-K$; replacing $K$ by $M$ only in the seed would not justify the exact projection identities. The magnitude of this last term is at most $5B_*/n^4$. Omitting it changes the required row incidence at $i\in U$ by $f_U t_i(1/M-1/K)$. Thus it is nonzero on every positive-weight locked incidence whenever $M\ne K$.

## 3. Exact inverse and star certificate

The mask preserves entries on intersecting triples and the diagonal. Compatibility is $(Q\widehat X)_{U,i}=0$ for $i\in U$; it follows for $Q_0=\mathcal M(I-B_{\rm rank}-P_0)$ from the exact locks. The operator $T(Q)=\mathcal M(EQE)$ preserves compatibility. As checked for the earlier weighted crossing, an outside incidence sum in a compatible column has at most $3(n-3)$ terms; the three inside sums are zero. Weighted incidence entries are at most one, so $\|\Pi Q\|_{\max}\le18n\|Q\|_{\max}/(a-6D)$, and absolute column sums of $\Pi$ are $\le6K/(a-6D)$. The contraction $\|I-T\|_{\max\to\max}\le99/320$ holds on the compatible subspace. Every zero-weight triple has zero row in $\Pi$ and $B_{\rm rank}$, and unit isolated diagonal in $R,P_0$; its defect row is identically zero. One can delete the zero rows before inversion or restrict the full inverse to zero-padded corrections. Both agree and neither takes $1/0$.

The full defect is at most $(30+5+2\cdot26+39)B_*/n^4=126B_*/n^4$. With $Q=T^{-1}Q_0$, the support of $Q$ has at most $3K$ entries per row, so $\|EQE\|_2\le3K\|Q\|_{\max}\le756/1105<0.685<\mu$. Thus $P=P_0+EQE$ has the locked entries and is PSD. The bridge Schur complement is $P+gg^\top/L+F_{\rm slack}$ and is PSD too. Scale each positive triple row and column by $1/f_U$ to obtain the physical reciprocal diagonal $1/p_U$; omit zeros rather than invert. The tight-star calculation returns $\mathbf1$, after rescaling, on the physical Gram. This is the usual weighted Cauchy–Schwarz certificate for theta $M$ on triples and bridges, attained by a maximizing old star.

## 4. Lower ranks, independent singletons, and older checkpoints

Full old singletons and pairs add exactly $n$ to the cap and every old star. For odd $n$, use pair endpoint sum modulo $n$ and put each singleton in its missing color. For even $n$, color pairs by a one-factorization and put all singletons in a separate class. Color classes are disjoint-face cliques in the disjointness graph. The weighted theta upper bound is subadditive; an exact old star attains the $s=n+M$ bound before new singletons.

In a cap-$s$ dual Gram, the bridge metric is $s\operatorname{diag}(1/w_e)-J$ irrespective of shared old anchors. A newcomer $b$ of singleton weight $u_b>0$ needs inner product $-1$ solely with its own bridge vectors. Writing $q_b=\sum_{e\ni b}w_e$, its needed projection is a linear combination with coefficients $-w_e/(s-q_b)$ on those bridges. These coefficients give the desired $-1$ and squared norm $q_b/(s-q_b)$; its full dual squared norm $s/u_b-1$ is at least this since $q_b+u_b\le n+1<s$. Supply an independent orthogonal residual for each new singleton. The resulting inner products with other marker singletons and their bridges need not be zero: those faces are disjoint even if their bridges share old anchors, so their dual entries are **unlocked**. A bridge-free singleton has zero projection, and a zero singleton is simply omitted. Every new-center star is at most $n+1$; the largest old star still supplies the matching lower bound.

In particular, the theorem includes common restored-triple weight $r\in[0,1]$ with maximum deletion degree at most $n^2/450$ and at most $n/400$ bridge weight per old point: $B_*\le(1-r)n^2/450+n/400\le n^2/400$ and $\alpha\le n^2/400$. An alternate earlier split into two ranges of $r$ also works; the replay retains its constant checks, but it is not needed for this theorem.

**Boundary distinction.** The narrower published common-r one-marker two-layer normalization is also contained **at its boundary**: if the original deletion maximum is $D_0\le n^2/400$, set $r=v/(u+v)$ for nonzero total. Effective $D=(1-r)D_0$, $W\le r$, so $B_*\le(1-r)D_0+r\le\max(D_0,1)\le n^2/400$ and $\alpha\le n$. The independent earlier weighted survivor-only residual has different hypotheses, including potentially larger individual bridge weights or $D+W$; it is not claimed subsumed. Nor does the new theorem cover every independently weighted one-marker model at $D=n^2/400$ with arbitrary positive $W$.

## 5. Independent exact controls

Run `node heterogeneous-replay.cjs`. The [standalone source](heterogeneous-replay.cjs) and [expected transcript](heterogeneous-expected.json) require no external dependencies. The fixture uses nine old points and old weights $0,1/9,1/4,4/9,9/16,1$, two zero-weight triples, six positive bridges over three markers, a fourth positive bridge-free singleton, and omitted zero bridge/singleton. Tight old site 8 has two bridges; slack site 6 has two. It tests both signs by rescaling bridges from 1 to $1/8$, giving $M-K=31/36$ and $-109/288$. These are independent exact algebra controls.

* Seed: 756 exact first-order identities, 81 quadratic incidence entries, 5,376 locked entries, 756 null entries of $P_0$, and 168 exact zero-row comparisons. An explicit heterogeneous weighted-complement trade plus one zero-row direction gives rayleigh numerator $8>36/5$, the exact gap lower bound.
* Each rank/marker control: 108 avoiding-fit equations, 246 actual positive-row compatibility locks (all fail if the signed offset is dropped), 88 tight outputs, 164 duplicated-column comparisons, ten disjoint bridge pairs, six singleton locks and twelve nonzero allowed cross-bridge entries. The computed cross-marker singleton projections are generally nonzero even with shared old anchors; no such entry is locked.
* The small fixtures prove identities, not real PSD at $n=9$. The same positive-center finite pattern padded to $n=1000$ is in theorem range without changing $D,W,\alpha$ or the strict maximizing center; its exact theorem-domain theta is $17982067/36$.

Closest scope references are the [one-marker proof](weighted-crossing.md), the previously audited weighted rank/mask construction, and [Ellis–Filmus–Friedgut §5](https://arxiv.org/html/2609.28404v3#S5). The new seed and multiple independently weighted markers are the scoped additions, not a priority or general H claim.
