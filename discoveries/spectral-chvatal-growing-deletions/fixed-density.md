# Fixed positive density of triple deletions

**Theorem.** The $n\ge1000,\ D\le n^2/400$ range is sufficient. The pure inner equality holds for every simple 3-uniform deletion pattern in that range. The weighted residual admits individual real weights $0\le w_i\le d_i-d_0$ at a common marker and needs the additional **actual weight budget** $\sum_i w_i\le n^2/8$. The two-layer lift uses this budget at its crossing, so the sufficient condition for a nonempty anchor set $A$ is $|A|\min_{i\in A}(d_i-d_0)\le n^2/8$. No extra restriction applies to an empty anchor set. This is not general Conjecture H, arbitrary downsets, nor arbitrary face weights.

## 1. Fit and tight stars with slacked positive bridges

Use $K=\binom{n-1}{2}, a=\binom{n-2}{2}, b_{\rm bin}=\binom{n-3}{2}, q=\binom{n-4}{2}$, and put $M=K-d_0,\ t_i=K-d_i,\ r_i=d_i-d_0$, where $d_0=\min_i d_i$. All surviving triples have weight 1. Zero-weight bridges are removed. On the old incidence matrix $X$, whose columns are $x_i$, write
$$
 \lambda_i=\frac{t_i}{M-w_i}\in(0,1], \quad
 \alpha=\sum_i w_i,\quad L=M-\alpha .
$$
At $w_i>0$, fit $C_i$ on avoiding triples to $C_i^\top x_j=\lambda_j$ for $j\ne i$, and set it to zero on incident triples. Start at $1/a$ on avoiding rows. Their incidence Gram begins at $b_{\rm bin}I+(n-3)J$, loses row absolute sum at most $3D$, and has inverse row norm at most $2/(b_{\rm bin}-6D)$. The demand gap is a difference between the avoidance-deletion count divided by $a$ and $1-\lambda_j=(r_j-w_j)/(M-w_j)$, both nonnegative. In particular it is bounded by the **maximum**, not sum, of those quantities. For every slack site $\lambda_i<1$, fit $V_i$ identically to the earlier proof with $V_i=1$ on incident triples and $X^\top V_i=t_i e_i$. Crucially, this includes a site having a positive but slack bridge, not just unbridged sites.

Put $g=\mathbf1-Cw$, $F=\sum_i (1-\lambda_i)V_i V_i^\top/t_i$ over slack sites, and $B=C\,\operatorname{diag}(w)C^\top+gg^\top/L+F$. The exact substitution uses $C^\top x_i=\lambda_i(\mathbf1-e_i)$ when $w_i>0$ and $\lambda_i\mathbf1$ otherwise. Thus $g^\top x_i=L\lambda_i$, and
$$
 Bx_i=\lambda_i(\mathbf1-w_i C_i)+(1-\lambda_i)V_i .
$$
Here an absent column contributes its zero-weight term only. This equals 1 on rows incident at $i$. With old block $B+P$ and bridge diagonal $\operatorname{diag}(1/w)$, Schur complement is $P+gg^\top/L+F$. The normalized Gram sends a tight star (any $\lambda_i=1$) to $\mathbf1$. When its bridge is present, the old-block contribution of the bridge vertex is $w_iC_i$ and cancels the term $-w_iC_i$; the bridge output changes $C^\top x_i=\mathbf1-e_i$ to $\mathbf1$. There is always a tight star, since a minimum-degree site has $w_i=0,t_i=M$. Its total weight is $M$. Cauchy–Schwarz in the normalized Gram then gives $M H-J\succeq0$, the weighted theta dual certificate after diagonal scaling. No Gram invertibility is needed for this Cauchy–Schwarz step.

## 2. Compatible maximum-entry inverse and real PSD

As in the earlier audit use $P_*=I+A_{\rm disj}/q-J/K$. It annihilates full old incidence and has gap $\mu=1-1/q$ on its orthogonal complement. Let $P_s$ be its survivor restriction, $\Pi=X(X^\top X)^{-1}X^\top,\ E=I-\Pi,\ P_0=EP_s E$. Padding survivor vectors by zeros preserves the same gap on $\ker X^\top$. Let $\mathcal M$ preserve symmetric entries at nonempty intersection, including the diagonal. Call a masked symmetric matrix $Q$ compatible if $(QX)_{U,i}=0$ whenever $i\in U$. The true defect $Q_0=\mathcal M(I-B-P_0)$ is compatible. Its usable expansion includes the signed term from the previous audit:
$$
 Q_0=\mathcal M[-(B-J/M)+(P_s-P_0)-(1/M-1/K)J].
$$
Omitting that term is not harmless. For $T(Q)=\mathcal M(EQE)$, each active incidence sum is already on the mask, so $T$ preserves compatibility.

For column $V$, $z=X^\top Q_{\cdot,V}$ has zero entries inside $V$ by symmetry and compatibility. For $i\notin V$, exactly $K-q=3(n-3)$ complete old triples through $i$ meet $V$; surviving ones are no more numerous. Hence $\|z\|_\infty\le3n\|Q\|_{\max}$. The full old incidence Gram is $aI+(n-2)J$; the deleted row-sum norm is at most $3D$. Its survivor inverse row norm is at most $2/(a-6D)$, and $\|x_U^\top(X^\top X)^{-1}\|_1\le6/(a-6D)$. Therefore
$$
 \|\Pi Q\|_{\max}\le \frac{18n}{a-6D}\|Q\|_{\max}.
$$
For the extra right multiplication by $\Pi$, absolute column sums obey
$$
 \sum_U|\Pi_{UV}|
 \le \sum_i t_i|[(X^\top X)^{-1}x_V]_i|
 \le \frac{6K}{a-6D}.
$$
Thus $I-T$ has induced max-entry norm at most $\rho=(2+6K/(a-6D))18n/(a-6D)$ **on the compatible subspace**. No such estimate on arbitrary masked matrices was assumed: the zero inside $V$ is essential. As $a-6D\ge(12/25)n^2$, $K\le n^2/2$, and $n\ge1000$, $\rho\le99/320<1$. Consequently the real geometric inverse gives a unique compatible $Q=T^{-1}Q_0$ with $\|Q\|_{\max}\le\|Q_0\|_{\max}/(1-\rho)$. Using symmetry and at most $3K$ locked entries per row,
$$
 \|EQE\|_2\le \|Q\|_2\le 3K\|Q\|_{\max}.
$$
This last step bounds the sparse $Q$, not a presumed sparse $EQE$. The corrected $P=P_0+EQE$ satisfies all locks and annihilates $X$. Its PSD is a real operator-norm consequence of the strict gap; a modular fit or a finite partial Neumann sum is not the certificate.

## 3. Coefficients and bounds

In the given range $K,a,b_{\rm bin},q\ge(99/200)n^2$. Also $M,t_i,M-w_i,a-6D,b_{\rm bin}-6D\ge(12/25)n^2$ since $M-w_i\ge t_i$; $L\ge n^2/3$. The $C$ deviation from $1/M$ on avoiding rows uses $|M-a|/(Ma)\le(n+D)/(Ma)$, not $M\ge a$, which can fail in this regime. The rational coefficient bounds at $n^{-3}$ and $D n^{-4}$ are $1250/297<5$ and $71875/2376<32$, respectively. The $V_i$ avoidance bound has coefficient $25/2<13$ at $n^{-1}$.

With $W=C-J/M$, rank cancellation gives $B-J/M=W\operatorname{diag}(w)W^\top+(Ww)(Ww)^\top/L+F$. Write $c=5/n^3+32D/n^4$. The rank bounds are
$$
 3D/M^2+\alpha c^2,\quad
 (3D/M+\alpha c)^2/L,\quad
 \frac{3D+nD(13/n)^2}{(12/25)^2n^4}.
$$
For $\alpha c^2$ use $\alpha\le nD$ on the first squared term and $\alpha\le n^2/8$ on the second; $3D/M+\alpha c\le16D/n^2$. These give $R\le[846293/28800]D/n^4<30D/n^4$ for $D>0$. The prior exact first-order/second-order projection bounds remain valid with $H_0=3D+m\omega\le6D,\ \omega=6/(n-1)\le7/n,\ \delta=1/q-1/K\le13/n^3$, using $m\le Dn/3$. Their coefficient estimates are $215074/16335<14$ for $A$ and $83005/2916<30$ for $U_2$. The deleted principal block norm is at most the full Kneser residual norm $1+(n-5)/q<2$; using $1+m/q$ instead is unnecessary and could lose the new regime. The constant offset is at most $(1250/297)D/n^4<5D/n^4$ at $D>0$. Thus
$$
 \|Q_0\|_{\max}\le(30+2\cdot14+2\cdot30+5)D/n^4=123D/n^4.
$$
At $D=0$ all perturbation terms vanish. The exact uniform correction bound is $738/1105$ using $\rho\le99/320$; the slightly looser claimed $123/184$ also holds. Both are below $0.99<\mu$, so the real Schur-complement PSD argument applies.

## 4. Full weights and marker

Let $S_0=n+K$, $b=\min_{i\in A}r_i$ for nonempty $A$. The unbridged minimum old star and minimum-surplus anchor are the two possible maximizing branches; if a minimum is itself anchored, $b=0$. Pure inner is the zero-bridge residual plus full lower ranks, value $S_0-d_0$. Pure outer is the earlier coloring and full triple layer, value $S_0+1$ if $A$ is nonempty, $S_0$ otherwise. At $b>0$, use $w_i=b$ on $A$ **directly**, not the full surpluses. The crossing $(u,v)=(1,b)$ then combines the residual with $b$ full triple layers and $b+1$ full lower layers for cap $s=(b+1)S_0-d_0$. The marker has total star weight $(|A|+1)b\le(1/8+1/400)n^2<s$. In the dual Gram the bridge block is $tI-J$, $t=s/b$. The squared length of the marker's forced projection is $|A|/(t-|A|)\le t-1$; complete with an independent vector to achieve norm squared $t-1$. All other old faces are disjoint from the marker and impose no fixed dual inner products. Interpolate the crossing and appropriate pure endpoint by positive homogeneity and subadditivity of weighted theta. Minimum-degree anchors at $b=0$ are common endpoint maximizers without a crossing; with no anchors both pure endpoints share a minimum-degree old center after ignoring the universal marker. The exact formula is $S_0(u+v)-d_0u+\max(0,v-bu)$ when $A\ne\varnothing$, omitting the max term when empty. The bridge budget is sufficient, not asserted necessary.

## 5. Independent replay and near-saturation product

Run `node fixed-density-replay.cjs`. The [replay](fixed-density-replay.cjs) has no external imports or data. The [expected transcript](fixed-density-expected.json) records the complete result.

* A small $n=12$ actual all-touched product-union fixture has fourteen deletions, minimum degree 3, both a tight bridge of weight 1 and a positive **slack** bridge of weight $1/2$. Exact rational Gauss–Jordan with actual avoiding incidences verifies 84 column-fit equations, 2,472 coupled coordinates and all 618 actual mask locks. It constructs the projected Kneser block and true compatible defect, and omitting its offset $3/2860$ fails all 618. For direct max-entry control an independently chosen signed four-cycle on four star triples with disjoint noncentral extras is a nonzero compatible $Q$, $\|Q\|_{\max}=1$. Exact projected products give $\|(I-T)Q\|_{\max}=19173412/207645525$; the nonzero first-/second-order products and maximum absolute column sum are independently tabulated. This is a diagnostic control, not a theorem at $n=12$.
* A near-saturating dense family: let $n=60s$, $s\ge17$, and form three bands of twenty disjoint blocks of $s$ old points each. For each $a\bmod20$ and $i,j\in\{0,1,2\}$, delete the complete three-way product $A_a\times B_{a+i}\times C_{a+j}$, except restore one whole block product. The 180 block triples are distinct and every block begins in nine products. Thus $m=179s^3$, $D=9s^2=n^2/400$, $d_0=8s^2>0$, with $57s$ old vertices of positive surplus $s^2$. Choosing 450 of them as anchors gives $b=s^2$ and $|A|b=n^2/8$ exactly. At $s=17$ this is 879,427 deletions out of 176,348,140 triples, $179/180$ of maximum degree-budget saturation; $M<a$ by 1,294. The exact crossing value is 150,708,078. The analytic and exact rational checks certify the stated regime without materializing the large Gram. The limiting deletion fraction is $179/36000$, approximately 0.49722%, close to the degree-budget ceiling 0.5%. This is fixed positive density, not density approaching one or all rank-three downsets.

The [linear-degree theorem](all-touched.md) used degree linear in $n$ and the older Frobenius estimate. The exact old rank, incidence, and marker lemmas are inherited rather than new. [Ellis–Filmus–Friedgut §5](https://arxiv.org/html/2609.28404v3#S5) distinguishes general projection-packing tightness from the proposed spectral H equality. The general Conjecture H is not settled. No priority claim or external peer review.
