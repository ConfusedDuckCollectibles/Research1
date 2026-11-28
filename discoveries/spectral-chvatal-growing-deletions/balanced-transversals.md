# Balanced weighted transversals and complete lower layers

This companion proves a top-layer weighted-theta theorem for balanced perturbations of a product. A direct spectral certificate adds complete lower layers for every part size in the theorem; divisibility is not needed. This is neither general spectral Chvátal H nor a perfect-matching factorization result.

## 1. Statement and conventions

Let $r\ge3$ and $s\ge100r^3$ be integers, with $r$ parts of size $s$. A transversal uses at most one point from each part; a top tuple uses all parts. Turn a family of faces into a graph by joining distinct disjoint faces. Give each top tuple $U$ weight $p_U\in[0,1]$, allowing zero and omitting such vertices at the end. Suppose the missing weight
$$
d=\sum_{U\ni i}(1-p_U)
$$
is the same for every point $i$, and $d\le K/(10r^3)$, where $K=s^{r-1}$. The weighted theta of this top graph is $M=K-d$. With all proper nonempty transversals added at unit weight, the weighted theta is exactly
$$
M+\sum_{k=1}^{r-1}\binom{r-1}{k-1}s^{k-1},
$$
with no divisibility restriction. There are no newcomer markers, no arbitrary lower-layer weights, and no guarantee for unbalanced top deficits. With binary weights, inclusion of all proper faces makes the positive support an actual downset after adjoining the empty face. Some proper faces may be in no allowed top tuple. The empty set is a universal vertex of weight 1 if included, which cannot increase an existing theta at least 1 (join rule).

Use the weighted primal with trace-one PSD matrix $Y$, entries zero across disjoint pairs, and objective $\langle\sqrt p\sqrt p^\top,Y\rangle$. A sufficient certificate is a symmetric $H$ with diagonal $1/p_U$, zeros on distinct intersecting pairs, and $H-J/M\succeq0$. Congruence by $\operatorname{diag}(\sqrt p)$ and testing against $Y$ proves the upper bound. A star of weight $M$ gives the opposite bound by its rank-one feasible primal. For separate layers with such $H_k,H_k-J/M_k\succeq0$, the block sum $H$ works for $M=\sum M_k$: $\sum (\mathbf 1^\top y_k)^2/M_k\ge(\sum \mathbf1^\top y_k)^2/\sum M_k$. The same center attains all layer masses here. This justifies subadditivity without an equality assumption across ranks.

## 2. Independent audit of the weighted top argument

**Product seed.** Let $X$ be tuple by point incidence and $A=(J_s-I_s)^{\otimes r}$. Put $q=(s-1)^{r-1}$, $P_*=I+A/q-J/K$, and $L=I-P_*$. Tensor factors decompose into constants and their orthogonal complements. The sector with $j$ mean-zero factors has $A$-eigenvalue $(-1)^j(s-1)^{r-j}$. The $j=0,1$ sectors exactly equal $\mathrm{im}\,X$ and are annihilated by $P_*$. The exact complementary minimum of $P_*$ is $1-(s-1)^{-2}$, from $j=3$ (the note's $\mu=1-(s-1)^{-1}$ is a safe lower bound, not exact). The smallest eigenvalue of $L$ is $-1/(s-1)$, from $j=2$; also $\lvert L_{UV}\rvert\le1/q$ directly. Bernoulli gives $q/K=(1-1/s)^{r-1}\ge1-(r-1)/s>.99$.

Write $F=\operatorname{diag}(\sqrt p)$, $f=F\mathbf1$, $\Delta=I-F^2$, $\widehat X=FX$, and $R=FP_*F+\Delta=I-FLF$. At the diagonal and distinct intersecting pairs, $R$ has exactly the entries of $I-ff^\top/K$. For $y\in\ker\widehat X^\top$, $Fy\perp\mathrm{im}\,X$. Hence $y^\top R y\ge\mu\|Fy\|^2+\|\Delta^{1/2}y\|^2\ge\mu\|y\|^2$, even when $F$ is singular.

**Persistent incidence kernel.** Let $Z$ be the part-constant vectors whose $r$ constants sum to zero. The complete Gram $G_f=X^\top X$ has diagonal part blocks $KI_s$, off-diagonal blocks $(K/s)J_s$. On the partwise mean-zero subspace it is $K I$, on the global constants $rK I$, and on $Z$ zero. Thus
$$
G_f^+=K^{-1}(I-\operatorname{blockdiag}(J_s/s))+(r^2sK)^{-1}J,\qquad \|G_f^+\|_\infty\le 3/K.
$$
The loss $B=X^\top\Delta X$ is symmetric, entrywise nonnegative, has row sums exactly $rd$, and kills $Z$. Hence $G=\widehat X^\top\widehat X=G_f-B$ has no additional kernel and has eigenvalues at least $K-rd>0$ on $Z^\perp$. The Neumann series $(I-G_f^+B)^{-1}G_f^+$ acts as the inverse on $Z^\perp$, zero on $Z$, and has absolute row-sum norm at most $(3/K)/(1-3rd/K)\le4/K$. In fact $G$ on the global constants is $rM I$, and leaves the partwise mean-zero space invariant, also giving a direct block derivation. The pseudoinverse statement does not require full column rank.

**Projection bounds.** Set $\Pi=\widehat XG^+\widehat X^\top$, $E=I-\Pi$. From $LX=X$,
$$
R\widehat X=FL\Delta X,\qquad \widehat X^\top R\widehat X=B-(\Delta X)^\top L(\Delta X).
$$
The first matrix has each entry at most $d/q$ in absolute value. Each row of $\widehat X$ has absolute sum at most $r$. Therefore $\|\Pi R\|_{\max}\le4rd/(Kq)\le5rd/K^2$. The second matrix is PSD and at most $rd(1+1/(s-1)) I$, since $(\Delta X)^\top(\Delta X)\preceq B$ and $\lambda_{\min}(L)=-1/(s-1)$. Rows of $\widehat XG^+$ have Euclidean norm at most $\sqrt r/(K-rd)$. Consequently $\|\Pi R\Pi\|_{\max}\le2r^2d/K^2$ at the stated cutoff. These estimates use no sign or sparsity assumption on $G^+$.

**Compatible correction.** Let $\mathcal M$ keep entries at nonempty intersection, including the diagonal, and zero all others; restrict to real symmetric matrices. For a masked $Q$ call it compatible if $(Q\widehat X)_{U,i}=0$ whenever $i\in U$. Take
$$
Q_0=\mathcal M(I-ff^\top/M-ERE).
$$
For $i\in U$, masking does not alter multiplication by column $\widehat X_i$ in row $U$. The term $ERE$ annihilates $\widehat X$, and $f^\top\widehat X_i=M$ by balance, so $Q_0$ is compatible. The signed offset $1/K-1/M$ costs at most $2d/K^2$ in max-entry norm. Hence $\|Q_0\|_{\max}\le(10r+2+2r^2)d/K^2\le6r^2d/K^2$.

For an outside point $i\notin V$, the entry $(\widehat X^\top Q)_{i,V}$ has contributions from at most $(r-1)K/s$ tuples: choose a point of $V$ not in $i$'s part and a tuple through these two. For inside points it is zero by compatibility. Thus $\|\Pi Q\|_{\max}\le4r^2\|Q\|_{\max}/s$. The absolute column sums of $\Pi$ are at most $4r$ by $\|\widehat X\|_1\le K$, $\|G^+\|_1\le4/K$, and $\|\widehat X^\top\|_1\le r$. So $T(Q)=\mathcal M(EQE)$ preserves compatible symmetry, and
$$
\|T-I\|_{\max\to\max}\le(2+4r)4r^2/s\le20r^3/s\le1/5.
$$
Preservation is exact: $EQE\widehat X=0$ before masking, and masking changes none of these incidence sums at the retained positions. Then $Q=T^{-1}Q_0$ has max-entry norm at most $\|Q_0\|_{\max}/(1-1/5)$. Sparsity belongs to $Q$, not $EQE$: every row of $Q$ has at most $rK$ entries. By symmetry $\|EQE\|_2\le\|Q\|_2\le rK\|Q\|_{\max}\le3/4<\mu$. Thus $P=ERE+EQE$ is PSD, annihilates $\widehat X$, and locks exactly to $I-ff^\top/M$. Its positive-support principal block, congruenced by $F^{-1}$ and augmented by $J/M$, is the required physical $H$. It also sends each weight-vector star $p\cdot X_i$ to $\mathbf1$. At $p_U=0$, projection rows are zero and $ERE$ has isolated unit diagonal; the compatible space with the correction rows zero is preserved by $T$. No zero reciprocal is taken at any stage. At $d=0$ correction vanishes; all inequalities are non-strict in the allowed upper boundary of $d$, with strict PSD margin on the orthogonal complement. The exact numerics below are controls, not the positivity proof.

## 3. Proper transversal layers: a divisibility-free weighted adjacency

The proposed clique works even under the weaker divisibility $k\mid rs$: write the parts in a cyclic list repeated $s$ times and break it into chunks of length $k\le r$, assigning each part's $s$ points once. If $k\nmid rs$ there cannot be a spanning matching of $k$-faces, so the fractional clique cover cannot prove this equality. The following alternative needs only $1\le k<r$ and $s\ge\max(2,k)$.

For a $k$-face $U$ denote its part support by $S$, and similarly use $T$ for $V$. Let
$$
 b=\binom{r-1}{k-1},\quad c=\binom{r-2}{k-2}\ (c=0\text{ if }k=1),\quad
 \epsilon=\frac{b}{s(b-c)+kc},\quad \gamma=\epsilon(b-c)>0 .
$$
Set $W_{S,T}=1-\epsilon |S\cap T|$. It is entrywise nonnegative because $s\ge k$. Define a symmetric weighted adjacency $A'$ on disjoint faces by
$$
 A'_{U,V}=\frac{W_{S,T}}{(s-1)^{|S\cap T|}s^{\,k-|S\cap T|}}
 \quad\text{if }U\cap V=\varnothing,\qquad A'_{U,V}=0\text{ otherwise}.
$$
For a row $U$, a specified point outside $U$ lies in a part either absent or present in $S$. In the first case the sum of the entries of $A'$ for all $V$ containing it is $(b-\epsilon kc)/s=\gamma$. In the second case the sum is $(b-\epsilon[b+(k-1)c])/(s-1)=\gamma$. Indeed, summing over disjoint label assignments for each $T$ removes its displayed denominator except the one forced point, and $\sum_{T\ni i}|S\cap T|=kc$ or $b+(k-1)c$ respectively. For a point already in $U$ the sum is zero. Thus $I+A'/\gamma$ sends every star to $\mathbf 1$.

It remains to prove this matrix PSD. Decompose functions on the labels of any support $S$ into products of constants and mean-zero label modes. For a fixed set $J_0$ of $j$ mean-zero parts, the only retained support indices contain $J_0$. A mode on any missing input/output part averages to zero; across shared parts the disjoint-label averaging matrix $(J_s-I_s)/(s-1)$ acts by $-1/(s-1)$ on mean-zero vectors and by 1 on constants. Thus on this sector the matrix is $(-1)^j/(s-1)^j$ times the principal support submatrix of $W$. Normalized constant label factors introduce no additional ratio because all supports have size $k$. At $j=0$, $W=J-\epsilon B_0B_0^\top$ with $B_0$ the $k$-part-support incidence matrix and $B_0^\top B_0=(b-c)I+cJ$. The part contrasts have eigenvalue $-\gamma$; the constant eigenvalue is a nonnegative row sum, and the remaining eigenvalues are zero. At $j\ge1$, the nonnegative restricted $W$ has row sum at most its row sum restricted only to one fixed part in $J_0$, namely $(s-1)\gamma$. All eigenvalues of the sector are therefore at least $-\gamma/(s-1)^{j-1}\ge-\gamma$. This proves $H_k=I+A'/\gamma\succeq0$. As $H_k x_i=\mathbf1$ for a star $x_i$ of mass $b s^{k-1}$, Gram Cauchy-Schwarz gives $H_k-J/(b s^{k-1})\succeq0$. This treats singletons as well ($H_1=J$).

For example, at $r=3,s=3,k=2$, this gives exact theta 6; the optimal fractional clique-cover cost is $27/4$. No clique can exceed 4 because there are 9 points, and a four-edge matching exists by grouping a cyclic part list into four pairs. Transitivity then averages the maximum clique to optimal fractional cover; the cardinality lower bound holds for any cover. The gap also occurs inside the top theorem, e.g. $r=3,s=2701$. More generally any complete $k$-layer with $k\nmid rs$ has fractional clique-cover lower bound $\binom rk s^k/\lfloor rs/k\rfloor>b s^{k-1}$.

## 4. Reproducibility and primary scope

Run `node transversal-replay.cjs` and compare the [expected output](transversal-expected.json). The [replay](transversal-replay.cjs) requires no dependencies. Core matrix replay uses rational Gauss-Jordan on $G+P_Z$, subtracting $P_Z$ afterwards, independently of the proof's Neumann argument. The $r=4,s=3$ weighted fixture has three pairwise distinct perfect matching classes with weight deficits $1,3/4,7/16$, exactly $d=35/16,M=397/16$, and 3 omitted rows. It verifies 144 pseudoinverse and 144 kernel entries, 6,561 projection idempotence entries, the 972 first and 144 second identities, and 324 initial and iterated compatibility positions. All 312 positive incidence positions fail if $K$ is used instead of $M$ in $Q_0$. Its three-factor tensor realizes exact seed gap $3/4$, not the conservative $1/2$. This fixture is below the theorem cutoff and does not claim a convergent correction there. Lower certificate replay uses direct exact rational semidefinite Schur elimination of the full 27-by-27 and 90-by-90 matrices, and separately all tensor support-sector ranks plus direct star sums for four layers (including nondivisible rank three).

A theorem-domain control: $r=3,s=2701$; the $\mathbb Z/s$ diagonal translation orbits of tuples are disjoint perfect matchings. Delete every orbit meeting a fixed two-part prefix $(0,0)$. There are exactly $s$ such orbits (one for each last label), so $d=s$ at every point and the chosen admissible prefix has **zero** extensions. This obeys the budget and gives top theta $7{,}292{,}700$, full theta $7{,}298{,}103$, with $2\nmid 3s$. This control deliberately has an elementary matching cover for the top; it demonstrates failure of high-codegree implications, not failure of all dense-product matching decompositions. Conversely mere balance cannot imply theta equality: the four even-parity triples in $[2]^3$ have all degrees 2 but no two disjoint triples, giving theta 4. Scaling weights by $1/2$ gives a fractional perfect matching in the hypergraph, but no convex combination of perfect matchings.

Closest checked primary work: the complete, unweighted signed-set **independence number** is classical (Bollobás–Leader, Deza–Frankl; see [Borg's theorem and attributions](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v14i1r41) and [Borg's injection paper](https://arxiv.org/abs/1912.10324)). It does not establish either the heterogeneous theta or the perturbation. Balanced multipartite graph clique decompositions [arXiv:2604.25206v2](https://arxiv.org/abs/2604.25206) decompose graph edges into small cliques, not transversal hyperedges into spanning matchings. [Cheng, Glock and Lane, arXiv:2611.27781v1](https://arxiv.org/html/2611.27781v1#S10) prove integral nonpartite 1-factorization at high codimension-one degree; their partite high-codegree analogue is itself Conjecture 49. Our theorem allows partite codegree zero and fractional tuple weights; it proves an SDP bound rather than any factorization. A matching cover of weight $M$ *would* imply the top conclusion, but in the balanced case necessarily decomposes $p/M$ into perfect matchings (every matching has at most $s$ tuples and the total tested tuple share is $sM$). Ordinary fractional vertex matching only satisfies point marginals, and does not do that. No exact identical prior emerged from the primary search; this is not an exhaustive novelty determination.

[Ellis–Filmus–Friedgut §5](https://arxiv.org/html/2609.28404v3#S5) formulate spectral H generally. The [earlier Research1 rank-three certificate](heterogeneous-weights.md) uses a near-complete unrestricted top layer; the previous higher-rank extension still uses that seed. Here product-transversal density relative to all $r$-sets tends to $r!/r^r$; for $r=3$ that is $2/9$, so nearly 78% of the unrestricted triples are already absent before any $d$-loss. The present theorem instead needs its structured product and regular near-density within that product. It is a valid scoped companion rather than a global advance on H or a replacement of the stronger near-complete rank-three constants.
