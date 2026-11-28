# Full theorem: minimum-degree normalization for growing triple deletions

**Theorem.** The equality holds in the claimed sufficient range $n\ge100,\ D\le n/8$, including patterns that touch every coordinate. The pivotal change is to normalize the fitted rank part at the largest surviving triple-star size $M$, while keeping the complete Kneser baseline $K$ in the projected comparison. Their explicit constant offset is necessary. No inference is made about multiple markers, arbitrary non-layer weights, general Conjecture H, or the optimal cutoff.

## Statement and residual certificate

On $n$ old coordinates the inner layer contains all nonempty old one- and two-faces and all old triples not in a simple deletion hypergraph. The outer layer restores the deleted triples, adds the marker singleton, and adds a bridge pair to the same marker at every old anchor in $A$. Give inner faces weight $u+v$, and outer-only faces $v$, with $u,v\ge0$. Write $d_i$ for deletion degrees, $d_0=\min_i d_i$, $D=\max_i d_i$, $K=\binom{n-1}{2}$, $M=K-d_0$, $r_i=d_i-d_0$, $t_i=M-r_i$, and $S_0=n+K$. For $A\ne\varnothing$, put $b=\min_{i\in A}r_i$; the exact theta and maximum star value is
$$
 S_0(u+v)-d_0u+\max(0,v-bu).
$$
For $A=\varnothing$, omit the correction. The old star at $i$ has value $S_0(u+v)-d_i u+v\,1_{i\in A}$. A minimum-degree anchor, when present, makes $b=0$ and attains the correction for all ratios. For $b>0$, a minimum-degree site is unbridged and a $b$-anchor attains the switched branch. Marker score $(|A|+1)v$ is never larger.

In the residual use surviving triples of weight one and keep only selected positive-surplus bridges of weight $w_i=r_i$. Let $X$ be the surviving triple incidence matrix, $x_i$ its columns, $\alpha=\sum w_i$, $L=M-\alpha$, and $\lambda_i=1$ at a retained bridge and $t_i/M$ otherwise. For the avoiding-coordinate fit $C_i$, keep starting value $1/a$, $a=\binom{n-2}{2}$. The adjustment at $j\ne i$ is a difference of the deleted avoiding-$i$ degree divided by $a$ and $1-\lambda_j$; both lie in $[0,D/a]$ since $M\ge a$ follows from $n-2-d_0\ge0$. Each unbridged $i$ with $r_i>0$ uses the unchanged fit $X^\top V_i=t_i e_i$, putting $V_i=1$ on incident triples. Zero bridge weights and zero surplus terms are omitted, not inverted.

Form $g=\mathbf1-Cw$, $F=\sum r_i V_iV_i^\top/(M t_i)$ over positive-surplus unbridged sites, and $B=C\,\operatorname{diag}(w)C^\top+gg^\top/L+F$. At a positive bridge $t_i+r_i=M$, otherwise $t_i=M\lambda_i$, so both give $g^\top x_i=L\lambda_i$ and
$$
 Bx_i=\begin{cases}
 \mathbf1-r_i C_i&\text{positive-weight bridge},\\
 \lambda_i\mathbf1+(r_i/M)V_i&\text{otherwise}.
 \end{cases}
$$
The $V_i$ term is zero at $r_i=0$. Thus every locked incidence $(Bx_i)_U$ with $i\in U$ equals one, not just incidences on tight stars. With corrected old block $B+P$, bridge diagonal $\operatorname{diag}(1/w)$, and old/bridge block $C$, the Schur complement is $P+gg^\top/L+F$. A minimum-degree coordinate has no positive bridge, surviving star weight $M$, and its incidence vector maps to $\mathbf1$ in the normalized Gram. Cauchy–Schwarz gives the weighted theta certificate at $M$, attained at that star. Reducing or dropping bridge weights preserves the lower bound, hence equality also holds for any individual bridge weights in $[0,r_i]$. In particular this proves the all-touched pure **inner** bound with no bridges. Merely restricting the complete graph would give $K$, not $M$.

## Exact offset and quantitative check

With $W=C-J/M$,
$$
 B-J/M=W\operatorname{diag}(w)W^\top+(Ww)(Ww)^\top/L+F.
$$
Use $M$ in the rank estimate, both $M,t_i$ in $F$, and bound surplus by $D$ pointwise and $3m$ in total. The bound off incidence is $|W_i(U)|\le(M-a)/(Ma)+6D/[a(b_{\rm bin}-6D)]$, $b_{\rm bin}=\binom{n-3}{2}$, at most $10/n^3$ using $0\le M-a\le n-2$, $M,t_i\ge11n^2/25$, $a\ge9n^2/20$, $b_{\rm bin}-6D\ge2n^2/5$. With $L\ge3n^2/10$, the exact rank coefficient bounds are $1875/484$ at $n^{-3}$ and $2003825/11616$ at $n^{-4}$; hence $R\le4/n^3+185/n^4$.

Keep $P_*=I+A_{\rm disj}/q-J/K$, $q=\binom{n-4}{2}$, because this annihilates **full** incidence. The survivor block $P_s$, compression $P_0$ off survivor incidence, masked inverse and first- and second-order estimates remain unchanged: $\mu=1-1/q$, $\kappa=1-6(n-4)/(a-3D)$. If $\mathcal M$ masks intersecting positions including diagonal, the actual defect is
$$
 Q_0=\mathcal M\big[-(B-J/M)+(P_s-P_0)-(1/M-1/K)J\big].
$$
The minus sign is exact since $P_s=I-J/K$ on the mask. Omitting that last term adds the spurious incident row sum $(1/M-1/K)t_i$, nonzero if $d_0>0$. Replacing $J/K$ by $J/M$ in the comparison would give the **same compressed matrix** $P_0$ because the projection off incidence kills $J$; it would not give the same before-compression nullspace identity used to bound the mixed deleted-row terms.

The offset magnitude $d_0/(KM)$ is at most $(125/198)n^{-3}$. The rank bound, unchanged projection terms, and at most $n^5/4$ ordered mask entries give $\|Q_0\|_F\le(14+260/n+900/n^2)/(2\sqrt n)$. For $n\ge100$, this is at most $1669/2000<76483/90000\le\mu\kappa$. The terminal rational gap is $689/45000$. Positive denominators and strict coercivity hold also at exact $D=n/8$ when integral. No claim follows about prior low-$D$ cutoffs below 100 without rechecking the changed normalization.

## Lift and edge cases

Pure inner at $(1,0)$ has value $S_0-d_0$ by the no-bridge residual and full lower ranks. Pure outer costs $K$ for full triples and $n+1$ for the remaining faces: color old pairs $i+j$ and bridges $2i$ modulo odd $N\in\{n,n+1\}$, then assign singletons missing colors in an $n+1$ palette. Every color is a disjoint class. With nonempty $A$ this is attained at an anchor. With empty $A$, the marker is disjoint from every old face and its theta contributes only a maximum with the old part, so the value is $S_0$.

For $b>0$, at $(u,v)=(1,b)$, add the residual with bridges reduced to $b$, $b$ full triple layers, and $b+1$ full lower-rank layers. This gives $s=(b+1)S_0-d_0$ before the marker. In the unscaled dual Gram the bridge block is $tI-J$, $t=s/b$. The marker vector needs inner product $-1$ with each bridge and squared norm $t-1$. Its prescribed projection has squared norm $|A|/(t-|A|)\le t-1$, since $s-(|A|+1)b=S_0-d_0+b(S_0-|A|-1)\ge0$. Other old faces are disjoint from the marker and impose no fixed dual entry. At the crossing a minimum-degree star attains $s$. For $v\le bu$, combine crossing coefficient $v/b$ and pure inner coefficient $u-v/b$; for $v\ge bu$, use crossing coefficient $u$ and pure outer coefficient $v-bu$. Weighted theta is homogeneous and subadditive by its common handle formulation. If $b=0$, a minimum-degree anchor maximizes both pure endpoints and the crossing is unnecessary. The all/one/proper/empty anchor cases and $(u,v)=(0,0)$ do not require reciprocal zero weights.

## Independent exact replay and scope

[All-touched replay](all-touched-replay.cjs) and [expected output](all-touched-expected.json) contain the source and transcript. Run `node all-touched-replay.cjs`. BigInt comparisons are exact; displayed squared ratios are approximations. Modular fits crosscheck identities and do not supply real PSD.

* **Boundary orders, all touched:** at $n=100$, delete all cyclic translates of $\{0,1,3\},\{0,1,4\},\{0,2,6\},\{0,2,7\}$ except the translation-zero triple of each of the first three types. The 397 distinct deletions have degree histogram $9:1,\ 10:1,\ 11:4,\ 12:94$. The 26 positive bridge anchors include three distinct degrees and $b=1$; 73 positive-surplus coordinates are unbridged. Fits from the actual 161,303 surviving triples pass 30,000 coupled coordinates including 900 locks. Exact squared sufficient ratio is about $0.2925327007$. Six anchor patterns (none, single min, single max, proper subset excluding/including a minimum, all) check 36 score cases and all corresponding color assignments.
* Before restoring the three triples the cyclic family is regular: 400 deletions, $d_0=D=12$, no positive-surplus bridges despite all anchors being chosen. The ratio is about $0.2689326330$; this exercises empty auxiliary columns at the maximal integer minimum degree allowed at $n=100$.
* At $n=10,d_0=5$ the direct projected-mask test gives 1,010 null incidence entries, 303 correct mask compatibilities, and **303 failures** in the counterfactual defect omitting the offset. This is a sign/necessity control, not a PSD result below the theorem range.

Compare the [original growing-deletion construction](proof.md) and [Ellis–Filmus–Friedgut §5](https://arxiv.org/html/2609.28404v3#S5). The underlying masked inverse is proved in [the original construction](proof.md). This removes an actual hypothesis in that construction, not the remaining general limitations.
