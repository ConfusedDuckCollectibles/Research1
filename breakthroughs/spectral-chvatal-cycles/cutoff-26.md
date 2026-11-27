# Triple-type refinement of the cycle PSD bound to n≥26

## Result

Use the notation of the published [Research1 cycle proof](https://github.com/ConfusedDuckCollectibles/Research1/tree/main/breakthroughs/spectral-chvatal-cycles), also available as [proof.md](proof.md). Its masked inverse, coupled Schur argument, and tail n≥32 stay the same. A refined Frobenius bound certifies every eligible n=26,…,31 and additional ranges at n=22–25. It follows that every simple 2-regular anchor graph with r≥5 is covered for n≥max(26,2r+1). In particular r≥13 is covered at every admissible n. The exact extra eligible pairs (n,r) are (22,5), (23,5–6), (24,5–8), and (25,5–10). This leaves 141 of the earlier 468 partition-and-size cases before subtracting independent small proofs.

## Refined bounds for an ordered pair of locked triple entries

Write K=binom(n−1,2), a=binom(n−2,2), b=binom(n−3,2), q=binom(n−4,2), L=K−2r, and ω=6/(n−1), δ=1/q−1/K. Let U have u anchors and t tips. Put h_U=2u+t and ρ_U=(1−ω)h_U+rω. These depend only on the type of U. Use c=(n−2)/(Ka)+12/[a(b−12)] and v=6(n−2)/(b−12) from the original proof.

The contribution of the coupled rank block B−J/K is bounded entrywise by the geometric mean of the following nonnegative row quantities. This follows either from its explicit Gram sum, or by Cauchy–Schwarz applied across all its rank-block terms:
    b_U = 2[u/K²+(r−u)c²] + (4/L)[u/K+(r−u)c]²
          + [t+(r−t)v²]/[K(K−1)].
This uses the exact counts u and t of incident anchor and tip columns.

For the terms from the projection of the Kneser survivor block, let X_D be the deleted-triple incidence and let x_U denote the incidence row of U. Put G=XᵀX, G_f=aI+(n−2)J, and u_U=x_U G^−1 X_Dᵀ. Its full-Gram approximation has components (ell_d−ω)/a where ell_d=|U∩d| for each deleted triple d. Since U survives, 0≤ell_d≤2. Since each anchor lies in two deletions and each tip lies in one, sum ell_d=h_U. For 0<ω<1, at least h_U/2 nonzero overlaps give
    sum_d |ell_d−ω| ≤ ρ_U,
    sum_d (ell_d−ω)² ≤ (2−2ω)h_U+rω².
The latter bound uses ell_d²≤2 ell_d.

By the inverse identity G^−1−G_f^−1=G_f^−1 X_DᵀX_D G^−1, the row difference x_U(G^−1−G_f^−1) has l1 norm at most 6ρ_U/[a(a−12)]. Indeed each row of X_D has l1 norm 3 and G^−1 has maximum absolute row sum at most 2/(a−12). Multiplication by X_Dᵀ has l1-to-l2 norm at most sqrt(2), since each coordinate belongs to at most two deleted triples. So the l2 error in u_U is at most e_U=9ρ_U/[a(a−12)]. Young's inequality, with ε=1/12, now gives
    ||u_U||² ≤ z_U = (13/12)[(2−2ω)h_U+rω²]/a² + 13 e_U².
The quadratic projection term is bounded by (1+r/q)sqrt(z_U z_W), using the old deleted principal block's row-sum norm bound. Combining it and the rank block by Cauchy–Schwarz gives a geometric mean with row quantities A_U=b_U+(1+r/q)z_U. For exact rational arithmetic bound this by (A_U+A_W)/2.

For a first-order projection between old survivors use the mixed block identity ΠP_s=−XG^−1 X_Dᵀ P*_(D,S). Its exact cross entries in the last factor are δ−1_(d meets W)/q. The constant term for the full-Gram approximation has absolute value δ|h_U−rω|/a. In the signed sum over deletions meeting W, the positive part is at most (1−ω/2)h_U, and the negative part at most ω h_W. Indeed ell_d−ω≤(1−ω/2)ell_d at positive overlaps, and there are at most h_W deletions meeting W. Bound the absolute signed sum by the larger of these two bounds, divided by aq. For the perturbation part use |(X_DᵀP*_(D,W))_i|≤2/q and the l1 row-error bound above. Thus
    c_U = δ|h_U−rω|/a +12ρ_U/[a(a−12)q],
    |(ΠP_s)_(U,W)| ≤ c_U + max((1−ω/2)h_U, ωh_W)/(aq).

The total locked error at U,W is bounded by
    F(U,W) = (A_U+A_W)/2 + c_U+c_W
             + [max((1−ω/2)h_U,ωh_W)+max((1−ω/2)h_W,ωh_U)]/(aq).
This follows by adding the two first-order terms, the quadratic term and the rank block on the mask. The bound depends only on the two triple types.

## Finite exact calculation

A type x=(u,t,l) has count N_x=binom(r,u)binom(r,t)binom(n−2r,l) in the complete old triple layer. For each type x and y=(u',t',l'), the number of ordered meeting pairs including the diagonal is exactly
    N_x [N_y−binom(r−u,u')binom(r−t,t')binom(n−2r−l,l')].
This counts the complete layer and therefore also bounds the corresponding count after deletions. Sum this count times F(x,y)². The resulting exact rational number bounds the squared Frobenius norm of the locked error. Compare it to μ²κ² from the published proof, with μ=1−1/q and κ=1−6(n−4)/(a−6). By that proof, strict inequality proves PSD of the complete coupled Gram and hence weighted theta equality. The extension to all two-layer ratios follows from the same crossing and endpoints.

Run `node cutoff-replay.cjs` from this directory. It checks every eligible (n,r) between n=22 and 31 and counts all cycle partitions. All eligible pairs at n=26,…,31 pass. The most crowded n=26 case r=12 has squared ratio about 0.99415918, with the comparison done as an exact integer inequality. The calculation closes 327 of the 468 previously listed partition-and-size cases. This refinement and the n≥32 tail use upper bounds, so failure of the sufficient inequality leaves a smaller case open. No claim is made for full H or for the independent square proof. [Ellis, Filmus and Friedgut, §5](https://arxiv.org/html/2609.28404v3#S5)