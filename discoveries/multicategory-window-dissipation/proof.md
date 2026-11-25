# A large same-data separation for finite-window multicategory dissipation inference
November 25, 2026. AI-assisted Research1 note on an exact constructed stationary continuous-time Markov example and family. No experiment or peer review. This is **not** a universal improvement over inference with the original uncensored path; the stronger existing excursion-rate comparator is quantified below.

**Companion:** For genuine first exits, the [completed-window refinement](first_exit_tent_audit.md) doubles the active-pair lower bound on the same table. This file gives the general classifier result, the model, and its different robustness guarantee.

## Claim and data contract

For an ordinary all-even stationary finite reciprocal CTMC with microscopic time-reversal entropy production rate sigma (natural logs, k_B=1), record two equally long label paths of half-duration t, beginning at the observed midpoint and read outward. On midpoint label H, independently apply the same classifier in both directions: return the first visible endpoint label j reached by time t, or 0 if none is reached. On all other midpoint labels return the common sentinel on both halves. The measured data are the midpoint probability p_H, half-duration t, and the **conditional finite-window joint categorical law M given midpoint H**. No observations of microscopic states, actual exit times, later jumps, or completed excursions crossing a window boundary are needed. There is no filtering on which halves succeed.

We compare Petkov's optimized finite-window CP deletion floor, allowing it the same observed midpoint stratification, to a pair-resolved convex affinity objective on the same Gram constraints. The two may be combined with a maximum for general applications; no universal dominance of the raw pair objective is needed. We construct M from a stationary irreducible reciprocal 15-state CTMC with six instantaneous labels. Its entire visible path is time reversible, and even every matched binary recoding of the same finite categorical table is blind. Nonetheless the new inferred lower bound is at least 2.27589 times the optimized Eq. (11) floor at a fixed witness and at least this comparison diverges logarithmically along the specified family.

## Affinity bound from the same finite table

At hidden midpoint x, let u and v be the conditional forward and reversed-path classifier probabilities. Conditional independence gives M=E_H[uv^T]. The usual stationary forward/reverse path KL in each direction, averaged over all midpoint microstates, is t sigma. The two-sided average of KL after classification therefore gives t sigma >= p_H E_H[J(u,v)/2], where J is the categorical Jeffreys divergence. Restricting to H drops nonnegative contributions only. This can be used even if holding times and total rates are unknown to the observer. At the physical level an entropy/heat interpretation supposes local detailed balance.

Put w_i=sqrt(u_i v_i) and form the completely positive Gram matrix G=E_H[ww^T] of [Petkov Eq. 13](https://hexagonmath.org/pdf/2611.00658v1). It satisfies G_ii=M_ii and G_ij<=sqrt(M_ij M_ji) by Cauchy-Schwarz. For b=u_i v_j, c=u_j v_i, the contribution of unordered pair i,j to J(u,v)/2 is (b-c)log(b/c)/2. This follows by writing D(uv^T||vu^T)=J(u,v) and grouping pairs. Write

    phi(x) = sqrt(1-x^2) atanh(sqrt(1-x^2)), 0<x<=1;
    (b-c)log(b/c)/2 = (b+c) phi(2 sqrt(bc)/(b+c)).

The extension is zero at x=1 and +infinity at 0 unless the pair weight is zero. The decreasing, convex scalar function and pair conditioning argument already appear for excursion measures in Petkov's work (see prior art below). Weighted Jensen separately on each unordered pair yields the finite-window certificate

    t sigma >= p_H inf_G SUM_{i<j} (M_ij+M_ji) phi(2G_ij/(M_ij+M_ji)).       (1)

The infimum obeys the specified CP, diagonal and entry caps. Any term of weight zero is omitted. No assertion that an optimizing CP G is realizable for an arbitrary observed table. A lower bound on this infimum is sufficient. The same midpoint-stratified [Petkov Eq. 11](https://hexagonmath.org/pdf/2611.00658v1) is p_H·2F(Delta_CP(M)), where F(d)=sqrt(d) atanh(sqrt(d)) and Delta_CP is one minus the mass of the largest entrywise CP minorant. The unconditional table's two observed strata give a weaker global Eq. 11 by convexity; the comparison below grants the stronger stratified version.

## Exact same-data separation

Take 0<e<=1/16, a=1/4, z=1/16, and use a five-vertex path 1--2--3--4--5. Define the 5x5 matrix A of diagonals (1,2,2,2,1), four path entries 1 in each direction, rare closing edge 1--5 with entries a, and all other entries zero. Let D=diag(e,1,e,1,e), P= D A D / Z, where Z=4+8e+(4+2a)e^2. P is normalized, and alpha=2a e^2/Z is its total rare-edge probability. Write b=P 1 and p=1-z. The exact six-category law, including first-exit censor, is

    M = [[z^2, zp b^T],
         [zp b, p^2 P ]].                                                    (2)

Positivity: y^T A y = Σ_i(y_i+y_{i+1})^2 + 2a y_1 y_5. If s_i=y_i+y_{i+1}, then y_1-y_5=s_1-s_2+s_3-s_4, so (y_1-y_5)^2<=4 Σ s_i^2. Therefore y^T A y >= (1-4a)Σ s_i^2 + a(y_1^2+y_5^2)>=0. Thus P is PSD. P-bb^T is PSD by Cauchy-Schwarz in the P inner product and 1^T P1=1, which proves PSD of (2). Every common deterministic or independent stochastic binary recoding is a congruence of this matrix, hence PSD with nonnegative entries; its binary Gram witness detects nothing. All independent pair terms also vanish if the shared CP coupling is removed, even with their individual two-by-two PSD caps.

**New lower bound.** From the active 5x5 submatrix of any feasible G form Gtilde=(Z/p^2)D^-1 G_active D^-1. It is CP, has diagonal (1,2,2,2,1), and support contained in the triangle-free five-cycle. Every nonnegative Gram vector has support at most two vertices, since nonedges have exact zeros. Its trace is at least twice its single edge entry, by (x-y)^2>=0. Summing gives Σ_edges Gtilde_ij<=4. If g0 denotes rare-edge retained weight, and g1..g4 the four path weights (with caps 1), then Σ(1-g_i)>=g0. Drop the nonnegative censor-to-other terms from (1). Put u=g0/a; the active part is at least

    (2p^2/Z) * [e Σ_i phi(g_i) + a e^2 phi(u)]
    >= p^2 alpha * [2u/e - log u]
    >= p^2 alpha * [1 + log(2/e)].                                            (3)

Here phi(v)>=2(1-v) by convexity and its limiting derivative -2 at 1, and phi(v)>=-log v. For the latter, with x=sqrt(1-v^2), the difference x atanh x +(log(1-x^2))/2 has derivative atanh x>=0 and vanishes at zero. The last minimization is exact at u=e/2; zero endpoint is covered by a limit.

**Upper bound on the optimized old lower bound.** Delete only the rare pair of entries in P to obtain P0, which is CP: it is the sum of four scaled nonnegative rank-one edge vectors f_i f_i^T. Denote their total retained mass 1-alpha. Extend the Gram factor for each f_i to (z·1^T f_i, p f_i). Add a single censor-only factor with value z sqrt(alpha) in category 0. This yields a CP minorant G0 of (2) with censor diagonal z^2 and cross row zp·P0 1. Its missing mass is exactly d=(1-z^2)alpha; this correctly accounts for nonzero censor and retains most censor-incident matrix entries. Thus the actual optimized old t sigma floor (before p_H) is at most 2F(d), not the value of an arbitrarily chosen old copositive score. It is positive: a CP representation of M would contradict the active trace bound. As e decreases with z and t fixed, the ratio of (3) to 2F(d) grows as [p^2/(2(1-z^2))] log(1/e). The comparison is made with numerical upper and lower constraints in opposite directions and does not require determining either optimizer.

## CTMC that produces precisely (2), and its true EP

Set theta=a e/2. For each of the four path edges start with its raw rank-one contribution to D A D: endpoint diagonal entries D_ii^2,D_jj^2 and off-diagonal e. On the two outer edges subtract e^2 theta from the diagonal assigned to its endpoint 1 or 5. For rare edge 1--5 use endpoint diagonals e^2 theta and off-diagonal a e^2. These five edge tables R_e sum exactly to D A D. Each is symmetric with strictly positive entries on its pair; its determinant is zero on the two interior path edges and negative on the others. For a raw matrix [[d_i,m],[m,d_j]], split into the average of the two nonnegative rank-one matrices with the same diagonal and off-diagonal entries m±sqrt(m^2-d_i d_j) interchanged. This is the two-factor decomposition from [Petkov's related excursion paper](https://hexagonmath.org/pdf/2611.00657v1), §9. They have the same total mass T_e=d_i+d_j+2m and are transposes of each other. Divide each by T_e and take its row/column sums to obtain product factors u_h,v_h. All supported entries are strictly positive. Label each of the resulting ten microscopic states h by H. The five visible states j retain their distinct labels. Given midpoint H, assign h weight T_e/(2Z), equally for the transpose pair.

Choose p_H=1/2, visible stationary masses (1-p_H)/5, and exit rate lambda=1 from each h. On a supported h--j pair set the stationary directed masses pi_h q_hj=lambda pi_h u_h(j) and pi_j q_jh=lambda pi_h v_h(j), with no other off-diagonal transitions. Set diagonals to the negative exit rates. Row sums balance at each h, and transpose pairing balances inflow and outflow at every visible state. Support is reciprocal and the underlying connected cycle makes the 15-state chain irreducible. Swap the two microscopic states in every transpose pair, leaving visible states fixed: this preserves pi and conjugates Q to Q^*, as seen directly from the stationary directed masses. The permutation preserves all readout labels. Thus the entire stationary observed path law equals its reverse, not just each matrix or selected timestamps.

Forward and reversed holding time at each h is exponential with the same rate lambda. There are no hidden-to-hidden jumps, so the first jump already arrives at a visible endpoint. The first-exit classifier has conditional vectors (z, p u_h), (z,p v_h) exactly, with z=exp(-lambda t). Conditional independence gives (2) for t=4 log 2. No limit, no substitution of excursion counts for a finite-window law, and no ignoring secondary jumps that could change the first exit. All five nonedges of the active observed matrix are exact zeros because each microscopic hidden midpoint connects to only two visible labels joined by one cycle edge.

The steady EPR is the sum of pair flux differences times log flux ratios, viz. 2 lambda p_H B0, where the rank-one construction gives

    B0 = (4e/Z) phi(sqrt(1-theta)) + alpha phi(e/2).

Thus B0/alpha = 1+log(4/e)+o(1) as e→0. Even along the growing certificate-ratio sequence the lower bound (3), after multiplying p_H/t, certifies a fraction tending to p^2/(2 lambda t) of actual dissipation (~15.85%). This is a constant fraction, not only a ratio of tiny old and new numbers that both lose all of the true signal. The absolute magnitude of the rate in this chosen time unit does tend to zero; empirical acquisition cost increases. For a=0, omit the rare states, keep all other rates and setting theta=0; the remaining 13-state connected chain is detailed balanced and both witnesses are zero.

## Locked comparison, e=1/32, natural rate units

| Input and certificate | EPR rate |
|---|---:|
| Same finite data M, all binary coarsenings and all independent-pair determinant/affinity floors | 0 |
| Same M and any direct time-reversal comparison of past/future | 0 |
| Even the full visible path: direct reversal relative entropy and all ordinary time-antisymmetric correlation expectations | 0 |
| Same finite data, optimized midpoint-stratified Petkov Eq. 11, with both strata (conditional H plus sentinel stratum) | at most 0.00004123477217913 |
| Same finite data, new Gram-pair certificate | at least 0.00009384585085908 |
| Actual microscopic EPR | 0.0006717183575408 |
| Richer uninterrupted visible record with **uncensored excursion occurrence rates** (not the same finite inputs), Petkov Appendix J.2 affinity method with the same copositive constraint | greater than 0.00059209033436929 |

The finite improvement ratio is at least 2.27589; the new bound certifies at least 13.97% of actual EPR. Old bounds with the unstratified M and deterministic or independent same-channel motif coarsenings have no larger CP deficit (coarsening monotonicity, Petkov Proposition 4.2); more informative hidden-state or observation contracts are not claimed to be dominated.

The richer excursion comparator can indeed do substantially more. Given uninterrupted visible paths, the observed symmetric complete-excursion rate matrix for this chain is lambda p_H P (visible entry flux into h times exit probabilities from h). It provides new statistics: counts per total elapsed time and each completed visible-to-visible segment through H, with no finite boundary censor. Petkov J.2 bounds the contribution of each excursion pair with one directional rate R_ij by 4 R_ij phi(G_ij/R_ij); its excursion Gram is CP with diagonal no bigger than the observed diagonal. The same path-trace and scalar estimates in (3) with no censor imply at least 2 lambda p_H alpha[1+log(2/e)], i.e. about 88.15% of true EPR here. Therefore this is not an advantage when that richer record is available, and canonical path-irreversibility statistics should not be conflated with the much stronger constrained excursion method. The distinct result is a factor separation within the compressed finite-window contract.

## Positive tolerated perturbations, conservatively

Let K=p^2/Z. If a nearby symmetric conditional matrix changes by at most eta entrywise, rescaling its active block by D^-1 and 1/K changes entries by at most mu=eta/(K e^2). For any CP matrix, the five-cycle inequality 2Σ_edges Gtilde_ij <= tr(Gtilde)+2Σ_nonedges Gtilde_ij applies without exact support. To see this for a nonnegative rank-one vector with sum 1, maximize its sum of edge products: transferring weight between any nonadjacent occupied vertices to the higher linear objective reduces support until a clique remains, with best value at most 1/4. Summing the resulting normalized inequality over Gram factors gives the statement (this is the five-cycle Horn witness used by Petkov).

The scaled feasible trace is at most 8+5mu and the five nonedge entries are at most mu each (unordered). The combined retained edge sum is therefore at most 4+7.5mu. Actual scaled heavy motif entries have sum at least 4-4mu. Thus the heavy loss is at least g0-11.5mu. With actual rare scaled entry a0>=a-mu and u=g0/a0, repeating the proof of (3) gives a conservative new lower on the objective of

    2K e^2 (a-mu) [1+log(2/e) - 23mu/(e(a-mu))],

provided the final bracket is nonnegative. In the prescribed matched independent label-noise channel, the half-classifier categories are passed through the same stochastic matrix, and each row changes its original category with probability at most eta/2. Then joint table error is at most eta by a union bound, and the old optimized CP deficit cannot increase: conjugating a CP minorant by the same channel preserves all its mass and minorant validity. For mu<=a e/100 the fixed-witness margin is still at least 2.1737 at eta≈1.57614e-8; the e^3 tolerance preserves the asymptotic separation. This is a sufficient, deliberately pessimistic calibration condition, not a sample-size or real-device precision claim. Without structural zeros or a comparable error bound, a purely empirical table would need careful finite-sample uncertainty and a justified stationarity model. Finite consecutive snippets are not automatically independent. General arbitrary additive errors would also require a separate bound on the old CP statistic; the preserved *ratio* here is asserted only for matched independent coarsening noise.

## Primary prior art and scope

1. [Petkov, Complete positivity and equilibrium realizations of reversible observations, Hexagon v1](https://hexagonmath.org/pdf/2611.00658v1), Thm 4.1, Eqs. (11)–(16): the matched finite-window contract, conditional Gram and generic mass-only floor already exist. The same paper §5 and J.1–J.3 gives a physical five-cycle with reversible observed law; it also already uses the pair-affinity objective for **uncensored excursion rates** in J.2, aggregated under homogeneous cycle rates by Jensen. Its explicit finite-window proof instead compresses to a single global deletion mass. Thus no claim to the five-cycle architecture, Gram, copositive inequality, weighted affinity calculation alone, or the first dissipation inference from reversible records.
2. [Petkov, Entropy Production under Hidden-State Reduction, Hexagon v1](https://hexagonmath.org/pdf/2611.00657v1), Nov 19, §9: already gives the scalar and paired rank-one construction for two individually observed states with richer marked excursion statistics. Use of these established ingredients is explicit. The previous three-binary-snapshot work's scalar is not redone or represented as new, and all binary readings of the same categorical table are blind here.
3. [Dobrynin–Stroujkova, Nov 18](https://arxiv.org/pdf/2611.24898) binary persistence and sensor kinetics, [Aznagulov, Nov 20 v2](https://arxiv.org/abs/2609.14623v2) architecture-independent lower bound via sojourn-density departure from complete monotonicity, and [Raux–Maier–Seifert v2](https://arxiv.org/abs/2608.03619v2) unequal two-time cross-correlations are distinct observation contracts. [Gu 2026](https://arxiv.org/abs/2609.39205v2) estimates using mark-specific actual waiting times/excursions; no claim to improve inference with such richer information.
4. Narrow candidate contribution: a genuinely joint **finite-window** affinity separation on a cut off five-category measurement with one rare heterogeneous motif, proved without assuming knowledge of the constructed architecture in the certificate. The exact law, all-binary blindness, comparison to the optimized prior finite-window bound, state reversal, constant true fraction, and robustness are verifiable. No matching finite-window logarithmic-gap theorem found in the bounded scan; no priority certification, experiment, field-scale generalization, or breakthrough claim.

## Replay

    sha256sum -c protocol.sha256
    python3 verify.py >/tmp/multicategory-replay.txt
    cmp /tmp/multicategory-replay.txt replay.txt
    sha256sum -c manifest.sha256

Standard library only. Verifier checks exact construction identities and rational table PSD, all 64 binary cuts, censor mass upper comparison, and effect thresholds using rational interval logarithms and square roots. Infinite families, exact stationarity and reversibility are algebraic proofs above rather than extrapolated numerical results.
