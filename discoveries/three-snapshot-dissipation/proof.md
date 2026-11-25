# A sharp universal hidden-dissipation certificate from three binary readings
Research note, 24 November 2026. One finite-time inference result, not a new conventional current TUR, not a thermodynamic measurement, not a priority or breakthrough claim.

## Relation to the closest prior calculation

Petkov’s [November 19 manuscript, §9, eqs. (9.1)–(9.8)](https://hexagonmath.org/pdf/2611.00657v1) already gives the same scalar formula and two rank-one optimizers for marked-excursion counts. For the conditional snapshot table `R=[[1−2s+r,s−r],[s−r,r]]`, use `w=2(s−r)` and take half his excursion count term `2w φ(2 min(√(R11 R22),√(R12 R21))/w)`. This is exactly `B(s,r)`, including the zero branch. The contribution here is using this expression on three fixed-lag binary observations, the finite-lag comparison on those same readings, and the heterogeneous Markov sharpness family. It is not a new scalar formula.

## Claim

Let X be a stationary finite homogeneous Markov jump chain, with stationary support pi_i>0, reciprocal positive rates on its allowed edges and ordinary all-even reversal. All logs are natural; set k_B=1. Let Z=h(X) be a fixed deterministic binary readout. Only the joint law at times -t, 0, +t is used; t>0, and intervening jumps are not counted. For each visible midpoint z with nonzero probability p_z, define

    s_z = P(Z(+t) != z | Z(0)=z) = P(Z(-t) != z | Z(0)=z),
    r_z = P(Z(-t) != z, Z(+t) != z | Z(0)=z).

The equality follows from stationarity and equality of the two off-diagonal entries in a stationary binary two-time table (no equilibrium assumption). Thus each conditional two-by-two matrix is symmetric. The process may even have an exactly time-reversible full visible trajectory. Define

    g=s^2-r;
    B(s,r) = 0                                      for g<=0,
    B(s,r) = 2 sqrt(g) atanh(sqrt(g)/(s-r))           for g>0.

The dimensionless steady-state entropy-production rate obeys

    sigma >= (1/t) sum_z p_z B(s_z,r_z).             (1)

A finite irreducible chain with both labels present has strictly positive transition probabilities at every t>0, so the argument of atanh is in (0,1) on nonzero blocks. Use extended lower-semicontinuous boundary values for degenerate Bernoulli data. Multiply the right side by k_B for physical entropy units; equating microscopic path divergence with physical entropy production assumes local detailed balance and the stated reversal. Reciprocal support can be dropped if the full rate is infinite, which is a trivial separate regime.

There are two distinct sharpness statements: (i) B is exactly the minimum of an explicit two-moment latent product Jeffreys-cost problem, with a two-state optimizer; (ii) the prefactor of (1) cannot be increased uniformly over the stated CTMC class (positive ratios approach equality). This does not prove a data-dependent CTMC optimum at an arbitrary fixed nondegenerate table or exact finite-rate equality.

## Proof

Use a row generator Q. Its stationary reverse has off-diagonal rates q^*_{ij}=pi_j q_{ji}/pi_i and exactly the same exit rates by stationarity. At microscopic midpoint i, let u_i and v_i be the probabilities of the future or reversed-past single half-window event of ending in the opposite visible label. Given i, the two half-paths are independent by the stationary Markov property. Thus the 2x2 two-event law given i is a product Bernoulli law, with marginals u_i and v_i; averaging within z gives the conditional observed table.

For half-path laws conditional on starting at i, relative entropy depends on sums of log rate ratios across jumps (holding terms cancel). Averaging the forward-to-reverse path KL over pi gives
t times sum_{i!=j} pi_i q_{ij} log(pi_i q_{ij}/(pi_j q_{ji})) = t sigma;
the reverse-to-forward KL gives the same value by index exchange. Data processing to the single binary endpoint event in both directions yields
t sigma >= sum_i pi_i J(u_i,v_i)/2, where J is the Bernoulli Jeffreys divergence (sum of the two KLs).

Here is the sharp conditional-product minimization. Within one midpoint block use expectation over its microscopic states, drop i, and write
a=(1-u)(1-v), r0=uv, b=u(1-v), c=v(1-u), w=b+c. Then bc=a r0, u-v=b-c, and

    J(u,v)/2 = (b-c) log(b/c)/2 = w phi(2 sqrt(a r0)/w),
    phi(x) = y atanh(y), y=sqrt(1-x^2), 0<=x<=1.

Terms with w=0 contribute zero. The affinity–Jeffreys function and its convexity are already given in [Petkov, Appendix F, eqs. (98)–(99)](https://hexagonmath.org/pdf/2611.00658v1). The application here is to the disagreement-conditioned Bernoulli product, followed by the observed-diagonal constraint. On (0,1), phi is decreasing, and
phi''(x) = (y/(1-y^2) - atanh(y))/y^3 > 0
by integrating the strictly increasing function 1/(1-y^2) from zero to y. Endpoint cases follow by limits. Jensen with weights proportional to w and then Cauchy-Schwarz give

    E[J/2] >= E[w] phi( 2 E[sqrt(a r0)] / E[w] )
           >= E[w] phi( min(1, 2 sqrt(E[a] E[r0]) / E[w]) ).

The observed inputs are E[r0]=r, E[a]=1-2s+r, E[w]=2(s-r). The determinant identity
(s-r)^2 - r(1-2s+r) = s^2-r
gives B. Weight the block inequalities by p_z to obtain (1).

Sharpness of the latent product problem: if g>0, take equal weights for (u,v)=(s+sqrt(g),s-sqrt(g)) and its swap. Entries are admissible because r>=0 and 1-2s+r>=0. The choices have exactly the same a, r0, w and phi argument, so both inequalities are equalities. For g<=0 and s,r>0, set u=v=r/s with probability s^2/r and u=v=0 otherwise, which gives exact first moment s and second moment r at cost zero. If s=0 the result is trivial. Even without equal forward and backward marginals the proof applies using the **symmetrized** two-event table; its optimality statement then applies to this symmetrized summary. No conditional probabilities of hidden states are measured.

## Comparison on identical binary inputs

- The midpoint/Pinsker lower bound is 8 sum_z p_z [s_z^2-r_z]_+/t. Apply Pinsker to the equality-realizing latent products (or Jensen directly): J/2 >= 2(u-v)^2. Thus (1) is at least as large. At fixed rare s with g small, B=2g/[s(1-s)]+O(g^2), so the gain need not be small. Existing finite-window bounds for actual **no-jump persistence** use a different observable; we compare their Pinsker step applied to the identical three-label data, not all bounds with richer records.
- Compare Petkov's author-circulated Theorem 4.1 for the same matched classifier (joint category: midpoint label, endpoint differs). On each binary block the CP deletion deficit is d=[2(s-r-sqrt(r(1-2s+r)))]_+. Indeed the two diagonals can be kept, and the off-diagonal must be at most their geometric mean; this necessary 2x2 PSD constraint is attainable by a nonnegative Gram matrix. Zero cross-block entries imply the total deficit is the weighted sum of these defects.
  Set w=2(s-r), x=2 sqrt(r(1-2s+r))/w, F(v)=sqrt(v)atanh(sqrt(v)). For a violation d=w(1-x), B=w F(1-x^2), and w<=1/(1+x) from E[a]+r>=2 sqrt(E[a]r). At w_max=1/(1+x), B=2F(d_max) by the tanh double-angle identity. At smaller w, strict convexity of F and F(0)=0 give B>=2F(d), strictly for g>0 and s!=1/2. A second Jensen over blocks gives (1) at least as large as the CP floor on this chosen classifier. This does not compare to CP using more labels or longer motifs.

## Physical sharpness and exact controls

Consider nodes 0,1,2,3 on a cycle with labels 0,0,1,1. On each edge from i to i+1 give rate H+j on the within-label edges and k+j on the two crossing edges. Opposite rates are H-j and k-j, respectively, with 0<j<min(k,H). The stationary law is uniform (every row and column sums to zero), and every microscopic exit rate is H+k. The label-preserving permutation i -> 1-i (mod 4) reverses all arrows while preserving edge types; it conjugates Q to Q^T. Therefore every observed finite-dimensional law is time reversible, not only the low-order table. The half-turn i -> i+2 interchanges visible labels and leaves dynamics invariant. In either block, the two microscopic forward crossing probabilities u and v are swapped under reversal, so s=(u+v)/2, r=uv and the latent-product step is exact already.

Elementary stationary edge fluxes yield
sigma = j log((k+j)/(k-j)) + j log((H+j)/(H-j)).
At fixed H, u = (k+j)t+O(t^2) and v = (k-j)t+O(t^2). Then B/t tends to the first (crossing-edge) term exactly. In the heterogeneous family with k=1, j=1/2, H->infinity and t=H^-3, the error in u/t or v/t is O((H+k)^2 t)=O(1/H), by uniformization including the terms of order two and higher and the single-jump factor. The hidden-edge entropy contribution tends to zero, and the ratio of the positive bound to sigma tends to one. Rescaling all rates puts any chosen positive lag on this sequence without changing the ratio; this is not saturation at fixed signal magnitude and bounded traffic.

For comparison, at fixed k,j,H and t->0, the CP floor approaches 4(k-sqrt(k^2-j^2)), and the Pinsker floor tends to zero. As j approaches k from below (H>k), the ratio of the new floor to the same-input CP floor can thus be arbitrarily large in finite reciprocal CTMC controls by taking t small enough at each finite j. This limiting gain is not a statistical-power statement; triple double-crossings then get rare.

The exact rational replay checks stationarity, the generator conjugacies, and enclosures for finite-time probabilities with a Poisson tail. Fixed example H=10,k=1,j=0.95,t=0.0001 (consistent inverse-time units): true sigma 3.66142952745; (1) 3.45679138912 (94.41%); same-input Petkov CP floor 2.73865902409; same-input Pinsker floor 0.000720413532. For H=10,k=1,j=0.5,t=0.001, respective values 0.599347873613, 0.537101967770 (89.61%), 0.524327001422, and 0.00195653264347. All displayed values rounded. Control j=0 gives exact zero. A strictly positive covariance product-mixture control is realizable at zero relaxed cost. Blind direction: an independent driven hidden ring tensor a reversible observed two-state chain has strictly positive total sigma but zero B exactly, since the two halves of the observed chain are conditionally iid given its visible midpoint. Thus zero certificate never proves reversible hidden dynamics.

## Originality audit and limits

1. [Horowitz–Gingrich, PRE 2017](https://doi.org/10.1103/PhysRevE.96.020103) treats current precision at finite times, and [Hasegawa–Van Vu, PRL 2019](https://doi.org/10.1103/PhysRevLett.123.110602) handles time-antisymmetric fluctuation-theorem observables; neither is a dissipation bound from a reversible binary triplet. [Buschmann–Seifert–Maier 2026](https://arxiv.org/abs/2608.14276) extends response/counting bounds, including kinetic traffic; [Vu et al., Nov 19, 2026](https://arxiv.org/abs/2611.26267) gives finite-time martingale response identities. These are genuine prior TUR/kinetic refinements, not claims this calculation improves their own information contract.
2. [Skinner–Dunkel, PNAS 2021](https://doi.org/10.1073/pnas.2024300118) already gives the same two-moment variational program (supplement A18–A20), its balanced-marginal case (A28–A29), and a finite-sampling reduction (A39–A50). The explicit general binary formula is also prior, as explained above. The [independent audit](audit.md) spells out the correspondence. [Van der Meer–Ertel–Seifert, PRX 2022](https://doi.org/10.1103/PhysRevX.12.041026) and [Gu 2026](https://arxiv.org/abs/2609.39205v2) infer from resolved transition waiting-time/excursion statistics. [Stroboscopic generator reconstruction, 2412.15642v2](https://arxiv.org/html/2412.15642v2) reconstructs or bounds from fully resolved state propagators, not a 2-label projection. [Raux–Maier–Seifert 2026 v2](https://arxiv.org/abs/2608.03619v2) bounds with asymmetric finite-lag correlations; those asymmetries vanish for the reversible observed control.
3. Closest indexed prior: [Dobrynin–Stroujkova, Nov 18, 2026](https://arxiv.org/pdf/2611.24898), Sec. S5.3 derives midpoint/Pinsker linear finite-window no-switch persistence and Sec. S8 gives the sign obstruction for finite-lag autocorrelation. Their short-time persistence rate bound already reaches the crossing-edge expression on our symmetric limiting family. We do not claim the midpoint construction, binary sign obstruction, four-state architecture, or infinitesimal coefficient as new. The application to three stroboscopic values and the saturating heterogeneous limit are the specific additions.
4. Closest general prior: [Samuil Petkov, author-circulated PDF](https://hexagonmath.org/pdf/2611.00658v1), Theorem 4.1 establishes the broader mirrored-class CP floor, and J.5 already uses the homogeneous four-cycle, the same label-preserving reflection, and a persistence deletion deficit. This is the author-hosted manuscript, not an arXiv identifier. The weighted conditional Bernoulli cost above strictly exceeds the explicit CP floor for rare asymmetric binary blocks, while that paper handles many categories and other forms of obstructions. [Aznagulov, Nov 20 revision](https://arxiv.org/abs/2609.14623v2) bounds observable sojourn distances using an arithmetic-reversible comparison; that is not this exact finite-lag product problem. An author-affiliated [Evidence Press technical release](https://www.evidencepress.org/releases/hidden-dissipation-floor/) advertises a persistence mutual-information certificate, also discussed in Petkov J.5. See the closest-prior clarification above for the identical scalar optimizer in the marked-excursion setting.
5. Need stationary data, no latent odd parity or time-dependent drive that changes the comparison, and readouts aligned at the same lag. Independent identical randomized endpoint classification is safe under data processing; coupled detector memory, missing-not-at-random selection or temporal blur need a different model. The theorem uses population probabilities. Selecting an informative lag on the same noise or treating overlapping triples as independent can cause false positives; rigorous application requires independent stationary windows or justified dependence control and multinomial confidence/coarsening treatment. Sparse double-crossings can make the best population floor impractical to estimate. No device design or empirical measurement is claimed.

## Replay

    python3 verify.py > /tmp/dissipation-replay.txt
    cmp replay.txt /tmp/dissipation-replay.txt

The code is standard-library only. The [independent audit](audit.md) checks the same inputs and proves that positive-dissipation finite models have a strict lower bound, with equality of the positive ratio attained only as a limit. Rational enclosures control the independently evaluated finite case; asymptotic optimality and dominance are the proofs above, not fits.
