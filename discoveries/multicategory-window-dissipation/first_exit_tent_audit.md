# Independent audit: doubling the active first-exit dissipation certificate
AI-assisted independent proof and replay, November 25, 2026. The generic matched-classifier theorem in proof.md is unchanged. The stronger result here applies to the stated **actual first exits from one observed block H**, not all classifiers. It uses the same censored finite-window categorical probabilities as the existing five-cycle example. No actual duration observations or uncensored trajectories are required for the bound.

## The exact specialization

Let X be a stationary finite, irreducible, reciprocal CTMC, every microscopic state even under physical reversal, with observed instantaneous labels. Let H consist of the microstates with one specified observable label, and V its complement. Nonzero stationary mass pi_H is observed. From a stationary midpoint in H, observe the first exit from H in each outward direction up to time tau, returning the corresponding visible label or a censor code if that half does not exit. Let M be the conditional table given midpoint H, and restrict its indices to **genuine active labels** (the active matrix is not renormalized). Multiple microstates may share any active visible label. Then

    tau*sigma >= 2*pi_H * inf_G sum_{i<j active}
                        (M_ij+M_ji) phi(2G_ij/(M_ij+M_ji)),                (A)

where phi(x)=sqrt(1-x^2)*atanh(sqrt(1-x^2)), positive weights only, and G is completely positive, with

    0 <= G_ii <= M_ii,       0 <= G_ij <= sqrt(M_ij M_ji).

Censor terms are not doubled. This is the active-only inequality at twice the prefactor of the arbitrary matched classifier result, with diagonal inequality instead of equality. In this five-cycle witness the trace argument only needs diagonal upper bounds, so the exact previous analytic lower bound doubles.

### 1. Microscopic excursion measure and sign audit

An oriented microscopic excursion starts with a V-to-H transition from state a to x1; it follows hidden states x1,...,xk with dwell lengths d1,...,dk>0 and ends at the next H-to-V transition into b, counting no holding time in visible states. Let L be the sum of the hidden durations. Its stationary intensity density on the duration coordinates, measured per unit physical time, is

    R(gamma)=pi_a q_{a,x1} [prod_{j<k} q_{xj,x(j+1)}] q_{xk,b}
                         * exp(-sum_j r_xj d_j),

where r_x is the total escape rate. The reversed excursion gamma^R has reversed states and dwell coordinates; let Rrev(gamma)=R(gamma^R) be the push-forward rate density of these reversed records in the same coordinates. Reciprocal support makes the forward and reverse measures mutually absolutely continuous. The exponential factors cancel in the ratio. By telescoping the stationary pi factors, log(R/Rrev) is exactly the sum of the edge log stationary-flux ratios along the excursion, with no missing waiting-time entropy or visible-endpoint term.

The stationary EPR is the mean rate of those edge scores. The hidden-containing excursions partition all jumps touching H; direct V-to-V jumps are treated separately. Thus the integrated oriented excursion cost plus the direct V jump cost equals sigma. Finite recurrent CTMCs have finite mean excursion jump count and bounded edge log ratios on their allowed rates, justifying the sums and integrals. Each unordered pair of distinct observed endpoint labels i,j has a **nonnegative** contribution

    integral_{gamma: i->j} (R(gamma)-Rrev(gamma))*log(R(gamma)/Rrev(gamma)) dgamma.

There is no factor 1/2 after one representative orientation i<j is chosen. The sum in that expression includes the reversed contribution formerly in the j->i integral. Diagonal endpoint classes are themselves closed under reversal, so their combined contribution is nonnegative by symmetrization; direct visible jumps also pair nonnegatively. Dropping classes is allowed only after this pairing, not on the signs of the individual oriented log terms.

### 2. Exact finite-window selection is a tent, not a sampled duration

A particular excursion of length L supports the joint observation of both endpoints precisely at midpoint locations s obeying 0<s<L, s<=tau and L-s<=tau (endpoint equality has measure zero). Their length is

    w_tau(L) = max(0, min(tau,L) - max(0,L-tau))
             = max(0,min(L,2*tau-L)), in [0,tau].

The stationary random-origin probability of the active table entry is the intensity weighted by this length:

    pi_H M_ij = integral_{gamma: i->j} w_tau(L) R(gamma) dgamma.

This follows by time-location integration of each occurrence interval, stationarity and Tonelli. It is an identity for the experiment, not a law substituted for the measurements. The same tilt multiplies Rrev, as excursion length is invariant under reversal. Since w_tau<=tau, the weighted paired costs of all distinct active pairs together cannot exceed tau*sigma. The comparison remains true for alternative compatible chains with arbitrarily many microscopic states and arbitrary hidden-to-hidden edges.

### 3. Why the tilted geometric mean is completely positive

For hidden midpoint microstate x, write P_x(i) and P_x^*(i) for **subprobability measures on microscopic stopped half-paths** from x to the first state of visible label i, with length at most tau, under Q and stationary reverse Q*. Let

    h_x(i) = integral sqrt(dP_x(i) dP_x^*(i)).

The full microscopic stopped word at midpoint x has conditional measure P_x^*(left)*P_x(right) under the forward stationary ensemble; reversing it changes it to P_x(left)*P_x^*(right), in the same outward coordinates. Multiply both by pi_x. The integral of their geometric mean over all marked stopped words of endpoint labels i,j therefore factors:

    integral_{gamma: i->j} w_tau(L) sqrt(R(gamma) Rrev(gamma)) dgamma
      = sum_{x in H} pi_x h_x(i)h_x(j)
      = pi_H G_ij.                                                        (B)

To check the identification rather than only a formal factorization, insert the midpoint inside a dwell of excursion gamma. Split that dwell into its left and right times. The map between the pair of half-path time coordinates and the full hidden durations plus s is linear with absolute Jacobian 1 on each such dwell cell. The product of the left reversed and right forward rate densities times pi_x telescopes to R(gamma), independent of s; the other orientation telescopes to Rrev. The allowed s set has precisely the length above, proving (B) including multiple hidden-to-hidden jumps and repeated hidden states.

The last Gram representation is CP. Cauchy-Schwarz on the full stopped-word records gives G_ij<=sqrt(M_ij M_ji), including G_ii<=M_ii. The diagonal need not saturate, because two microscopic records with the same endpoints may have different forward/reverse densities. Jensen with the pair weight w_tau(L)(R+Rrev), using (v-u)log(v/u)=2(v+u)phi(2sqrt(vu)/(v+u)), gives (A). Zero-weight terms are omitted and a positive-weight infinite limit is handled by monotone extension.

### 4. Strict same-data comparison and scope

The five-cycle trace/copositivity estimate of proof.md needs only Gtilde_ii at most (1,2,2,2,1), not equality. Therefore the lower bound (3) for active pairs becomes 2*p_H*p^2*alpha*(1+log(2/e))/tau on the EPR rate. The old midpoint-stratified Petkov Eq. 11 ceiling does not change. No rate parameter or microscopic graph is assumed known when constructing the lower certificate from the table. The actual 15-state model and its EPR are unchanged.

| Finite same-data comparison at e=1/32 | Audited conservative value |
|---|---:|
| New specialized EPR-rate lower bound | at least 0.0001876917017181 |
| Optimized midpoint-stratified prior Petkov Eq. 11 EPR-rate ceiling | at most 0.00004123477217914 |
| Guaranteed factor over prior | at least 4.55178219253 |
| True microscopic dissipation fraction certified | at least 27.9420235595% |
| Original path reversal and binary cuts | 0, unchanged |
| Existing uncensored excursion-rate bound on a richer raw record | at least 88.1456% of true rate, unchanged |

As e->0 at fixed tau=4 log 2, the true EPR fraction certified by the analytic bound tends to p^2/(lambda tau)=0.31700... (about 31.70%). The logarithmic improvement ratio's leading coefficient doubles. A physical entropy interpretation still requires the usual local detailed balance and ordinary even-state reversal assumptions.

**Noise caveat:** The generic matched-classifier perturbation calculation in proof.md remains correct, but (A) is not automatically robust to arbitrary relabeling of the censor code into an apparent active exit. Such a record need not sit inside a full excursion and would not have the stated exact weight. Censor-creating deletion or equal relabeling of *genuine active* endpoints can be included via a common reversal-compatible mark kernel; no additional doubled numerical noise claim is needed or provided. Nor is (A) a statement about simple endpoint snapshots that can exit and then reenter H between sampled times.

### 5. Prior and independent check

[Petkov Hexagon 2611.00658v1](https://hexagonmath.org/pdf/2611.00658v1) Eq. 13 gives the conditional midpoint Gram for general finite classifiers, and Appendix J.2 uses nonnegative paired affinities and a CP geometric Gram for **uncensored complete-excursion rates**. [Petkov Hexagon 2611.00657v1](https://hexagonmath.org/pdf/2611.00657v1), §10, already uses a duration tilt (Laplace) for complete excursions, and §9 gives the scalar equality/factorizations. The independent derivation uses those ingredients but checks the triangular `w_tau(L)` exactly from the unchanged compressed finite-window observation and factors the geometric measure by marking and splitting at the stationary midpoint. A narrow search did not identify that same finite first-exit specialization; this is not a priority claim.

The accompanying standard-library replay verifies the scale inequalities by rational enclosures and independently checks the finite-window/tent identity coefficient by coefficient in a distinct nonuniform four-state ring with internal H--H transitions and direct V--V jumps. In that ring the ordinary killed-path, stationary reverse, and geometric killed generators are rational and nonsymmetric in the ordinary Euclidean inner product. This test complements rather than substitutes for the exact path proof. The general matched-classifier and specialized certificates are distinct statements.
