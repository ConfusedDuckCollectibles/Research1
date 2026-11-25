# Independent audit of the final size-78 closure

**Verdict:** the linked-layer argument and its last finite gate exclude every size-78 spectral pair over F13³ when combined with the independently verified earlier branch reductions. No publication action. This is not a theorem for other primes, and the full ambient F13³ implication additionally uses the audited separate cardinality ledger.

## Labelled six-point graph atlas, independently regenerated

Given an affine abscissa support S of six distinct field points, take its least 13-bit mask under the 156 affine horizontal maps. Of all binom(13,6)=1,716 six supports, this leaves fourteen orbit representatives. Our graph DFS does not use the author's list of five-direction palettes or constraints from an origin secant. In each support, subtract the affine height through its first two (distinct) abscissas, giving first heights 0,0. Normalize the first subsequent nonzero residual to 1 or allow all zero. Build successive remaining heights, maintaining the set of slopes of all pairs already inserted; only prune upon more than five slopes. All fixtures are labelled by their fixed ordered abscissas. An all-zero residual is collinear and invalid. The earlier cap-five zero gate excludes a graph with fewer than five directions; the present gate concerns exactly five.

Exactly five support masks have such a nonzero graph, and each has one reduced residual shape:

| Support mask | Increasing abscissas | Reduced heights | Slope mask |
|---:|---|---|---:|
| 111 | 0,1,2,3,5,6 | 0,0,1,0,9,3 | 6275 |
| 187 | 0,1,3,4,5,7 | 0,0,1,8,10,1 | 2693 |
| 215 | 0,1,2,4,6,7 | 0,0,1,0,0,3 | 203 |
| 311 | 0,1,2,4,5,8 | 0,0,1,4,1,4 | 1411 |
| 427 | 0,1,3,5,7,8 | 0,0,1,10,1,10 | 2693 |

Each normalized height array h(0)=0 is uniquely a⋅h_reduced+b⋅s with nonzero a and arbitrary b. These 156 heights project to 156 palettes with one height above each for 111,215,311; and to 78 palettes with exactly two heights each for 187,427. In the two-model case every palette has the same labelled difference h(5)−h(3) across both models. The independent script tests this directly for every palette and computes every stable pair: (1,7),(3,5) for 187; (0,8),(1,7),(3,5) for 427. One fixed ordered pair can therefore be chosen depending only on the support, not on which palette is later observed. For the three one-model masks any pair works. Optional algebraic explanation: the common reduced palette 2693 is {0,2,7,9,11}; its nontrivial affine symmetry is slope reflection k↦9−k, fixing 11. The reduced secant from 3 to 5 has slope 11 in both templates, explaining why its vertical difference survives palette ambiguity.

An actual support S need not be the representative. Fix once an affine map s=a v+b from its representative S₀ to S. Then the graph over v has slope palette a I if the graph over s had palette I, because the height difference is unchanged and the denominator is divided by a. Subtract the height at the image of zero before comparing to the atlas. The selected stable pair transfers to the same two actual levels for all palettes. No per-plane changes of horizontal scale, point relabelings, or independent layer translations occur here.

## Link-to-line contradiction

Obtain perpendicular layer normals from the verified residual lemma. In E's vector coordinates, take first basis vector equal to the normal of A’s layered family, choose a second vector in the kernel of the normal covector of E's layered family, and take the third basis vector where that covector is one. Then take the dual basis for A. A layered normal is nonuniform for its own set, hence is not a difference direction of the opposite set. These omissions make both projections injective along the respective within-layer axes, so each layer is a full graph over a field coordinate, not a partial graph. Write graphs as E=(f_s(y),y,s), A=(t,u,g_t(u)), with six-element level sets S,T.

For r∈F13 let I_r be the set of all i such that (i,r,1) is the normalized direction of a cross-layer E difference. Orthogonality forces
\[
\sum_{t\in T} \zeta^{it} C_t(r)=0\qquad(i\in I_r),\qquad C_t(r)=\sum_{u}\zeta^{ru+g_t(u)}.
\]
If the original E pair has layer gap d≠0, applying the cyclotomic automorphism ζ↦ζ^{d^{-1}} gives exactly this common array for every d. Thus all differences, not just those coming from a single layer pair, are covered. At an active r, at least one C_t(r)≠0; with a(r) active coefficients, the prime Fourier nonzero-minor theorem bounds |I_r|≤a(r)−1≤5. This also follows from its uncertainty formulation. At least two distinct r are active: a cap-five 13-point A graph cannot be affine and determines at least five directions (in fact ≥10 by the companion planar gates). A graph slope −r occurs exactly when u↦ru+g_t(u) fails to be injective, and for thirteen inputs the cyclotomic sum vanishes exactly when the thirteen values are a permutation. The sharper earlier bound of twelve active slopes is valid but unnecessary for the final contradiction.

For each active r and every y₀, the affine plane y−rs=y₀ in E contains exactly one point from each of the six layers. Viewed as a graph over s∈S, its height is h(s)=f_s(y₀+rs), and all its slopes lie in I_r. By the six-point direction gate and cap five, it has at least five and thus exactly five directions, giving the *same* palette I_r for every y₀. Consequently S is one of the five affine support classes. Pick the stable pair s,s′ according to S as above. For each active r,
\[
 f_s(y_0+rs)-f_{s'}(y_0+rs')=c_r\qquad\text{for every }y_0.
\]
The constant may depend on r (not on y₀); the pair does not. Put v=y₀+rs and subtract the conclusions for two distinct active r. This shows f_{s′}(v+δ)−f_{s′}(v)=c for every v, with δ=(r_1−r_2)(s′−s)≠0. Over the prime field the successive multiples of δ run through the whole field, so f_{s′}(v)=f_{s′}(0)+(c/δ)v. Even if c=0 the result is affine. The corresponding thirteen points lie on one affine line, contradicting cap five. This finishes size 78.

## Dependency checklist and full field scope

No circular use of size78 nonexistence: (1) plane cap thirteen, no K14, universal line cap six, and K≥12⇒K13 with a layer normal follow by elementary projection/pencil arguments; (2) presence of any six-line forces both layered, giving cap-five for the neither/one-layer analyses; (3) the independently replayed cap-11, K10, K11 atlas and weight/degree argument exclude neither-layer; (4) the independently replayed AP chart and deficit argument exclude exactly one layered; (5) universal no-six-line and existence of a perpendicular pair were independently checked, including three small cap-two gates and all 21 group-size partitions; (6) link coordinates, Fourier-minor union, and labelled atlas contradiction above close the only remaining branch. The optional bound of 12 or 13 directions within each linked layer orientation is not required for this final step. Earlier proofs/readback replays are in the separate nonlayer78 and mixed78 private directories.

For the complete F13³ ledger, the classical finite-field size/tiling reductions (e.g. [Ferguson–Iosevich §8](https://arxiv.org/html/2610.30712v1) and [Aten et al.](https://arxiv.org/abs/1509.01090)) leave only 13m for 2≤m≤12. Cardinalities 1,13,169,2197 are the standard common spectral/tiling cases; tiling implies spectrality in dimension three by the same sources. [Malikiosis Theorem 1.4](https://doi.org/10.1007/s43670-023-00078-2) excludes 130,143,156 via 169−12√13<126. The classical Brock/de Launey determinant parity removes 65 and 91 by powers of 5 or 7 congruent to −1 mod13. Separately reconstructed exact arguments remove 26,39,52,104,117. The present argument removes 78. Thus the complete F13³ conclusion follows only when those separate proofs and the standard reductions are included, not from this six-point check on its own.

Two source cautions: [Ferguson–Iosevich Lemmas 8.3–8.5 and Proposition 8.6](https://arxiv.org/html/2610.30712v1) directly supply the large-pencil/weighted/moment framework (Lemma 8.3 already holds for general p), while the odd-prime complete theorem in that paper is over F11. [Tao’s proof of Chebotarev’s prime Fourier-minor/uncertainty theorem](https://arxiv.org/abs/math/0308286) is classical primary precedent, not new. Targeted public research located no matching primary 78-point exclusion or complete F13³ proof, but is not a priority certification.

Independent replay: use independent.js as plain JavaScript and compare stdout to independent-expected.json. It regenerates all fourteen support classes, all successful reduced arrays, every affine vertical lift and every palette-stable pair, rather than trusting a stored count.
