# Exact ordered optimum for two amplitude-damped phase calls

## Claim, resource and prior

There is one labeled call to each of two independent qubit channels, with survival \(a\) then \(s\), where
\[
 E_{\theta,v}(\rho)=e^{i\theta P}(B_{0,v}\rho B_{0,v}^{\dagger}
                +B_{1,v}\rho B_{1,v}^{\dagger})e^{-i\theta P},\quad
 P=|1\rangle\langle1|,\quad B_{0,v}=\operatorname{diag}(1,\sqrt v),\quad
 B_{1,v}=\sqrt{1-v}|0\rangle\langle1|.
\]
The single unknown scalar is a local phase, both durations normalized to one, the two damping environments are inaccessible, and arbitrary finite-dimensional noiseless ancillary registers and phase-independent intermediate control are allowed. The order is fixed before the experiment. Every run has two uses, with no postselection. We optimize pointwise symmetric-logarithmic-derivative (SLD) QFI; local controls may be designed at the same nominal phase, not using the unknown deviation. Take the nominal phase to be zero by known rotations. The two phase and damping *channels* commute at each slot, and the signal sign/global phase used in [Liu, Hu, Yuan and Yang, PRL 130, 070803 (2023)](https://doi.org/10.1103/PhysRevLett.130.070803) is equivalent at evolution time one.

**Theorem.** For \(0<a,s\leq1\), with \(p\) the first input excitation probability,
\[
\begin{aligned}
 P_0&=1-(1-a)p,&&P_1=1-P_0,\qquad d=\frac{\sqrt{ap(1-p)}}{P_0},\\
 J_1(s)&=\frac{4s}{(1+\sqrt{s})^2},&
 V_s(d)&=\max_{0\leq x\leq1}\frac{4s(d+\sqrt{x(1-x)})^2}{1-(1-s)x},\\
 J(a\to s)&=\max_{0\leq p\leq1}\{P_0 V_s(d)+P_1J_1(s)\}.
\end{aligned}
\]
A physical two-qubit protocol (one routed sensor plus a noiseless tag, and a stored binary midpoint outcome) attains this optimum. In the inner maximization put \(x=1/(1+t^2)\). Then the maximizer is the unique positive root \(s-t^4-d(1-s)t(1+t^2)=0\). The computation is a scalar maximization with a monotone quartic evaluation, not a numerical optimization over quantum combs. If either survival is zero the optimum equals the other one's single-use optimum \(J_1(v)\), with \(J_1(0)=0\). **Order:** \(J(1\to s)<J(s\to1)\) for \(0<s<1\); identical strengths and zero-survival edges tie. For two unequal partially damped channels the formula gives certificates at tested pairs, not yet a proven uniform sign or exhaustive interior tie classification.

The earlier Research1 result establishes strict adaptive-vs-parallel separation, not the ordered optimum. [Liu et al. 2023, Supplement E.2](https://arxiv.org/html/2203.09758) give exact sequential SDP problems and numerically extract optima; [Liu et al. 2024 §§5.1.1, 5.1.3](https://arxiv.org/html/2409.07068) include unequal cases on the ray of decay probabilities \(p_1=2p_2\) and plot the maximum over the two orders. [Li–Tong–Bao 2026, Fig. 4/Eq. 13](https://arxiv.org/html/2610.31249v2) use conditional loss recycling in continuous-time three-level erasure/Ramsey sensing; the general reset idea is not new. The contribution here is an algebraic fixed-order optimum with unrestricted ancillary/control benchmark, an explicit physical saturation and certified order comparisons. No claim about optimal indefinite order, transverse noise, adaptive dynamics *within* calls, global estimation, or hardware cost.

## 1. Upper bound with ideal flags and arbitrary controls

Purify initial preparation, and take a Stinespring isometry for any phase-independent intermediate processing, retaining its ancillary environment as allowed memory. Discarding added memory cannot raise QFI. Separately refine each damping channel by an *ideal classical* Kraus outcome register. The unflagged strategy is recovered by ignoring these outcomes, and using them for controls can only increase the supremum. This is a comparison operation, not an assumption of physical access to the damping environment. Conditional on ideal outcomes, the augmented quantum state is pure.

An arbitrary purified input decomposes as \(\sqrt{1-p}|0\rangle|\alpha\rangle+\sqrt p|1\rangle|\beta\rangle\); normalized ancillary states \(\alpha,\beta\) need not be orthogonal. The no-jump sector has two orthogonal vectors because the sensor bits differ, and relative phase \(\theta\); its normalized tangent after subtracting a global phase has norm \(d=\sqrt{ap(1-p)}/P_0\). The jump branch has no phase except an irrelevant global phase under a before-damping representation. Its weight and the no-jump weight are phase independent. The jump contribution is bounded by \(P_1J_1(s)\).

For the second call after no jump write normalized input and horizontal tangent at the nominal phase as \(\psi,\ \psi'=i d\phi\), where \(\|\psi\|=\|\phi\|=1,\langle\psi|\phi\rangle=0\). This follows after removing a local scalar phase; for \(d=0\), any orthogonal \(\phi\) suffices. A phase-independent isometry may position these two vectors, but cannot enlarge \(d\). Allow all possible placements for the upper bound; a two-dimensional logical rotation attains it below. Let \(Q=P\otimes I_{\rm mem}\), \(x=\langle\psi|Q|\psi\rangle\), \(z=\Re\langle\psi|Q|\phi\rangle\), \(\ell=1-s\), and \(W=1-\ell x\). Because \(Q\) is a projection and \(\phi\perp\psi\), \(|z|\leq\|(Q-x)\psi\|=\sqrt{x(1-x)}\).

For an unnormalized pure branch \(v(\theta)\) of positive weight \(q=\|v\|^2\), its contribution to *total flagged QFI*, including the classical Fisher term, is
\[
 (q')^2/q+qF_Q(v/\sqrt q)
   =4\left(\|v'\|^2-\frac{[\operatorname{Im}\langle v|v'\rangle]^2}{q}\right).
\]
This follows by differentiating the normalization and using \(q'=2\Re\langle v|v'\rangle\). It does not discard information in changing probabilities, which need not be constant for arbitrary admissible controls. For the two branches \(v_j=B_{j,s}\psi\) and \(v_j'=iQv_j+i d B_{j,s}\phi\) with \(Qv_1=0\). Directly,
\[
 \|v_0'\|^2+\|v_1'\|^2=d^2+sx+2sdz,\quad
 \operatorname{Im}\langle v_0|v_0'\rangle=sx-\ell dz,\qquad
 \operatorname{Im}\langle v_1|v_1'\rangle=\ell dz .
\]
For \(0<s<1\), \(0<x\leq1\), the total QFI divided by four is therefore
\[
 d^2+\frac{sx(1-x)+2sdz-\ell d^2z^2/x}{W}. \tag{1}
\]
The imaginary part of \(\langle\psi|Q|\phi\rangle\) can affect individual flag derivatives, but cancels from the total formula after including classical Fisher information.

**Zero-weight cases.** At \(x=0\) the absent branch has \(v_1=0\) and hence zero density derivative even if \(v_1'\ne0\); one must omit it for pointwise SLD. The no-loss contribution is \(4d^2(1-\ell\langle\phi|Q|\phi\rangle)\leq4d^2\). A zero-weight sector can contribute to continuous Bures QFI, a distinct definition. Even granting the first-order purification norm \(4\sum_j\|v_j'\|^2\) at \(x=0\) bounds it by \(4d^2\), so this exceptional point cannot improve this optimum under either convention. At \(x=1\), \(z=0\), and (1) is \(d^2\). If \(s=1\), the jump is identically absent and the pure-state expression is \(4[d^2+x(1-x)+2dz]\); the maximum is \((1+2d)^2\). If \(s=0\), the second call is a phase-independent replacer; use a dummy input and preserve the first output in noiseless memory. The \(p=0,1\) first-call endpoints produce no phase-sensitive branch, so their value is \(J_1(s)\). At the saturating controls for positive strengths, a Kraus weight is either identically zero or \(x\) is interior; no choice of interpretation at a rank change is needed for the optimum.

## 2. Scalar saturation and physical lower bound

For \(0<s<1,d>0\), maximizing (1) over \(|z|\leq\sqrt{x(1-x)}\) gives \(z=\min\{\sqrt{x(1-x)},sx/(\ell d)\}\). The unsaturated region is \(x\leq \ell^2d^2/(s^2+\ell^2d^2)\), and its value is \(d^2+sx/\ell\), increasing up to the shared endpoint. Thus the overall maximum can be attained with saturated \(z=\sqrt{x(1-x)}\), where (1) is \(s(d+\sqrt{x(1-x)})^2/W\). If \(d=0\), this identity gives \(J_1(s)\) directly at \(x=1/(1+\sqrt s)\); if \(s=1\), the pure maximum. The \(x=0,1\) upper value \(4d^2\) poses no gap: the saturated objective at \(x=1\) is \(4d^2\), and strictly higher at its interior maximizer.

To attain the ideal flag bound, begin with \(\sqrt{1-p}|0,0\rangle+\sqrt p |1,1\rangle\) and measure **only** sensor/tag parity after the first call. This does not distinguish the logical amplitudes within no jump. Keep the midpoint flag. In the first no-loss branch rotate its nominal horizontal pair to
\[
 \psi=\sqrt{1-x}|0,0\rangle+\sqrt{x}|1,1\rangle,\quad
 \phi=-\sqrt{x}|0,0\rangle+\sqrt{1-x}|1,1\rangle.
\]
This is a parameter-independent unitary on the code at the same nominal phase; a global tangent phase drops out. In the first loss branch prepare a tagged single-call optimum. In both branches use the second labeled channel once. Both the pre-call state and its tangent lie in the tag code, and their two final Kraus sectors lie in orthogonal physical parity subspaces. Consequently the ideal-flag QFI is attained without either environment, an unknown-phase control, or postselection. One tag can stay with the same sensor routed through the two slots, or its logical sensor can be swapped into a distinct second device.

With \(t=\sqrt{(1-x)/x}\), the positive square root of \(V_s(d)/(4s)\) is
\[
 \frac{d(1+t^2)+t}{\sqrt{(1+t^2)(s+t^2)}}.
\]
Its derivative has the sign of \(s-t^4-d(1-s)t(1+t^2)\), going from positive at zero to negative at infinity and strictly decreasing for \(t>0\). Thus the unique optimizer has \(0<t\leq s^{1/4}\leq1\) and interior \(x\).

The outer \(f_{a,s}(p)=P_0 V_s(d)+P_1J_1(s)\) is concave. For each \(p_i\), coherently tag the canonical input with an additional orthogonal label \(i\), weight it by \(\sqrt{w_i}\), and superpose. The aggregate excitation is \(\sum_i w_i p_i\); reading the label after the first call and executing the corresponding continuation gives QFI \(\sum_i w_i f_{a,s}(p_i)\), which the analytic upper bound at aggregate \(p\) already covers. This proves concavity without assuming that a supremum of concave functions is concave. For \(0<p<1\) the envelope is real analytic because the quartic root is simple and interior. Each such \(p\) yields \(d>0\) and \(V_s(d)>J_1(s)\) by evaluating at the single-call maximizer, so the endpoint values are strictly smaller. The outer maximum is unique: if a concave function achieved it at two interior points it would be constant between them, contradicting analyticity and the strictly smaller endpoint values. Rational bisection isolates the quartic root and bounds the objective and its \(p\) derivative. Positive/negative slope enclosures bracket every global optimum; a supporting tangent then bounds the maximum within the bracket. No floating-point global optimizer is used in the interval results.

## 3. Order edges and exact benchmarks

The formula computes either fixed order. If the first use is noiseless, \(P_0=1\) and \(d=\sqrt{p(1-p)}\) is maximized at \(p=1/2\); \(V_s\) increases with \(d\). With \(r=\sqrt{s}\) the value is the maximum over \(t\geq0\) of
\[
 \frac{r^2(1+t)^4}{(1+t^2)(r^2+t^2)}\leq(1+r)^2 .
\]
Indeed \((1+t^2)(1+t^2/r^2)\geq(1+t^2/r)^2\), and Cauchy–Schwarz gives \((1+t)^2\leq(1+r)(1+t^2/r)\). For \(s\to1\), the feasible \(p=1/2\) gives \(1+2r+2r^2/(1+r^2)\), strictly greater than \((1+r)^2\) for \(0<r<1\). This classifies the edge with either survival one; both corners are included by equality for \(s=0,1\). A zero-survival slot is a phase-independent replacer, so the orders tie at the one-call optimum; equal survivals tie trivially. A uniform sign for *all* unequal, partially lossy pairs is not claimed: a grid or single certified comparison cannot prove the interior tie boundary.

At \(a=16/25,s=81/100\), choose rational \(p=1/2,t=9/10\), with first weight \(41/50\), tangent \(d=20/41\), and the optimal jump-branch reset. The exact physical QFI is \(203724139/66974525\), over **8.4937%** above the exact arbitrary-ancilla parallel optimum \(5184/1849\) from the earlier proof. Both branches, total probability one, and the accessible SLD are replayed directly. The global optima obey these outward brackets:
| Order | Certified optimum QFI interval |
|---|---:|
| 0.64 then 0.81 | [3.04612987458271, 3.04612987458272] |
| 0.81 then 0.64 | [2.95638347987805, 2.95638347987806] |
| 0.60 then 0.60 | [2.55737781219901, 2.55737781219902] |
| 0.50 then 0.50 | [2.17926669798591, 2.17926669798592] |
| 0.25 then 1.00 | [2.41835635011737, 2.41835635011738] |
| 1.00 then 0.25 | [2.03112583027372, 2.03112583027373] |

The unequal order gap is in [0.08974639470466, 0.08974639470467]. These are directed rational interval results, not floating-point error bars. Independent direct SLD controls for an admissible complex overlap give total QFI \(5679/4525\) and nonzero flag CFI \(171/4525\); another control has a zero-mass branch with nonzero ket derivative, omitted under SLD as required.

The [pinned author sequential comb SDP dataset](https://github.com/qiushi-liu/strategies_in_metrology/blob/07748a4e94d7b29dc888d2cbf41d0b0f421c58dc/data/QFI_seq_rangep_0_1_N_2.txt) is pinned by commit and independently checked by SHA-256. The author solver uses \(N=2\), decay \(p_{\rm decay}=1-s\), phase evolution time one and nominal phase one; absorb the latter into known input/control rotations, preserving optimized local QFI. The inner variable \(t\) here is not their evolution time. Rows are decay \(j/200\); author code was not executed. Our independently evaluated scalar curve agrees with the 199 interior source SDP rows to maximum absolute \(7.962\times10^{-10}\), including the rounded PRL 2.179 at decay .50. This is a numeric comparison to author outputs, not a fresh MOSEK solve or a proof by numbers. The 2024 unequal figure runs only at their specific decay ratio and takes the maximum over orders; it is not a direct external verification of this particular unequal pair.

## Replay

From this directory run `python3 verify_ordered.py` and `python3 physical_sld.py` (Python standard library only). The exact interval replay is the first nine lines of `replay.txt`. To also reproduce the comparison against the authors' published numerical values, run `python3 fetch_sdp_reference.py` and then `python3 verify_ordered.py --sdp`. Only the pinned data file is downloaded, and its SHA-256 is checked before use. No third-party code is executed.
