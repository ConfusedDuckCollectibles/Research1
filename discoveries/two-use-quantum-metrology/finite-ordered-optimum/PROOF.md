# Arbitrary finite ordered amplitude-damped phase sensing

## Claims and model

For each labeled call with known survival \(s\in[0,1]\), the qubit experiences
\[
 \mathcal E_{\theta,s}(\rho)=e^{i\theta P}
     (B_0\rho B_0^\dagger+B_1\rho B_1^\dagger)e^{-i\theta P},
\quad P=|1\rangle\langle1|,\quad
 B_0=\operatorname{diag}(1,\sqrt s),\quad B_1=\sqrt{1-s}|0\rangle\langle1|.
\]
All calls encode the same scalar phase, have independent inaccessible damping environments, and are used once in the given fixed order. Controls outside the slots are parameter-independent CPTP operations, may use ideal noiseless memory, and may be locally designed at nominal \(\theta=0\). Measurements and CPTP resets with transient entropy sinks are allowed. We optimize **pointwise SLD QFI**, with no postselection, neither access to the model's environment nor controls inside a slot. This is not a model of a closed finite-dimensional unitary sensor or of globally unknown phase.

1. **Exact local recursion, finite horizon.** The suffix needs only the incoming horizontal squared tangent \(q\), not a high-dimensional density matrix. A two-by-two compression of its next excitation projector suffices, and the optimum reduces to two real scalar controls. It is attained using a routed sensor, two reusable noiseless qubits and classical feedback. This is an upper construction, not a necessity bound on ancilla count.
2. **The two-call saturation rule fails strictly in later genuine horizons.** For any incoming \(q>0\), a partially lossy next call \(0<s<1\), and any later suffix containing a positive-survival call, maximal real overlap \(z=\sqrt{x(1-x)}\), which erases all incoming phase tangent on loss, is strictly suboptimal. In fact an optimal real reduced control has \(0<z<\sqrt{x(1-x)}\) and nonzero conditional tangent after loss. An exact three-call witness with ordered survivals \(1,1/4,1\) exceeds every middle-stage saturated strategy by more than \(0.0337\) QFI, without computing the unrestricted optimum.
3. **Constructive known-slope attainment with bounded quantum memory.** For \(N\) identical partially lossy calls at fixed \(s\), a deterministic local first-jet update uses at most three total reusable quantum qubits, with no past outcome record carried across calls, and achieves \(4Ns/(1-s)-O_s(\sqrt N)\) QFI. The slope itself is the known channel-extension coefficient, not new. The deficit estimate is an achievable bound, not a claim that \(\sqrt N\) is unavoidable or optimal.

**Priority and exact scope.** [Liu et al. 2023](https://doi.org/10.1103/PhysRevLett.130.070803) and [Liu et al. 2024](https://doi.org/10.1002/qute.202400094) give general finite comb SDPs and numerical independent amplitude-damping examples. [Kurdziałek et al. 2023](https://arxiv.org/abs/2212.08106) give iterative adaptive upper bounds and asymptotic equivalence with the optimal parallel coefficient; that recurrence is a bound rather than the attained Bellman identity below. [Kurdziałek et al. 2025](https://arxiv.org/html/2403.04854v2#S5.SS4), §V.4 and Appendix D.4, numerically optimize this same axial noise with restricted ancillary size and display simpler two-call code-sector/conditional-preparation strategies. Their exact recursive protocol in §V.3/D.3 is for **perpendicular** signal/noise, a different channel. [Li–Tong–Bao](https://arxiv.org/html/2610.31249v2) already exploit heralded reset for continuous-time erasure sensing. [Compressing Quantum Fisher Information](https://arxiv.org/abs/2602.09358) treats pure-state phase QFI compression including a sequential example; no novelty attributed to retaining a pure local state's tangent span. [Šafránek](https://doi.org/10.1103/PhysRevA.95.052320) distinguishes SLD from continuous Bures information at rank changes.

For bounded memory, [Liu–Yang, Quantum 8,1571 (2024), §3.2](https://doi.org/10.22331/q-2024-12-18-1571) report numerical near-linear scaling with one ancilla and expressly leave the analytical SQL question for one ancilla open. Our two-ancilla upper construction does not settle that **one**-ancilla question; it does prove a fixed-size-memory result with the exact coefficient. [Zhou, PRL 133,170801 (2024)](https://arxiv.org/abs/2402.18765) studies restricted controls (e.g. no noiseless ancilla, unitary/unital conditions) rather than the ideal tag and nonunital CPTP alignment/reset used here. [Zhou–Manes–Jiang, PRA 109,042406 (2024)](https://arxiv.org/abs/2303.00881) attain metrological limits for continuous-time Markovian Hamiltonian sensing via error correction on many noisy probes and short recovery intervals, not fixed noninfinitesimal black-box slots with three reusable qubits. Their Appendix A, Eq. (A10), explicitly discusses the same axial channel in the separate setting of parallel probes, quotes the known assisted asymptotic coefficient, and explains why continuous frequent recovery is a different access model. [Liu et al. 2026, §V](https://arxiv.org/html/2608.01114v1#S5) provide finite-memory numerical approximations for simultaneous multi-parameter magnetometry with amplitude damping, including step-independent repeated controls; their iteration is not our exact scalar recursion or fixed-slot slope certificate. [Demkowicz-Dobrzański et al. 2012](https://arxiv.org/abs/1201.1693) already give the channel-extension rate for amplitude damping, and [Zhou–Jiang 2021](https://doi.org/10.1103/PRXQuantum.2.010343) prove general asymptotic attainment with resource assumptions different from this bounded-register construction. The claim is exact ordered scalar recursion, strict non-saturation with future uses, and constructive bounded-memory saturation of the already known coefficient; not general adaptive priority or hardware advantage.

## 1. Exact flagged recursion and physical implementation

Represent an incoming normalized pure local family by \(\psi,\ \psi'=i\sqrt q\,\phi\) at \(\theta=0\), with \(\phi\perp\psi\), \(\|\phi\|=1\). A scalar phase gauge removed the component along \(\psi\); for \(q=0\) orientation is irrelevant. Let \(F_{\mathcal T}(q)\) be the suffix's best local QFI divided by four for this incoming jet. It is independent of the input dimension and orientation because prior free phase-independent controls can isometrically rotate the nominal pair. Define \(F_\varnothing(q)=q\). Locally identical density jets (nominal density and first derivative) can always substitute for one another; later parameter-dependent channels propagate equality of both through the product rule. No claim about representing the full family away from zero is needed.

For an upper bound purify initial states and intermediate CPTP operations, keep their ideal auxiliary outputs, and append ideal classical labels of each \(B_j\), with phase-independent controls allowed to consult them. Ignoring added registers reproduces any original protocol, hence cannot give a larger optimal QFI. Conditional on channel histories the enlarged state's branches are pure. At the next stage, all dependence on state placement is captured by the positive contraction
\[
 M=\begin{pmatrix}x&c\\\bar c&y\end{pmatrix},\qquad 0\preceq M\preceq I,
 \qquad c=\langle\psi|P\otimes I|\phi\rangle=z+iw .
\]
Namely only the Gram entries of the excitation projection on the nominal pair appear. Any \(M\) is realized as a compression of this projection by the isometry
\[
 V u=|0\rangle_S\otimes \sqrt{I-M}\,u
               +|1\rangle_S\otimes\sqrt M\,u ,
\]
for a two-dimensional logical \(u\) and one qubit of temporary memory. Coherently tag the sensor basis (append \(A=0\), apply controlled copy \(0\to0,1\to1\)). The no-loss output has even \(S,A\) parity and the loss output odd parity exactly, including tangents. This tag does not reveal individual basis amplitudes within no loss. Therefore one can recover the ideal *classical* Kraus labels in the actual sensor circuit with no access to its environment. Each positive-probability output branch can recanonicalize state and horizontal tangent into one logical qubit by a fixed nominal isometry and discard the remaining registers; there is no first-jet loss because both span vectors are retained. Then reuse the sensor, temporary memory and tag. A classical history or equivalent controller state may be retained for the general exact finite optimum.

Put \(\ell=1-s,\ W=1-\ell x,\ D=\sqrt q\), and use weights \(w_0=W,w_1=\ell x\) (distinguished from \(w=\Im c\)). The classical derivatives are \(w_0'=2\ell D w=-w_1'\), sign for the chosen \(iD\phi\) convention. With \(v_j\) the unnormalized output and \(\nu_j=v_j'\), the inner products and squared norms are
\[
\begin{array}{c|cc}
 j&\langle v_j|\nu_j\rangle&\|\nu_j\|^2\\ \hline
 0&i(sx-\ell D c)&sx+2sD z+q(1-\ell y)\\
 1&i\ell D c&\ell q y.
\end{array}
\]
When \(w_j>0\) the output horizontal squared tangent is
\[
 q_j=\|\nu_j\|^2/w_j-|\langle v_j|\nu_j\rangle|^2/w_j^2.
\]
The exact recursion is
\[
 F_{s,\mathcal T}(q)=
 \sup_{0\preceq M\preceq I}\sum_{j:w_j>0}
       \left[\frac{(w_j')^2}{4w_j}+w_j F_{\mathcal T}(q_j)\right]. \tag{1}
\]
For a direct sum classical-plus-quantum block the QFI is the probability-weighted conditional SLD QFI plus the probabilities' classical FI. Applying this recursively is valid also when histories have nonzero score, since each next conditional classical score has conditional mean zero, or equivalently one can use the direct-sum identity at every stage. A nominal zero-mass branch has zero density derivative, even if the ket derivative is nonzero, so it and its descendants contribute nothing under pointwise SLD; finite subsequent differentiable channels do not change that first jet. The upper and physically flagged lower match by recursion (first as suprema; attainment established below). The compression is valid for the **local** pure-state tangent, not a rank-changing continuous Bures extension.

## 2. Shape of the continuation and elimination of controls

Every \(F_{\mathcal T}\) is concave as a function of \(q\). For a convex combination \(q=\sum_b\lambda_b q_b\), form the nominal pair using an orthogonal extra noiseless label \(b\) with state \(\sum_b\sqrt{\lambda_b}|b\rangle\psi_b\) and horizontal derivative \(i\sum_b\sqrt{\lambda_b q_b}|b\rangle\phi_b\). Their inner product is zero and the derivative squared norm is \(q\); this is a valid fixed re-encoding of the incoming pair. Reading the label and applying a suffix strategy for each branch gives the weighted QFI. Parking the incoming tangent and ignoring all calls gives \(F_{\mathcal T}(q)\geq q\). Consequently the concave function \(H(q)=F_{\mathcal T}(q)-q\) is nondecreasing on \(q\geq0\): a negative secant slope would force a linear decrease to negative values at large \(q\). In particular \(F_{\mathcal T}\) increases by at least an increment in \(q\).

For continuity also at zero, let \(d\in\mathbb R\) be the signed input horizontal amplitude. For any fixed control the output nominal state is independent of \(d\) and its density derivative affine in \(d\), so the fixed-state SLD quadratic form is convex in \(d\). The supremum is convex, even by replacing the logical sign of \(\phi\), and finite. Finiteness follows also by coherently dilating all controls/noise and centering each phase generator with operator norm \(1/2\), giving \(F_{\mathcal T}(d^2)\leq(|d|+|\mathcal T|/2)^2\) by the derivative triangle bound. A finite convex function on the real line is continuous; composition with \(\sqrt q\) gives endpoint continuity.

Suppose now \(0<s<1\) and \(0<x<1\). The total immediate QFI divided by four is
\[
 T=\sum_j (w_j')^2/(4w_j)+w_j q_j
   =q+\frac{s x(1-x)+2sD z-\ell q z^2/x}{W}. \tag{2}
\]
For fixed \(x,y,z\), reducing \(|\Im c|\) to zero increases both \(q_j\) (the squared overlaps in their expressions decrease) and remains feasible. The immediate \(T\) does not change. The objective in (1) is \(T+\sum_j w_jH(q_j)\); since \(H\) is nondecreasing, a real control is always at least as good. For real \(c=z<0\), flipping its sign increases \(q_0\) and preserves \(q_1\) and the weights; feasibility is unchanged. Thus an optimum can always use real \(z\geq0\), and its conditional flag score is zero. This assertion is not that arbitrary admissible controls have zero classical Fisher information.

For real \(c=z\), feasibility is
\[
 0\leq z\leq e=\sqrt{x(1-x)},\qquad z^2/x\leq y\leq1-z^2/(1-x).
\]
As \(y\) increases through that interval, \(q_1=q(y/x-z^2/x^2)\) fills \([0,R]\), with
\[
 R=\frac{q[x(1-x)-z^2]}{x^2(1-x)},\qquad q_0=(T-\ell x q_1)/W.
\]
The weighted mean \(T=Wq_0+\ell xq_1\) is fixed. By concavity of \(F_{\mathcal T}\), the best split moves toward equal branch \(q\)'s: if \(R\geq T\), choose \(q_1=q_0=T\); otherwise choose \(q_1=R\). For proof, when \(q_1<T\), increasing \(q_1\) toward \(T\) and decreasing \(q_0\) with the fixed weighted mean is a convex mixture with the balanced pair, whose value is at least as high by Jensen. Thus (1) is the explicit two-real-variable maximization
\[
 F_{s,\mathcal T}(q)=\max_{x,z}\left[
    W F_{\mathcal T}\!\left(\frac{T-\ell x\,b}{W}\right)+\ell x F_{\mathcal T}(b)\right],
 \quad b=\min(T,R),\quad 0<x<1,\quad 0\leq z\leq e, \tag{3}
\]
with the endpoint parking value \(F_{\mathcal T}(q)\) admitted. At \(x=0\) or \(1\), all feasible controls have objective at most that value: with \(x=0\) the loss term is absent and remaining tangent is at most \(q\); with \(x=1\) the mean is \(q\) and Jensen applies. They attain it by parking. As \(x\to0\) or \(1\), uniformly in \(z\), (2) gives \(\limsup T\leq q\); Jensen bounds (3) by \(F_{\mathcal T}(T)\), and continuity handles any rank-boundary limit. Thus no larger value lives solely in such a limit. In fact \(x\in(0,1),M=xI\) gives \(q_1=q\), \(q_0=q+sx(1-x)/W^2>q\), so an optimum lies strictly in the interior of the \(x\) interval. Compactness away from endpoints and continuity give attainment. For \(s=0\) the use is a parameter-independent reset, so the value is \(F_{\mathcal T}(q)\). For \(s=1\), tangent norm grows by at most \(1/2\) and that sum is attainable by aligned balanced sensing, giving \(F_{\mathcal T}((\sqrt q+1/2)^2)\).

## 3. Why maximal overlap is wrong before a later active call

If the suffix \(\mathcal T\) contains any positive-survival call, skip initial zero-strength slots and write its first active survival as \(u>0\). At \(q=0\) choose its optimal \(x\in(0,1)\) (interiority also for \(u=1\)); one has \(q_0=u e^2/W^2,\ q_1=0\). Keep that \(x\) as \(q\) rises, and at that active slot choose saturation \(z=e,\ y=1-x\). Then \(q_0=u(\sqrt q+e)^2/W^2,\ q_1=0\), with unchanged weights. The next suffix has \(F(q')-q'\) nondecreasing, so this feasible choice shows
\[
 F_{\mathcal T}(q)\ \geq\ F_{\mathcal T}(0)+ (2u e/W)\sqrt q .
\]
The constant is strictly positive. Thus its right slope at zero is infinite, quantitatively, including when the remaining active slot is noiseless.

At an intermediate \(0<s<1,\ q>0\), fix interior \(x\) and start from saturated \(z=e\) with \(y=1-x\). The loss has \(q_1=0\) and no loss has \(q_0>0\). Perturb \(z=e-\epsilon\) and choose \(y=1-z^2/(1-x)\). The feasible loss tangent now has \(q_1=R=\Theta(\epsilon)>0\); the no-loss tangent changes by only \(O(\epsilon)\) at fixed \(x\). The active tail increases the loss term by \(\Omega(\sqrt\epsilon)\) with weight \(\ell x>0\), whereas the no-loss cost is \(O(\epsilon)\) since a finite concave function is locally Lipschitz at positive \(q_0\). This strictly beats every saturated control at an interior \(x\); the endpoints cannot maximize by §2. Moreover \(z=0\) is not optimal for \(q>0\): the balanced split at \(z=0\) has both branch tangents positive, since \(T>q\), \(R=q/x>0\). For small \(z>0\), \(T\) increases to first order and \(R\) decreases only quadratically. Keep the old loss tangent or decrease it only enough to stay within the new capacity. The no-loss tangent then strictly increases linearly; its gain is at least this increase times its weight, while the loss cost is quadratic by local Lipschitz continuity. Thus \(0<z<e\) at an optimal real control, and (3) gives positive tangent in both branches. This is an optimum obstruction, not just a numerical example or a test of one special final design. It applies at the middle step of three genuinely sensing calls when the first no-loss branch brings in positive tangent.

**Exact three-call comparison.** On ordered survivals \(1,1/4,1\), any first unitary can give at most \(q=1/4\) on a purified input, and this is attained. The last suffix is \(F(q)=q+\sqrt q+1/4\). At a saturated middle slot its continuation value is
\[
 T+\sqrt{s}(\sqrt q+e)+1/4,\qquad s=1/4,
\]
because the loss tangent is zero. Completing the square in (2) gives \(T-q\leq sx/(1-s)\leq 1/3\), while \(e\leq1/2\). Thus the continuation is at most \(q+\sqrt q/2+5/6\) for each conditioned pure branch, also at the parking endpoints. This comparison allows extra filtering before the middle slot: if additional accessible labels have QFI contribution \(C\) in quarter units and pure conditional tangents \(q_b\) with weights \(p_b\), then \(C+\sum_b p_bq_b\leq1/4\) by data processing from the first noiseless call. Jensen gives \(\sum_b p_b\sqrt {q_b}\leq 1/2\). In total the quarter-QFI is at most \(1/4+1/4+5/6=4/3\), or QFI \(16/3\), for a policy saturated on each such conditioned middle input. For an allowed nonsaturated compression \(x=3/5,z=1/4,y=27/32\) and \(q=1/4\), the weights are \(11/20,9/20\), \(q_0=8407/15488,\ q_1=355/1152,\ T=1539/3520\). Exact square bounds \(\sqrt{q_0}>736/1000,\ \sqrt{q_1}>555/1000\) prove its QFI exceeds \(590377/110000\), which is \(11131/330000>0.0337\) above the uniform saturated upper bound. Actual QFI of this witness is about \(5.36894\); no claim that it is the unrestricted optimum.

## 4. Known asymptotic slope via deterministic bounded-memory tangent flow

This section uses the same unrestricted ideal inter-slot CPTP controls, but does **not** require a growing classical outcome history. Let all \(N\) survivals equal \(0<s<1\), write \(\kappa=s/(1-s)\). From (2) for \(q>0\),
\[
 T-q=\kappa x-\frac{\ell q}{xW}
       \left(z-\frac{\kappa x}{\sqrt q}\right)^2\leq\kappa .
\]
At \(q=0\), \(T=sx(1-x)/W\leq\kappa\); endpoints are already covered. By induction in (1) or (3), \(F_N(0)\leq N\kappa\), the familiar rate in this normalization.

**Deterministic seeding through one call.** Start in a tagged balanced sensor superposition. Align normalized no-loss nominal state and horizontal tangent to the logical \(|0\rangle,|1\rangle\); align the loss state's nominal vector to \(|0\rangle\), with zero derivative there. Erase the tag/flag after correction, with no selected outcome discarded from the data. At nominal phase the resulting density is pure, and its first derivative is that of a pure logical tangent of amplitude \(W \sqrt s/(1+s)=\sqrt s/2\), where \(W=(1+s)/2\); so \(q_\ast=s/4>0\). The actual family away from the nominal point can be mixed; equality of first jets is all that the next slot or SLD readout needs.

**Balanced increments without a stored history.** With real control, if \(R\geq T\), choose \(y\) so that both post-call tangents equal \(T\). In either parity apply a conditional phase-independent state-and-tangent alignment to the same nominal logical pair. These two density first jets coincide and their weights have zero first derivative. Thus the unconditional recovered jet is exactly the same pure-state jet, and erasing the flag costs no pointwise SLD. This can be one CPTP map with parity-conditioned Kraus operators instead of a persistent classical register. Reuse the three qubits. The control at each slot depends only on the known schedule's current \(q\), not a realized record. Measure/reset implicitly uses an external transient entropy sink; no claim of a globally closed unitary circuit on three qubits.

For a deterministic warmup at \(q\geq q_\ast\), use \(z=0\), a fixed small \(x\in(0,1/2]\) with \(2s x^2\leq q_\ast\), and balance. Indeed \(R=q/x\geq T=q+sx(1-x)/W\) iff \(q\geq s x^2/W\), satisfied since \(W\geq1/2\). Constant positive increments reach any fixed threshold in finitely many steps depending only on \(s\).

For large \(q\), put \(C=2\kappa\), choose \(x=1-C/\sqrt q\), \(z=\kappa x/\sqrt q\). Whenever \(q\geq\max\{1,4C^2\}\) one has \(x\in[1/2,1)\), \(T=q+\kappa x\), and
\[
 R=q/x-\kappa^2/(1-x)\ \geq\ T .
\]
Indeed \(q(1-x)/x\geq C\sqrt q\), whereas \(\kappa^2/(1-x)+\kappa x\leq(\kappa^2/C)\sqrt q+\kappa\leq C\sqrt q\). Also \(R>0\) proves that \(z^2<x(1-x)\), so the chosen compression is feasible and balancing is allowed. Each step increments by at least \(\kappa/2\), so the accrued shortfall \(\sum \kappa C/\sqrt{q_i}\) is \(O_s(\sqrt N)\) by an elementary sum along \(q_i\geq q_{\rm warm}+i\kappa/2\). The final Pauli-\(Y\) measurement on the carrier has nominal probabilities \(1/2,1/2\) and derivatives \(\sqrt q,-\sqrt q\), yielding the achieved QFI and requiring no nominal zero-probability outcome. This is a constructive memory bound with a known leading coefficient and an explicit achievable sublinear deficit, not a statement about the optimal deficit or a separation from asymptotically optimal parallel strategies.

**Uniform heterogeneous variant.** The same bound works for a fixed ordered sequence all of whose survivals lie in a fixed compact strip \(0<a\leq s_i\leq b<1\). The upper accumulates \(4\sum_i\kappa_i\), \(\kappa_i=s_i/(1-s_i)\), by the single-slot identity. The first call deterministically seeds \(q_\ast=s_1/4\geq a/4\). Choose one warmup \(x\in(0,1/2]\) with \(2b x^2\leq a/4\); the warmup increments are at least \(a x(1-x)\). Use the uniform threshold \(\max\{1,16[b/(1-b)]^2\}\), thereafter choose \(C_i=2\kappa_i\) and the same balanced update with that slot's survival. The per-slot gains are at least \(a/[2(1-a)]\), and losses from the corresponding upper increments are at most \(2[b/(1-b)]^2/\sqrt q\). Thus the three-qubit, no-record construction attains \(4\sum_i\kappa_i-O_{a,b}(\sqrt N)\) uniformly over the given order. This does not characterize finite-\(N\) order preference or sequences approaching the noiseless endpoint.

## Replay

Run python verify.py for exact symbolic identities, the rational witness and direct sensor–memory–tag physical controls; see replay.txt and manifest.sha256. A separate comparison, python compare_comb.py --fetch, evaluates the Bellman formula for three identical calls at decay probabilities 0.10, 0.40 and 0.80 and matches the numerical author comb data to \(9.2\times10^{-11}\) absolute or better at these points; see comb-replay.txt and the pinned original source URL and hash in the script. This is a finite numerical smoke test, not an independent certificate of global numerical optimization or the mixed-strength witness. The main proof does not rely on the numerical optimizer or on any claim about globally preferred order over all finite strengths.
