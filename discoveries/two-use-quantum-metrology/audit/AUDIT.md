# Independent audit: two-use phase metrology with amplitude damping

**Assessment.** The theorem is correct in the stated fixed-call local-QFI model. The analytic strictness across the two-survival square, clean boundaries, and one-syndrome construction are the contribution. The equal-decay hierarchy, passive equal-decay objective and threshold, and the broad idea of heralded-loss recycling have prior antecedents. No first hierarchy, global adaptive optimum, experimental advantage or finite-shot precision claim.

## Operational model, flags, and exact resource match

For each known survival $s_i$, the named channel is ordinary qubit amplitude damping followed by $e^{i\theta |1\rangle\langle1|}$, with a common unknown phase, fresh inaccessible environments, and one call to **each** of two named channels. These are not two freely selectable calls to the stronger channel. If up to one per name is allowed, padding with an unused ground-state input produces the same maximal QFI. Parallel operations may use arbitrary initial joint ancillas and joint readout, but may not feed the first output back into the other input. Sequential operations compare predetermined orders; no indefinite order, free additional phase encoding, energy/sensing-time or practical latency comparison.

For an arbitrary purified input $\sum_x|x\rangle|a_x\rangle$ set $p_x=\|a_x\|^2$. Purification only enlarges the admissible parallel ancilla, since a partial trace cannot increase QFI. Let $L$ be a subset of lost bits; Kraus weight at $x$ is
$$
 w_L(x)=\mathbf1_{\{L\subseteq x\}}\prod_{i\in L}(1-s_i)\prod_{i\in x\setminus L}s_i .
$$
In its pure sector the images of distinct $x\supseteq L$ are orthogonal already on the sensing bits, regardless of overlap of the $|a_x\rangle$. Their phase labels differ as $K=|x|$ since subtracting $|L|$ is a branchwise common phase. The sector weights are phase independent. Adding the unattained environmental label gives an upper bound $4\sum_L P_L\operatorname{Var}(K\mid L)$ by trace monotonicity. It does **not** assume experimental access to the environment.

For every distribution $p$, implement the upper bound with the coherent tagged state $\sum_x\sqrt{p_x}|x\rangle|x\rangle$. The jointly accessible basis pairs $(x\setminus L,x)$ for distinct $(L,x)$ are mutually orthogonal. Their bitwise parity is $L$; it can be read by a nondemolition measurement of parity only. Thus jumping versus not jumping is identified with sensing/tag output data and no loss of within-sector coherence. This proves the exact parallel variational objective against arbitrary ancillas. Its form as a sum of infima over posterior means gives concavity; continuity at a zero-weight branch follows from the bounded value of $K\in[0,2]$. In the same-strength case averaging the laws under swap is admissible as one coherent tagged probe, not a classical mixture that might lose QFI.

For an explicit local measurement in any nonzero-weight pure sector, write its normalized state at the chosen phase as $|\psi\rangle$, generator $K$ (up to scalar), variance $\sigma^2>0$, and normalized tangent $|e\rangle=i(K-\langle K\rangle)|\psi\rangle/\sigma$. Project onto $(|\psi\rangle\pm|e\rangle)/\sqrt2$ and the complement, after resolving the orthogonal parity sectors. The two main outcomes have probabilities $P_L/2$, slopes $\pm P_L\sigma$, and total FI $4P_L\sigma^2$. A zero-variance sector contributes zero. This confirms a local achievable information calculation; a globally optimal phase-independent measurement or estimator is not claimed. The feed-forward control itself is phase independent.

The first jump freezes its own sensor/tag to $|0,1\rangle$, and the other sensor/tag is a known two-level pure state independent of phase. If phase were encoded before damping, the first-jump phase would only be global. Reset only this remaining pair to its optimum for its own last call; because the pair's old state and target are known, a conditional unitary on the two-dimensional tagged code is sufficient. The first no-jump branch undergoes exactly the parallel continuation, and the initial parity observation leaves that branch unchanged. Both branches have the second call once, no runs are discarded, and first flag probabilities do not depend on phase.

## Strictness and active boundaries

Write $p=(a,b,c,d)$ in order $00,10,01,11$ and $s_1=r^2,s_2=t^2$. One use with excitation $y$ gives $Q_s(y)=4s y(1-y)/(1-y+sy)$. For $s>0$, completing the square gives
$$
 Q_s^*-Q_s(y)=\frac{4s(y-(1+\sqrt s)^{-1})^2}{1-y+sy}.
$$
Thus the reset adds the first-jump probability times that strict deficit, with input $y=d/(b+d)$ when channel 1 is first and $y=d/(c+d)$ when channel 2 is first. If both $s_i$ are interior and $d=0$, a usable parallel optimum has $b>0$ or $c>0$; that corresponding first-jump branch has $y=0$ and improves. If $d>0$ and one of $b,c$ is zero, the corresponding $y=1$ improves. All remaining no-gain candidates therefore have $d,b,c>0$ and
$$
 b=td,\qquad c=rd .
$$
All four loss patterns then have positive mass, including for $a=0$. Their posterior mean losses are $h_1=1+t/(1+t)$, $h_2=1+r/(1+r)$. If $g_x=\sum_Lw_L(x)(|x|-\mu_L)^2$ at the actual means, $g_x$ is the partial derivative of the exact QFI divided by four. Transfer of mass from an active to any other coordinate gives equal active scores and no larger inactive one. This uses differentiability only where the relevant pattern probability is positive. In the one-noiseless case a loss of the noiseless bit is identically impossible and its term can be removed before differentiating. No KKT claim at a vanishing but potentially active pattern is needed.

Let $m$ be the no-loss mean. When both hypothetical jump continuations already optimize, the alternating score sum cancels the single loss terms:
$$
 D(m)=g_{11}-g_{10}-g_{01}+g_{00}
      =s_1s_2(2-m)^2-(s_1+s_2)(1-m)^2+m^2 .
$$
If $a=0$, $m>1$ from $d>0$, and $D(m)>0$ since $D(1)=1+s_1s_2>0$ and both its quadratic and linear coefficients are positive for interior $s_i$. But $b,c,d$ are active and $a$ inactive, so $D(m)\le0$.

If $a>0$, equality of the scores of $00,10,01$ gives
$$
 m^2=\frac{t^2}{(1+t)^2}+\frac{r^2}{1-r^2}(1-2m)
     =\frac{r^2}{(1+r)^2}+\frac{t^2}{1-t^2}(1-2m).
$$
For $r\ne t$, subtract and divide by the nonzero difference of $r^2/(1-r^2)$ and $t^2/(1-t^2)$. This forces $z=1-2m=(1-r)(1-t)(r+t+2rt)/[(1+r)(1+t)(r+t)]$. Substitution gives a positive residual
$$
 m^2-\frac{t^2}{(1+t)^2}-\frac{r^2z}{1-r^2}
 =\frac{r^2t^2(1-r)(1-t)}{(1+r)(1+t)(r+t)^2}>0.
$$
For a direct algebra check, put $S=r+t,P=rt,C=(1+r)(1+t)$. Then $m=(S^2+SP-P-P^2)/(SC)$; clearing $S^2C^2$ in the residual gives $(S^2+SP-P-P^2)^2-t^2(1+r)^2S^2-r^2(1-t^2)S(S+2P)=P^2(1-r^2)(1-t^2)$.

For $r=t=\sqrt s$ with $a>0$, instead use $D(m)=[(1-s)m+2s]^2-2s=0$, so its bracket is $\sqrt{2s}$. Equality $g_{00}=g_{10}$ gives $m[(1-s)m+2s]=2s/(1+r)$. Substitution forces $1-\sqrt{2s}=1-r$, contradicting $r>0$. Taking an unequal-strength limit would not have established strictness, so the separate diagonal argument matters.

If channel 1 is lossy with $0<s_1<1$ and channel 2 is noiseless, use channel 1 first. At a parallel optimum it must have $b+d>0$: otherwise the objective is at most 1, while independent one-use optimal probes give $1+Q_{s_1}^*>1$. The only possible equality after the intervention would be $b=d>0$, so $h_1=3/2$. Equal active scores $g_{10}=g_{11}$ force $m=3/2$. Then these scores equal $1/4$, but $g_{00}=9/4$; contradiction whether $a$ is active or inactive. Relabel for the other edge. If either survival is zero, its call is a phase-independent replacer and can be absorbed into pre- or postprocessing of the one informative call. If both survivals are one, a dilation of phase-independent controls gives an effective generator which is a sum of two conjugate projectors; its width is at most two, giving QFI at most four even before tracing out dilation registers. Parallel GHZ attains four.

## Independent exact witness with opposite order preference

Use $s_1=4/25,s_2=9/25$ and $p=(33,41,19,107)/200$. These are unequal, strong damping and are not on the 2024 ray $1-s_1=2(1-s_2)$ or its swap. Trial posterior means $m=489/1000,h_1=1485/1000,h_2=1477/1000$ give the globally valid parallel upper bound $U=149825737/156250000$ by the maximum of the four scores. The chosen probe's unmodified parallel QFI is $1173318793/1224558300$. Using the first-jump parity and optimal reset for the second channel gives

| Predetermined order on this one prepared input | QFI |
|---|---:|
| More noisy first ($1\to2$, reset $y=5/8$) | $3360283141/3449460000$ |
| Less noisy first ($2\to1$, reset $y=5/7$) | $9944872779/10204652500$ |

The second minus the first is positive, $101407/256452000$. The second exceeds $U$ by $9986716030943/637790781250000>0.01565$. Both exceed $U$. This reverses the noisier-first preference of the specified GHZ-based procedure in its high-survival region, **only as a fixed-input constructive comparison**, not as a ranking of the two globally optimal sequential orders.

The integer-rational replay directly propagates weights and integer phase labels along each Kraus path. It checks that nonzero terms occupy disjoint physical sensor/tag basis pairs with parity labels, reconstructs both adaptive interventions before the second call, and evaluates the two-outcome local measurement above via pairwise mass products, not the closed gain equation. It also verifies the author's diagonal/unequal GHZ and high-loss diagonal certificates, a new noiseless/damped-edge exact control (parallel $144/49$, lossy-first $625/196$), the zero and noiseless equality controls, and the exact diagonal untagged-GHZ negative control $27/14<9/4$. The latter follows by direct 00/11 density-block calculation: the double-loss population overlaps the untagged informative block, making its QFI $16|\rho_{00,11}|^2/(\rho_{00,00}+\rho_{11,11})$.

    python3 -B verify.py | diff -u replay.txt -
    sha256sum -c SHA256SUMS

## Priority checks and claim boundary

- [Liu–Hu–Yuan–Yang, 2023 PRL, Fig. 2 and Supplement E/I](https://arxiv.org/html/2203.09758) explicitly state the strict parallel–sequential–indefinite hierarchy for all equal nontrivial losses. Their SDPs are exact characterizations; the displayed amplitude-damping hierarchy relies on numerical evaluations/circuit synthesis. The full claim is not a first claim on the diagonal or at survival $0.6$.
- [Liu et al., 2024 §5.1.3, Fig. 8](https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/qute.202400094) already numerically treat the same *nonidentical* channels with known decays $p_1=2p_2$, $p_i=1-s_i$, and choose the better of the two sequential orders. Do not claim first unequal-strength examples, nor claim they covered the full two-dimensional parameter domain. Their §5.1.1 has the equal-strength benchmarks too.
- [Demkowicz-Dobrzański et al., 2009 §§IVB, VA, VI and VII B](https://arxiv.org/pdf/0904.0456) already give the passive one-lossy-arm variance formula, the time-bin binary-loss-pattern objective (Eqs. 55–56, exact if the other arm is lossless), the swap/concavity reduction (Eq. 57), and the local sector SLD measurement (Eq. 50). Eq. 47 at $N=2$ has exactly the bracket $1+2s-7s^2$; its positive zero is the diagonal threshold. The optical time bins have common transmissivity, and this is not the heterogeneous adaptive comparison. Orthogonal loss sectors in that passive example correspond physically to loss of a known time-bin photon; use no novelty claim on the broad flag idea.
- [Li–Tong–Bao, November 2026 v2, Heralded recycling, Fig. 4 and Eq. (13)](https://arxiv.org/html/2610.31249v2) already use conditional re-preparation of phase-independent heralded losses to improve frequency QFI. Their continuously operating three-level erasure model loses *both* sensing levels to an exterior state and compares to no-reset Ramsey, not two black-box qubit amplitude-damping calls with arbitrary parallel entanglement and ancillas. Thus heralded recycling as a broad mechanism is explicitly prior. The exact bounded-call theorem and particular tag/parity realization are narrower.
- [Kurdziałek et al., 2023 PRL](https://arxiv.org/abs/2212.08106) concern asymptotic equivalence in large call count. [Meyer et al., 2025 finite-sample metrology](https://arxiv.org/abs/2307.06370) emphasize different finite-sample success criteria; no success-probability or global phase-estimation result follows here. This is also not an error-correction protocol that restores a lost first-call phase: it repurposes a still-unused second call, not lost signal.

No proof of the entire two-survival claim, with this two-call ancilla-inclusive benchmark and exact strictness argument, was located in this bounded check. A carefully attributed analytic result is supportable. It is not a generic metrology paradigm discovery.
