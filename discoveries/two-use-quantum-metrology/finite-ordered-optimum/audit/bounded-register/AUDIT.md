# Fixed-register asymptotic optimality for amplitude-damped phase calls

**Independent audit: pass.** Fix a survival $s\in(0,1)$. For $N$ independent calls to a qubit channel consisting of amplitude damping of survival $s$ and phase $e^{i\theta P}$, $P=|1\rangle\langle1|$, a single routed sensor plus two noiseless qubits suffices to attain
\[
 0\le \frac{4s}{1-s}N-F_N\le O_s(\sqrt N)
\]
for the final pointwise SLD QFI. Only a known nominal phase is used in designing the phase-independent controls; one may set it to zero. These are ideal general CPTP, not unitary-only, controls between fixed-duration calls. No outcomes, probes, or damping environments are retained across calls other than the three-qubit register. No shot is rejected. Constants need not be uniform as $s$ tends to 0 or 1. No global-phase or fixed-control synthesis claim.

## Achieving the slope with one carrier

The state carried between uses needs only the following nominal density and derivative on a known two-dimensional subspace:
\[
 \rho_0=|0\rangle\langle0|,\qquad
 \dot\rho_0=i d(|1\rangle\langle0|-|0\rangle\langle1|),\quad d\ge0.
\]
One cannot assume the nearby family is exactly pure. That does not matter for the theorem: the density and derivative evolve exactly by the product rule under each differentiable phase slot and each fixed control. The final Pauli-$Y$ readout on the carrier has probabilities $1/2,1/2$ and derivatives $\pm d$, so its pointwise classical FI is already $4d^2$, identical to the SLD QFI. This readout itself has no zero-mass ambiguity.

Let $q=d^2$, $a=s/(1-s)$, $\ell=1-s$. As in the exact finite control reduction, realize a real symmetric excitation compression on the local state/tangent pair by the isometry
\[
 W_M=|0,0\rangle\otimes\sqrt{I-M}
       +|1,1\rangle\otimes\sqrt M,\qquad
 M=\begin{pmatrix}x&z\\z&y\end{pmatrix},\qquad 0\preceq M\preceq I .
\]
The first two factors are the routed sensor and an accessible parity tag; the last is a logical qubit. Extend $W_M$ from the carrier to a unitary on the three-qubit register. The even and odd accessible parity sectors contain, respectively, the no-jump and jump nominal state and derivatives. For $0<x<1$ both masses $p_0=W=1-\ell x$, $p_1=\ell x$ are positive. Real $M$ makes both mass derivatives zero. Direct Kraus differentiation gives the conditional squared tangents $D_0,D_1$ with
\[
 T=W D_0+\ell xD_1=q+\frac{sx(1-x)+2sdz-\ell qz^2/x}{W},\qquad
 D_1=q(y/x-z^2/x^2).
\]
For example, the no-jump squared derivative norm before normalization is $sx+q(1-\ell y)+2sdz$, its state/derivative inner product is $i(sx-\ell dz)$, and the corresponding jump quantities are $\ell qy$ and $i\ell dz$. Subtracting the latter overlaps after normalization gives the formulas.

Whenever $z^2\le x(1-x)$, physical feasibility is precisely $y\in[L,U]=[z^2/x,1-z^2/(1-x)]$. Over this interval, $D_1$ increases from zero to
\[
 D_1^{\max}=q/x-qz^2/[x^2(1-x)].
\]
Thus if $q>0$ and $D_1^{\max}\ge T$, choose $y=L+xT/q$ to make both conditional squared tangents equal to $T$. Independently rotate each branch's normalized state and horizontal tangent to the same two canonical carrier directions, then forget parity. Concretely, if $\Pi_j$ are accessible parity projections and each $U_j$ is the chosen fixed unitary extension, $R(\rho)=\sum_j U_j\Pi_j\rho\Pi_j U_j^\dagger$ is a CPTP map on the same register. Two local pairs are aligned with the same sign; each branch weight has zero derivative. Hence this unconditioned map returns the nominal carrier density and $d_{\rm new}=\sqrt T$. The two branch operations are internal to one fixed control: only ephemeral syndrome processing or a fixed local reset bath is needed, not a stored string of outcomes. This is not a wholly closed three-qubit unitary implementation. One can precompute the deterministic control schedule from $s$ and call index; no claim is made about memory for specifying analog settings or their digital precision.

**Seed.** Input the balanced tagged state $(|0,0\rangle+|1,1\rangle)/\sqrt2$. The raw no-jump weight is $P_0=(1+s)/2$, its horizontal norm is $\sqrt s/(1+s)$, and the jump branch has zero phase derivative. Align the two nominal states to the same $|0\rangle$, preserving the no-jump derivative toward $|1\rangle$, and erase parity. By linearity the output derivative has norm $P_0\sqrt s/(1+s)=\sqrt s/2$. Thus $q=s/4>0$.

**Finite startup.** Set $x=1/4,z=0$. Then $T=q+\Delta$, with $\Delta=3s/[4(3+s)]>0$. The proposed $y=xT/q$ is at most one exactly when $q\ge sx^2/W$. Here $W\ge3/4$, so $sx^2/W\le s/12<s/4$, and the seed and all later startup points work. Repeat finitely until $r=\sqrt q\ge \max(4a,1)$. For fixed $s$ this time is independent of the final horizon, but can be large.

**Asymptotic steps.** Take $x=1-2a/r\in[1/2,1)$ and $z=ax/r$. The positivity condition is immediate from $a x\le2r$. Substitution gives
\[
 T=q+a x,\qquad
 D_1^{\max}-T
 =q/x-a^2/(1-x)-q-a x
 =a[r(2/x-1/2)-x]>0.
\]
In the last step $r\ge1$ and $1/2\le x<1$, so the bracket is at least $1/2$. The chosen $y$ is admissible and the controlled carrier obeys the exact recursion
\[
 q_{j+1}=q_j+a-\frac{2a^2}{\sqrt{q_j}},\qquad q_{j+1}-q_j\ge a/2.
\]
After $k$ updates, $q_j\ge q_{\rm start}+(a/2)j$, so the cumulative deficit from $ak$ is bounded by
\[
 2a^2\sum_{j=0}^{k-1}(q_{\rm start}+aj/2)^{-1/2}
 \le \frac{2a^2}{\sqrt{q_{\rm start}}}
       +8a\big(\sqrt{q_{\rm start}+ak/2}-\sqrt{q_{\rm start}}\big).
\]
The fixed startup contributes only an $s$-dependent constant to the deficit, and $F_N=4q_N$. This completes the constructive lower bound for fixed $s$ and arbitrary large finite $N$. The $\sqrt N$ deficit here is an upper estimate for this policy's loss, not an optimal convergence-rate lower bound.

## Unrestricted upper, checked without accessing the bath

Use the two equivalent Kraus representatives
\[
 K_0(\theta)=|0\rangle\langle0|+\sqrt{s}e^{i\theta}|1\rangle\langle1|,
 \qquad K_1(\theta)=\sqrt\ell\,e^{-ia\theta}|0\rangle\langle1|.
\]
The jump phase is a gauge change, not part of the actual controls. Their isometry $V$ into a purifying environment satisfies $V^\dagger\dot V=0$ and $\dot V^\dagger\dot V=aP$ at the nominal phase (in fact throughout). For an arbitrary purified state/tangent incident on a slot, the squared derivative norm after $V$ is the previous squared norm plus $\|\dot V\psi\|^2\le a$; the cross term vanishes. Phase-independent intervening control isometries preserve the previous norm. From parameter-independent preparation the final global derivative norm squared is at most $aN$. Pure-state QFI is no greater than four times that norm squared, and tracing away the fictitious environments cannot increase SLD QFI. So arbitrary finite ancillas and adaptive fixed-order control obey $F_N\le4aN$. Both bounds use the same discrete-call resource. This is the standard zero-$\beta$ channel-extension argument in this case, not a claim to a new upper bound.

## Closest published precedents and replay

- [Zhou–Jiang, PRX Quantum 2, 010343 (2021), §VII.2 Eqs. 52–58, §III and Fig. 2](https://arxiv.org/html/2003.10559) give the known optimal coefficient $4s/(1-s)$ and even a three-qubit (two ancilla) **local code per probe** producing a weak effective dephasing channel. Their asymptotic protocol then uses a growing spin-squeezed register of logical probes. The claim here is the same coefficient with three qubits **total** between uses, by a different, local tangent alignment and erasure mechanism.
- [Liu–Yang, Quantum 8, 1571 (2024), §3.2 Eq. 23, Table 2 and Fig. 4](https://arxiv.org/html/2403.09519v3) address precisely the same discrete channel and fixed-dimensional CPTP controls. They report nearly linear finite numerical behavior with **one** noiseless ancilla, and explicitly leave theoretical achievability for that one-ancilla case open. This two-ancilla theorem supplies a rigorous fixed-size construction and proves the known full coefficient; it does not solve their one-ancilla question or prove minimality. Their discussion of previously known unbounded routes is not an impossibility proof for bounded memory.
- [Zhou, PRL 133, 170801 (2024), Theorems 1 and 5 and Fig. 1](https://arxiv.org/html/2402.18765) gives an $O(1)$ no-go without noiseless ancilla for strictly contractive channels, a class including this damping. That does not prohibit two noiseless qubits with nonunital reset. [Sekatski et al.](https://arxiv.org/abs/1603.08944) treats dissimilar fast-control continuous-time resources, and [Li–Tong–Bao](https://arxiv.org/abs/2610.31249) studies dissipative reset including an erasure-assisted Ramsey setting. Neither should be confused with this fixed-survival discrete protocol. No exhaustive priority claim.

Run `python3 verify.py` with SymPy (exact test runtime about one second; full log and checksums included). The replay constructs 8-dimensional pre-use isometries with Cholesky factors (instead of matrix square roots), the two physical amplitude-damping Kraus matrices and an explicit trace-preserving erasure. It verifies the unflagged density and antisymmetric derivative after erasure, three seeds, an actual startup transition, a non-diagonal large-mode transition, and exact rational feasibility/startup threshold checks for representative survivals. These controls do not replace the continuum estimates.

**Recommendation:** Sound and materially distinct from existing asymptotic bounds and direct near-linear fixed-size numerics. Use a narrow lead such as “The known optimal phase-sensing slope under discrete amplitude damping is achievable with a three-qubit reusable register and no retained syndrome history.” Cite both closest predecessors and preserve the assumptions of local pointwise estimation and arbitrary reset-enabled controls.
