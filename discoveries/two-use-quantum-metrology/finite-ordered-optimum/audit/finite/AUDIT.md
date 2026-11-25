# Exact local tangent recursion for any finite fixed order

**Audit passes.** For any prescribed finite sequence of independently amplitude-damped qubit phase calls, pointwise SLD QFI optimizes by an exact Bellman recursion in one nonnegative tangent norm. A two-by-two positive contraction is the complete instantaneous control variable; real symmetric contractions with nonnegative off-diagonal suffice. The maxima are attained. An implementation uses one routed sensor, two noiseless qubits and at most one classical bit per call, without reading a damping environment. This extends the two-use result by a new reduction, not by a claim that exact finite adaptive optimization did not previously exist.

## Local problem and exact transition

Fix a nominal phase, taken to be zero. For an ordered suffix $\mathbf v=(v_1,\ldots,v_n)$, let $J_{\mathbf v}(d)$ be the supremal **total** pointwise SLD QFI when the incoming normalized state is pure and its horizontal ket derivative has norm $d$, before arbitrary phase-independent processing. There is just one unknown real phase, shared by all calls. The $v_i\in[0,1]$ are known survival probabilities; each scheduled call is used exactly once, independently, with no postselection. Only phase-independent free processing, ancillary memory and a local design at the nominal phase are admitted. The actual schedule objective is $J_{\mathbf v}(0)$. For nonzero $d$, arbitrary normalized vectors with horizontal derivative are unitarily equivalent as local state-and-derivative pairs, so the value depends only on $d$.

Write the state and a horizontal representative of its derivative as $\psi$ and $i d\phi$, with $\psi,\phi$ orthonormal. At $d=0$ choose an arbitrary auxiliary orthogonal direction. For a call of survival $v$, let $\ell=1-v$ and $Q=|1\rangle\langle1|\otimes I$ on the routed sensor and memory. Compression in basis $(\psi,\phi)$ gives
\[
 M=\begin{pmatrix}x&z+iw\\z-iw&y\end{pmatrix},\qquad 0\preceq M\preceq I.
\]
Temporarily expose the two Kraus outcomes. For each unnormalized output ket $a_j$, put $q_j=\|a_j\|^2$, $b_j=\langle a_j,\dot a_j\rangle$, $N_j=\|\dot a_j\|^2$. The Kraus operators are $B_0=\operatorname{diag}(1,\sqrt v)$ and $B_1=\sqrt\ell\,|0\rangle\langle1|$, each followed by phase rotation on the sensor; the jump-output excitation is zero. Thus $\dot a_j=i Q a_j+i d B_j\phi$. It follows directly, with $m=z+iw$, that
\[
\begin{array}{lll}
q_0=1-\ell x,&b_0=i(vx-\ell d m),&N_0=vx+d^2(1-\ell y)+2vdz,\\
q_1=\ell x,&b_1=i\ell d m,&N_1=\ell d^2 y.
\end{array}
\]
Hence $\dot q_j=2\Re b_j$, in particular $(\dot q_0,\dot q_1)=(2\ell d w,-2\ell d w)$. At $q_j>0$, dividing by $\sqrt{q_j}$ and subtracting the remaining global-phase derivative gives
\[
 d_j^2=N_j/q_j-|b_j|^2/q_j^2.
\]
The formula is independent of the global-phase choice. The recurrence, with the existence of a maximum justified below, is
\[
 J_{\varnothing}(d)=4d^2,\qquad
 J_{(v,\mathbf u)}(d)=\max_{0\preceq M\preceq I}
 \sum_{j:q_j>0}\left(\frac{\dot q_j^2}{q_j}+q_j J_{\mathbf u}(d_j)\right).
\]

### Upper, chain rule, and physical attainment at fixed controls

Purify preparations and each phase-independent intermediate CPTP control, retaining its auxiliary system. Refining each noisy call by a classical Kraus label only strengthens the model because one can ignore all labels and recover the original circuit. Every positive-weight conditional state is now pure, so its nominal state and one horizontal derivative fit in their two-dimensional span regardless of all previous memory. The resulting optimization cannot exceed the displayed transition with arbitrary $M$, followed by conditionally optimal suffixes.

For a positive block $\rho_j=q_j\sigma_j$ of a flagged state, a conditional SLD $L_j$ for $\sigma_j$ gives block SLD $L_j+(\dot q_j/q_j)I$, and $\operatorname{Tr}(\sigma_j L_j)=0$. This proves the exact classical + conditional split and therefore the classical-history chain rule when applied recursively. If $q_j=0$, then $a_j=0$ and $\dot\rho_j=|\dot a_j\rangle\langle a_j|+|a_j\rangle\langle\dot a_j|=0$, even if $\dot a_j\ne0$. A subsequent differentiable parameter-dependent linear map acts by the product rule on a zero density and zero derivative and cannot revive this branch at first order. There is no conditional normalization or Fisher score at zero weight under pointwise SLD. A continuous Bures or Hellinger convention is not the theorem.

For the reverse inequality, **every** effect $M$ can be physically realized on the two-dimensional span. Use the isometry into sensor, one tag, and one logical qubit
\[
 W_M=|0,0\rangle\otimes\sqrt{I-M}+|1,1\rangle\otimes\sqrt M.
\]
It has $W_M^\dagger W_M=I$ and compresses sensor excitation to $M$. The no-jump operator maps the even sensor/tag code into even parity and the jump operator into odd parity; the same holds for their derivatives under this phase. Measuring only accessible parity therefore preserves the full flagged nominal density and density derivative. The actual family is not required to remain in the encoded span away from the nominal phase: the input ket and its first derivative are in that span at the nominal phase, which suffices by differentiating any subsequent fixed controls and smooth slots. Each chosen two-vector isometry extends to a fixed unitary on this same eight-dimensional sensor/tag/logical register (complete both orthonormal pairs to bases), separately for each prior recorded outcome. The qubits can be reused; no register that scales with the protocol depth is required except classical flag history (at most $n$ bits). Keep all flag branches and make the final local measurement conditionally. There is no unphysical access to a bath trajectory, no $\theta$-dependent feed-forward, and no postselection. At this step one can first reason about suprema, then prove actual existence.

## Real controls and attainment

Three general facts, for an $n$-call suffix:
1. **Tangent splitting:** for $d\ge d_0$, re-embed the pure local pair into a product of a logical tangent of norm $d_0$ and a parked, noiseless ancillary tangent of norm $\sqrt{d^2-d_0^2}$. The orthogonal product derivative has total norm $d$ and defines an isometry on the given pair. Apply any suffix policy to the first factor and preserve the second. At the nominal point the output and derivative equal those of the corresponding product family (which suffices for SLD), and its QFI is additive. Thus
\[
 J_{\mathbf u}(d)\ge J_{\mathbf u}(d_0)+4(d^2-d_0^2).
\]
In particular $J_{\mathbf u}(d)-4d^2$ is nondecreasing. This is a control argument, not preparation of a new unknown-phase resource: the parked derivative is transferred from the incoming tangent.

2. **Continuity:** on the canonical local pair the nominal density is independent of $d$. For each fixed strategy the derivative is affine in $d$, and the difference at $d,d'$ is the input tangent difference through the nominal strategy's CPTP channel. The square root of SLD QFI is a seminorm of the derivative at fixed state (spectral expression with the zero denominator terms omitted), and contractivity bounds this difference by $2|d-d'|$. Hence the supremum of these square roots is $2$-Lipschitz.

3. **Growth:** purify a fixed strategy and each channel, retaining the environments just for this bound. Use the equivalent centered phase generator $P-I/2$ with operator norm $1/2$. The derivative of the purified output has norm at most $d+n/2$ by the product rule, the triangle inequality and isometries. Thus $4d^2\le J_{\mathbf u}(d)\le4(d+n/2)^2$, with the lower bound obtained by preserving the input and sending dummy probes.

Replacing $z+iw$ by $z$ preserves both PSD conditions (it reduces the off-diagonal magnitude) and preserves $q_j,N_j,\Im b_j$. Removing $\Re b_j$ deletes $\dot q_j^2/q_j$, while increasing $4q_jd_j^2$ by exactly that amount. The other part of $q_jJ_{\mathbf u}(d_j)$ cannot decrease by Fact 1. In a zero-probability branch the omitted off-diagonal is already zero unless the outcome is identically absent by $v=1$, in which case the active branch has no classical score. Replacing $z$ by $|z|$ then preserves feasibility, keeps the positive jump tangent norm, and cannot lower the no-jump tangent norm: in $q_0^2 d_0^2$ its sign-dependent term is precisely $2vdz$. This conclusion is also valid at $v=0$ and $v=1$ when a vanishing block is omitted. Thus the maximum may be restricted to three real coordinates $(x,y,z\ge0)$ with $z^2\le xy$ and $z^2\le(1-x)(1-y)$.

For attainment with $0<v<1$, $q_0\ge v>0$, and the only excluded mass is $q_1=\ell x=0$. Let effects with $q_1>0$ converge to $x=0$ and limiting $y$. Positivity forces $m\to0$, so the no-jump contribution tends to $J_{\mathbf u}(d_0)$ for $d_0=d\sqrt{1-\ell y}$. From Fact 3 and the branch formulas, the jump contribution is bounded by
\[
 \frac{\dot q_1^2}{q_1}+q_1 J_{\mathbf u}(d_1)
 \le 4N_1+4n\sqrt{q_1N_1}+n^2q_1
\]
because $\dot q_1^2/q_1+4q_1d_1^2=4(N_1-(\Im b_1)^2/q_1)$ and $q_1 d_1\le\sqrt{q_1 N_1}$. Its possible limiting contribution is at most $4\ell d^2 y$. Fact 1 therefore bounds such a singular limiting objective by $J_{\mathbf u}(d)$, already feasible at $M=0$ by parking the incoming tangent. Sequences contained in the omitted-mass face similarly have objective at most $J_{\mathbf u}(d)$ by Fact 1. Compactness of the effect set and continuity elsewhere now give an attained maximum. At $v=1$ the jump is identically absent and the objective is continuous. At $v=0$ the unflagged slot is a phase-independent replacer and equivalent to free processing plus a mandatory dummy call, so $J_{(0,\mathbf u)}(d)=J_{\mathbf u}(d)$, attained with $M=0$. The already-proved physical realization of every effect prevents a larger value from the ideal-flag formula even at this endpoint. Backward induction supplies an attaining conditional protocol for every finite schedule.

## Independent controls and literature

Run `python3 verify.py` (Python + SymPy); see `replay.txt`. The independent construction uses explicit 8-dimensional party/tag/logical density matrices and directly solves the rank-one-block SLD equation. It checks a complex effect with both a classical score and a negative real off-diagonal, two distinct conditional second-call controls with explicit amplitude differentiation, a discontinuous zero-mass control, and both endpoint survivals. This is a consistency replay, not the proof of all dimensions or an algorithmic complexity bound. A continuous scalar value function with an exact local max need not be computationally cheap as $n$ grows; no efficient approximation or memory lower bound is asserted.

- [Liu et al., PRL and supplements E.2/G.1](https://arxiv.org/html/2203.09758) and [their 2024 review, §§3.2, 4 and 5.1](https://arxiv.org/html/2409.07068) already give exact finite strategy optimization by quantum-comb SDP, including $N=3$ amplitude-damping numerics and nonidentical-channel examples. This note adds a specialized scalar-state, two-by-two-control recursion and constant-quantum-register construction, not the first exact fixed-order formalism.
- [Tang et al. § II](https://arxiv.org/html/2602.09358) already compress pure unitary phase-encoded QFI into one qubit and classical outcomes by an energy-basis measurement, with sequential compression of equatorial qubit copies. That compresses an already encoded phase family, not this call-dependent noisy control optimization. Here the same nominal state and density derivative are preserved through each local step via control on the first-order span; arbitrary unknown-phase family transport is not claimed.
- [Albarelli et al.](https://quantum-journal.org/papers/q-2018-12-03-110/) and [Fallani et al., App. A](https://doi.org/10.1103/PRXQuantum.3.020310) already use classical trajectory FI plus average conditional QFI for monitored systems. That identity and square-root effect dilation are standard. The special contribution is their exact use at each step of this unmonitored-channel optimization via an accessible tagged parity implementation, including arbitrary compressed effects and strict rank-change/attainment control. The targeted prior search is not exhaustive.
