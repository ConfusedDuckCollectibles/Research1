# Independent audit: all-path CHSH corner profile

**Recommendation.** The original all-path result is correct under its stated model. The alternative proof in `uniform-proof.md` strengthens the remainder to an explicit $t/5$ normalized error simultaneously over all approach ratios. This supports a narrowly framed mathematical note about the margin in an already solved optimization problem. It is not a new loss threshold, exact CHSH optimization result, device-independent self-test, or certification claim.

## Model and reduction

With independent local erasures of fixed, state/input/outcome-independent efficiencies and every no click assigned $+1$, the one-sided erasure terms are $2\eta_A(1-\eta_B)\langle A_0\rangle$ and $2(1-\eta_A)\eta_B\langle B_0\rangle$ in the canonical CHSH sum; erasing both contributes $2(1-\eta_A)(1-\eta_B)$. The nonerased contribution is $\eta_A\eta_B\,\mathrm{CHSH}$. This gives the stated tilted and recorded scores, including their normalization. It is a model of a quantum prediction, not an assumption about the local alternatives under a properly designed, full-data Bell test. No outcomes are discarded. The usual Bell assumptions are still necessary; this does not address their implementation or finite-trial significance. It also says nothing about correlated erasures, using a third outcome, optimizing other tests, or quantitative entanglement certification.

For two binary POVMs, use the self-adjoint unitary dilation $R(T)=\left(\begin{smallmatrix}T&\sqrt{I-T^2}\\ \sqrt{I-T^2}&-T\end{smallmatrix}\right)$ for each contraction $T$, with the same embedding of the original space in the first summand for both settings. All terms here use at most one local operator per party, so the dilations preserve the scores. Then Jordan decomposition and averaging isolate reflection blocks of size at most two on each party. Finite-dimensional reduction also bounds the supremum over unrestricted ordinary bipartite dimensions by finite-rank approximations. If one party's block is one-dimensional, that party is deterministic and the other party's means enter linearly, hence the maximum is local. Within nontrivial blocks, rotate $A_0=B_0=Z$ and choose each off-diagonal axis real. The parameters $1-u$ and $1-v$ are the cosines and the accessible box is **all** $u,v\in[0,2]$, with sine magnitudes $\sqrt{u(2-u)}$, $\sqrt{v(2-v)}$. This includes obtuse angles and commuting endpoints. Local bound $L=2+a+b$ follows from CHSH at most 2 and the separate positive-tilt terms, all saturated by all-positive outputs.

## Independent determinant and far-corner control

The code independently constructs the four CHSH terms and both tilts by tensor products. To do exact rational manipulations, the physical qubit matrix with off-diagonal $\sin\theta$ is diagonally conjugated by $\operatorname{diag}(1,\sin\theta)$: its off-diagonals become $\sin^2\theta$ and 1. The determinant, and all principal minors of this diagonally conjugated tensor matrix, equal the physical ones; the vanishing-sine endpoints follow by polynomial continuity. Expanding via signed permutations gives exactly the candidate determinant. This avoids importing the candidate's formula as matrix input or computing any floating eigensystem.

For $0<a<t\le1/100$ the lower three-dimensional compression of $M(0)$ has diagonal part bounded below by $(4-2s)I$. The Kronecker perturbation has norm $2\sqrt{uv}$, making that compression positive definite when $uv<(2-s)^2$. Otherwise $uv\ge(2-s)^2\ge(1.99)^2$ gives $e=2-u,f=2-v\in[0,.02]$. Its internal block on the $01,11$ coordinates has diagonal bounds $3.98,7.94$ and off-diagonal magnitude at most $.4$. For a $2\times2$ symmetric matrix with positive determinant, each diagonal inverse entry increases with the magnitude of the off-diagonal and decreases with each diagonal, so the substituted rational determinant bounds give the stated upper bounds $.30,.13,.02$ for inverse entries (off-diagonal in magnitude). The coupling from coordinate $10$ has squares at most $.08f$ and $8f$. The upper inverse quadratic form is therefore at most
$$
(.30\cdot.08+.13\cdot8+2(.02)\sqrt{.08\cdot8})f=1.096f<1.1f\quad(f>0);
$$
it vanishes if $f=0$. The Schur complement is at least $2a+2e+.9f-ef>0$, since $a>0$ and $e\le .02$. This covers the opposite-angle corner where a small-angle series alone would be insufficient. Adding $dI$ preserves positivity of the lower compression, so the full determinant sign diagnoses the unique potentially nonpositive eigenvalue.

## Asymptotic proof audit

The original first localization uses only the full square: $A-2v\ge-a^2$ makes $v\le 2h[1+O(t)]+O(t)$ on positive-gap blocks. On any fixed $h\le 1-\epsilon$ this makes $A-2v$ uniformly positive, yielding the indispensable sharper $v=O(h)$ and $u=O_\epsilon(s)$. Dividing $F(0)$ by $s^2h$ then has controlled error even if $h$ goes to zero faster than any power of $t$. The derivative bound at the root, ensuing profile, and matching admissible angle choice all check. The original convex squeeze at fixed $s$ is valid near $h=1$: the supremum of affine functions in $a$ is convex, party exchange gives midpoint symmetry, and on the left half it is nonincreasing. Therefore the original endpoint uniformity was proved, without presuming uniform constants in its interior expansion. The alternative exact envelope additionally proves a uniform rate and avoids switching arguments between two scales.

## Exact replay and relation to the known maximum

`exact_audit.py` independently certifies all 17 supplied rational gap intervals. Lower endpoints use the direct rational determinant at the supplied physical angles. For upper endpoints put $J=sC-ab\,d$. If positive gain is possible then $v<(sC+2a^2)/B<A/2$ and $J-Bv>0$. Completing the square in $u$ reduces strict upper certification to positivity of the cubic
$$
P(d)(A-2v)-2v(J-Bv)^2 \quad\text{on }[0,J/B].
$$
The program checks positive endpoints and absence of any roots by a rational Sturm chain, independent of the supplied square-root envelope calculation. It also checks six new $t/5$ bracketing controls, including $t=1/100$ with ratios within $10^{-12}$ of each endpoint and $t=10^{-40}$, and twelve near-opposite-angle compression controls by exact principal minors. These are audits of the analytic proof, not its substitute.

The determinant is **symbolically identical** to [Gigena et al., appendix C, Eq. (46)](https://arxiv.org/html/2405.08743v2#A3) after writing $c_A=1-u,c_B=1-v,\lambda=L+d$. Further, at an interior stationary point of the gain term the two derivative equations give $Au=Bv$. With $p=AB,z=v/A$ they reduce, together with zero determinant, to
$$
4pz^2-3pz+J=0,\qquad P(d)-8p^2z^3(1-2z)=0.
$$
An exact symbolic resultant, substituting $J=sC-ab\,d$ and $\lambda=L+d$, is $16p^4(ab+2\lambda)^2$ times the sextic of [Gigena et al., Eqs. (42)–(43)](https://arxiv.org/html/2405.08743v2#A3) with $\tau_5=-4ab$ as corrected by [Pakhunov §5](https://arxiv.org/html/2607.13774v1#S5). All 17 numerical enclosures isolate the largest real root of that corrected polynomial by an independent Sturm check, and none contains a root with the original printed sign. Thus the determinant/equations compute the same known optimum; the all-angle proof isn't publicity for a new polynomial.

## Priority and scope

- [Gigena et al., 2025 journal](https://doi.org/10.1038/s41534-025-01029-6) and their latest [arXiv v2](https://arxiv.org/abs/2405.08743) already give the loss-to-tilt map, critical region, optimal strategies and the exact maximum. The versions checked do not state a uniform double-scaling profile at $(a,b)=(0,2)$.
- [Pakhunov §2](https://arxiv.org/html/2607.13774v1#S2) states the cubic coefficient $(4-\beta)(2+\beta)/[54\beta(2-\beta)]$ for **fixed** interior $\beta$, and explicitly notes its divergence and crossover to the known quadratic law. The coefficient is $(8+ab)/(54ab)$ on $a+b=2$. Its $a\to0$ behavior agrees with $2/(27a)$; it does not by itself justify exchange of limits or give $g(h)$ at intermediate paths. Pakhunov's §6 scaling theorem is on the symmetric line $\alpha=\beta=1-s/2$, not this endpoint. The [companion paper](https://arxiv.org/html/2607.13762v1) deals with NPA overshoot along that symmetric line. No matching double-scale formula was located there.
- [Bamps–Pironio](https://arxiv.org/abs/1504.06960) already give the singly tilted quantum bound, hence the $1/8$ endpoint. [Wu, Xie and Ramanathan, appendix D.1](https://arxiv.org/html/2611.27413v1#A4.SS1), dated November 2026, compute fixed-score CHSH under a **common** $\eta$, plus compare entropy witnesses. Their explicit loss formula is the symmetric old tilt, not the asymmetric $(\eta_A,\eta_B)\to(1/2,1)$ approach ratio here. Neither Eberhard's nor [Brunner et al.'s](https://arxiv.org/abs/quant-ph/0702130) detection thresholds are changed.

Targeted primary checks found no statement of the explicit all-path $g(h)$ or the quantified error. Priority is not established. Its value is the compact, globally justified asymptotic readout of already known exact quantum scores.

To reproduce from this folder, run `python3 exact_audit.py` with Python 3 and SymPy. The proof of the new remainder does not depend on this program.
