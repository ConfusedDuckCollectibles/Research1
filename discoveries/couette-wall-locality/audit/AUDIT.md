# Finite width versus one-wall scalar propagators: independent audit

**Audit passes, with one corrected constant.** For $\partial_s u=\partial_y^2u+i a(s)y\,u$ with real $a\in L^\infty(0,t)$ and diffusivity one, let $U_{L}$ be its $L^2$ propagator on $(0,L)$ with either homogeneous Dirichlet or Neumann independently at each wall. Let $U_{+,D}$ and $U_{+,N}$ be the corresponding half-line propagators, and let $B_0,B_L$ be the two finite wall types. Uniform absolute numerical constants give
\[
 \big|\|U_L\|-\max(\|U_{+,B_0}\|,\|U_{+,B_L}\|)\big|
 \le C t/L^2+C e^{-c L^2/t}.
\]
For a right wall, reflection changes $a$ to $-a$ and adds only a scalar phase; complex conjugation gives the same norm. No positivity or self-adjointness of the imaginary-potential evolution is used. This is an **additive**, not relative, approximation. It applies to arbitrary bounded time histories, with a bound independent of their amplitude or derivatives. The case $t=0$ is exactly identity.

Each half-line essential norm is the known free-space Kelvin norm $M$, and in fact the essential spectrum of $U_+^\dagger U_+$ is $[0,M^2]$. Every half-line singular pair of singular value $\sigma>M$ is exponentially localized, with any rate
\[
 \gamma<\sqrt{\log(\sigma/M)/t}
\]
for the $L^2$ weight $e^{\gamma y}$. If $\|U_+\|>M$ it is attained, but strict wall enhancement is not asserted for each history. The originally suggested bound $2p_t(x+y)$ on the Neumann-minus-free kernel is **false** uniformly over bounded histories; the correct uniform sharp bound is $3p_t(x+y)$, where $p_t(r)=(4\pi t)^{-1/2}e^{-r^2/(4t)}$.

## Norm comparison and mixed walls

Use Feynman–Kac with the phase of modulus one. The reference process is Brownian motion of variance $2s$ with normally reflecting Neumann walls and killed at any Dirichlet wall. On the interval this is periodic 1-Lipschitz folding of Brownian motion, optionally killed at the first visit to specified walls. On the half-line use absolute value for the reflecting wall. If the evolution kernel is written with the Brownian path started at the output coordinate, use the reversed gradient history $a(t-s)$ in the phase. This convention does not affect the estimates. Existence for bounded measurable $a$ also follows by step-history approximation and dominated convergence, since each sample path on a compact time interval is bounded. Killing or folding and the unit-modulus phase yield the Gaussian second-moment estimate
\[
 \sup_x\int |x-y|^2 |K(x,y)|\,dy\le 2t,\qquad
 \sup_y\int |x-y|^2 |K(x,y)|\,dx\le 2t .
\]
The first follows from displacement of the folded process no larger than that of its Brownian driver; the second uses symmetry of the underlying diffusion path measure and the time-reversed history. The same domination gives $L^2$ contraction by Schur.

### Quadratic cutoff identity

Let $Vf=(\chi_j f)_j$ for a real quadratic partition $\sum_j\chi_j^2=1$, and let $\operatorname{Lip}(x\mapsto(\chi_j(x))_j)\le A$. At the **propagator** level itself,
\[
 U-V^\dagger\operatorname{diag}(U) V
 \quad\hbox{has kernel}\quad
 \frac12\sum_j(\chi_j(x)-\chi_j(y))^2 K(x,y).
\]
This uses only the pointwise quadratic partition identity, with no derivative of the gradient history, no generator commutator integration, and no self-adjointness. Both Schur bounds are $\le A^2 t$. More generally the comparison blocks may be different operators, each compressed to the appropriate cutoff support: the norm of their recombination is at most the supremum block norm as $V$ is isometric.

If two geometries agree until hitting a boundary at distance at least $h$ from a cutoff support, couple the two paths with the same Brownian driver before that first unmatched hit. The integrals and positions agree on good paths. For one fixed endpoint in that support, the row total variation of the difference kernel restricted to the other endpoint in the support is at most $2P(\text{unmatched hit by }t)$ (and analogously for columns by reversed histories). As folding is 1-Lipschitz and killing cannot make a hit more likely, this is $\le8e^{-h^2/(4t)}$ using the Brownian maximum inequality. Schur gives the corresponding operator bound on that support. This remains valid for a matched Dirichlet wall that kills both paths before the unmatched wall.

* **Finite upper:** take a sine/cosine pair constant up to $L/4$ and after $3L/4$, with left and right supports separated from the other wall by $L/4$. Compare each block to the half-line retaining its matching wall. The vector cutoff Lipschitz constant is $O(1/L)$.
* **Finite lower:** localize each half-line with mesh $h=L/4$. On each $[jh,(j+1)h]$ use exactly $\chi_j=\cos(\pi(x-jh)/(2h))$ and $\chi_{j+1}=\sin(\pi(x-jh)/(2h))$; then extend by zero on the other cells. Individual cutoffs have support $[(j-1)h,(j+1)h]$ intersected with the half-line, and the vector Lipschitz constant is $\pi/(2h)$. The first two support blocks are within $[0,L/2]$ and compare directly to the matching finite wall. All subsequent support blocks are at distance at least $h$ from the half-line wall, so compare to the full line, translate the support of length $2h$ into $[h,3h]$, then compare to the finite interval in its interior. Translation of a full-line linear phase adds only $e^{i c\int a}$ for translation distance $c$, preserving the block norm. Both comparisons incur only Gaussian rare-hit error. Even infinitely many half-line blocks recombine with the **supremum** error via the isometry above, not a sum of per-block errors. This yields $\|U_{+,B}\|\le\|U_L\|+O(t/L^2+e^{-cL^2/t})$ for either matched wall, without assuming attainment or a wall-localized maximizer.

The quadratic rate cannot uniformly be replaced by $o(t/L^2)$ over all allowed walls: at $a=0$ and D–D, $\|U_L\|=e^{-\pi^2 t/L^2}$ while $\|U_{+,D}\|=1$.

## Essential spectrum and weighted singular vectors

Set $A(s)=\int_s^t a(r)\,dr$ and $\overline A=t^{-1}\int_0^t A(s)\,ds$. Fourier characteristics give
\[
 \widehat{U_{\mathbb R}f}(k)=
 M e^{-t(k-\overline A)^2}\,\widehat f(k-A(0)),\qquad
 M=\exp\{-\int_0^t (A(s)-\overline A)^2\,ds\}.
\]
This verifies the exact magnitude bound $|K_{\mathbb R}(x,y)|=M p_t(x-y)$, not just a diffusion domination. Write the half-line bulk compression $T_0=P_+U_{\mathbb R}E_+$, where $E_+$ extends by zero and $P_+$ restricts. A Brownian path joining positive $x,y$ that hits zero has unweighted endpoint density $p_t(x+y)$ by reflection. Dirichlet omits this contribution, so its kernel difference from $T_0$ is $\le p_t(x+y)$. Neumann via the absolute-value path has the additional endpoint $-y$, at density $p_t(x+y)$, and the phase histories for positive-endpoint paths that hit zero may differ by a unit-phase difference of modulus at most 2. Hence the bound is **three**, and the correction $D=U_+-T_0$ is Hilbert–Schmidt, with $\|D\|_{\mathrm{HS}}\le3/\sqrt{4\pi}$ (Dirichlet uses 1), since $\int_0^\infty z p_t(z)^2\,dz=1/(4\pi)$. Left weighting by $e^{\gamma x}$ preserves the Hilbert–Schmidt property for each fixed $\gamma$; the same holds for $D^\dagger$.

Because the full kernel phase separates in $x,y$, $T_0=M A_1 H_+ A_2$ for unitary phase multipliers $A_1,A_2$ and $H_+=P_+e^{t\partial^2}E_+$, a positive self-adjoint contraction. For any $e^{-t\xi^2}\in(0,1]$, normalized packets $n^{-1/2}\phi(x/n)e^{i\xi x}$ with $\phi$ supported in $[1,2]$ are weakly null approximate eigenvectors of $H_+$ by full-line Fourier concentration. Hence the essential spectrum of $H_+$ is $[0,1]$, and the compact boundary difference gives $\sigma_\mathrm{ess}(U_+^\dagger U_+)=[0,M^2]$. One may alternatively obtain only the essential norm by translating any compactly supported near-maximizing input of $U_{\mathbb R}$ arbitrarily far from the wall.

For the weight, use the self-adjoint singular block $\mathcal H=\left(\begin{smallmatrix}0&U_+^\dagger\\U_+&0\end{smallmatrix}\right)$ and its bulk part $\mathcal H_0$. A singular pair $v$ satisfies $(\sigma-\mathcal H_0)v=(\mathcal H-\mathcal H_0)v$. The weighted-conjugate bulk kernels have magnitude $M e^{\gamma(x-y)}p_t(x-y)$, whose full-line convolution bound is $M e^{\gamma^2t}$ by the Gaussian moment. Thus solve the weighted block equation in $L^2$ by a Neumann-series resolvent whenever $M e^{\gamma^2t}<\sigma$; its right side is in $L^2$ by the weighted HS estimate. Multiplication by bounded $e^{-\gamma x}$ intertwines the bounded kernels and returns the unweighted equation; uniqueness proves the weighted solution is exactly $e^{\gamma x}v$. This avoids applying an unbounded exponential weight to an unproved localized vector. All $\sigma>M$ are isolated eigenvalues of finite multiplicity in the singular block; in particular strict above-bulk operator norm is attained. The result is finite-time singular-vector localization, not localization of Bloch–Torrey generator eigenvectors.

## Sharp Neumann factor, with bounded histories

Take horizon one and the free bridge $\beta$ from 0 to 0 with diffusion variance $2s$. At the continuous boundary kernel value,
\[
 \frac{K_N(0,0)-K_{\mathbb R}(0,0)}{p_1(0)}
 =2\mathbb E e^{i\int a|\beta|}-\mathbb E e^{i\int a\beta}.
\]
Write $a$ in Brownian time, reversing the chosen experimental history if necessary. Let $s_j=4^{-j-1}$, $\sigma_j=\sqrt{2s_j(1-s_j)}$, $\xi_j=\beta(s_j)/\sigma_j$, $\mu=\sqrt{2/\pi}$. The normalized bridge coordinates are centered unit Gaussians with correlations $0\le\rho_{jk}\le2^{-|j-k|}$. Define ideal pulse phases
\[
 G=\frac{\pi}{m\mu}\sum_{j=1}^m\xi_j,\qquad
 H=\frac{\pi}{m\mu}\sum_{j=1}^m|\xi_j|.
\]
Their means are 0 and $\pi$. Their scaled averages have variances at most $3/m,4/m$ respectively: for absolute values, couple $Y=\rho X+\sqrt{1-\rho^2}Z$ to independent $Z$, giving $|\operatorname{Cov}(|X|,|Y|)|\le\sqrt2|\rho|$ by $L^2$ Lipschitz and Cauchy–Schwarz. Sum the geometric correlations. The standard Taylor inequality for $e^{ix}$, taking expectation of the centered variable, gives
\[
 \big|2\mathbb E e^{iH}-\mathbb E e^{iG}+3\big|
 \le\operatorname{Var}H+\tfrac12\operatorname{Var}G
 \le\frac{11\pi^3}{4m}.
\]
Replace each point pulse by a nonoverlapping rectangular pulse on $[s_j,s_j+\eta^2s_j(1-s_j)]$ of the same mass $\pi/(m\mu\sigma_j)$. These are bounded for every fixed finite $m$ and $\eta>0$. Since $\mathbb E|\beta(u)-\beta(s_j)|^2\le2|u-s_j|$, this changes each full phase expectation by at most $\pi\eta/\mu$, also for absolute value. The combination changes by at most $3\pi\eta/\mu$.

The formula extends to strictly positive endpoints using bridge means affine in endpoints and $(x,y)\mapsto(x,-y)$ for the second reflected branch. At $x=y=\epsilon$ changing any bridge expectation from the boundary value costs at most $\epsilon\|a\|_1$. The ratio $p_1(0)/p_1(2\epsilon)=e^{\epsilon^2}$ costs at most $4\epsilon^2$ in the normalized kernel difference for $\epsilon^2<1/2$. Taking $m=256$, $\eta=1/100$, $\epsilon=2^{-276}$ gives total error below $1/2$ from $-3$, using $\pi<22/7$, $\mu>3/4$, and $\sigma_j\ge\sqrt{s_j}$ to bound the finite masses. See exact rational replay. Choosing $m\to\infty$, $\eta\to0$ and suitably tiny positive endpoints proves the supremum factor exactly 3, and continuity shows failure of the proposed factor 2 for a positive-measure neighborhood with a fixed bounded history. This sharpness is uniform across all bounded histories; it makes no claim about a class with a prescribed common amplitude budget.

## Priority and scope

- [Stoller–Happer–Dyson, Phys. Rev. A 44, 7459 (1991)](https://doi.org/10.1103/PhysRevA.44.7459) already solves the constant-gradient reflective half-line/interval Torrey problems by complex Airy analysis. [Grebenkov, localization regime (2014), Eqs. 5 and 13](https://pmc.polytechnique.fr/pagesperso/dg/publi/2014_10.pdf) records the exact free Gaussian constant-gradient propagator, explicitly describes the far-wall recovery of bulk behavior, and works with barriers/echo signals. Do not present the constant-gradient equation, free $e^{-a^2t^3/12}$ law, or existence of physical boundary localization as new.
- [Herberthson–Özarslan–Knutsson–Westin, *Dynamics of local magnetization in the eigenbasis of the Bloch–Torrey operator*, §§II–III](https://doi.org/10.1063/1.4978621) treat general gradient waveforms including slab/echo examples, chiefly to compute the magnetization and integrated signal. That signal is not the largest singular value over all $L^2$ initial fields. [Özarslan et al., *Characterizing magnetic resonance signal decay due to Gaussian diffusion*, §3](https://pmc.ncbi.nlm.nih.gov/articles/PMC4864615/) treat the arbitrary-waveform free result by path integrals.
- [van den Berg, *Heat equation and the principle of not feeling the boundary* (1989)](https://people.maths.bris.ac.uk/~mamvdb/pdfs/Publication27.pdf) is a close tool-level precedent for Gaussian remote-wall error. Quadratic IMS localization, Schur estimates, and compact-perturbation arguments are standard. [Grebenkov–Helffer–Henry](https://arxiv.org/abs/1603.06992), [Helffer–Sjöstrand](https://arxiv.org/abs/1001.4171), and [Averseng–Frantz–Hérau–Raymond](https://arxiv.org/abs/2512.20202) address related constant-gradient resolvents, semigroup decay, and/or eigenfunction localization, not the specific waveform-uniform finite-time width comparison or the singular essential edge established here. The targeted search did not locate a direct equivalent; no exhaustive priority claim.
- The claim is a scalar PDE/operator statement for D/N/no-flux ideal endpoints, not a model of a complex laboratory or medical acquisition. It may be presented narrowly as a uniform operator-norm wall comparison, an exact essential singular spectrum/localization bound, and the sharp 3 correction. It is not a quantitative claim about tissue, reach, or practical ultra-concentrated pulses.

Replay: `python3 sharp_constant.py | diff - replay.txt` (standard-library rational inequalities only). It certifies the finite bounded-history counterexample through the analytic Gaussian bounds above; it is not a Monte Carlo test.
