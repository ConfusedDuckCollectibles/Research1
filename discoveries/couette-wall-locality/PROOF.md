# Single-wall limits of finite-time mixing in modulated Couette channels

## Model, comparison and theorem

Let \(S_I^a(t)\) be the propagator on \(L^2(I)\) from 0 to \(t>0\) for
\[
 \partial_s f(s,x)=\partial_x^2 f(s,x)-i a(s)x f(s,x),
 \qquad a:[0,t]\to\mathbb R
\]
bounded and measurable, with homogeneous Dirichlet or Neumann conditions, set independently at each finite endpoint. Write \(I_L=(0,L)\), \(N_L=\|S_{I_L}^a(t)\|_{2\to2}\), and \(N_D,N_N\) for the corresponding one-wall norms on \((0,\infty)\). Put \(m=\max(N_{\rm left},N_{\rm right})\), choosing each one-wall type as on the corresponding finite channel. Reflection conjugates the right-wall propagator up to a scalar phase, so its norm agrees with this convention.

**Single-wall norm bridge.** For every \(t,L>0\), every such history \(a\), and \(0<d<L/2\),
\[
 |N_L-m|\ \leq\ \frac{\pi^2 t}{(L-2d)^2}+16e^{-d^2/(4t)}. \tag{1}
\]
This holds independently of the amplitude, switching times and sign reversals of the shear. The leading constant is sharp uniformly over this class: the supremum discrepancy over all histories and D/N wall types is \((\pi^2+o(1))\,t/L^2\) as \(t/L^2\to0\). For the upper, use \(d=\sqrt{12t\log(L/\sqrt t)}\) when \(L/\sqrt t\) is large; for the lower, at \(a=0\) with two Dirichlet walls, \(m=1\) and \(N_L=\exp(-\pi^2t/L^2)\).

**Bulk threshold and boundary singular modes.** With \(A(s)=\int_0^s a(r)\,dr\), \(\bar A=t^{-1}\int_0^t A(s)\,ds\), put
\[
 B=\exp\!\left[-\int_0^t(A(s)-\bar A)^2\,ds\right].
\]
The free whole-line norm is \(B\). For either wall type the essential spectrum of \(S_+^{a\,*} S_+^a\) is exactly \([0,B^2]\), and \(\|S_+^a\|_{\rm ess}=B\). Every singular value above \(B\) is attained with finite multiplicity and its right singular vectors decay exponentially in \(L^2\) away from the wall. If \(m>B\) at fixed \(a,t\), then (1) strengthens to
\[
 |N_L-m|\le C(a,t,\mathrm{wall\ types})\,e^{-c(a,t,\mathrm{wall\ types})L}
 \qquad(L\ \text{sufficiently large}), \tag{2}
\]
and the finite-channel maximizing right singular vectors are uniformly exponentially localized near the two walls. This last estimate is for a fixed positive gap; it is not a derivative expansion near the essential spectrum edge.

**Application to passive scalar.** In \((x,y)\in\mathbb T\times(0,H)\), scalar diffusivity \(\nu>0\), given affine flow \(u=(\gamma(s)y,0)\), a nonzero horizontal Fourier mode \(k\) solves this equation after factoring out ordinary horizontal heat decay \(e^{-\nu k^2 T}\). For an arbitrary reference \(\gamma_0>0\), set
\[
 \ell=(\nu/(|k|\gamma_0))^{1/3},\qquad
 \tau=(\nu k^2\gamma_0^2)^{-1/3},\qquad
 L=H/\ell,\quad t=T/\tau,\quad
 a(r)=\mathrm{sgn}(k)\,\gamma(\tau r)/\gamma_0.
\]
Thus \(t/L^2=\nu T/H^2\); at fixed enhanced time \(t\) and decreasing viscosity the wall widths tend to infinity and (1) controls exactly the limiting worst-case mixing on the finite physical channel. This is an input-to-output operator norm, not a decay statement about every initial mode, a no-slip hydrodynamic stability estimate, or nonlinear Navier–Stokes. The two real Fourier partners give the identical norm if the scalar is required to be real.

**Known baselines and prior art.** The constant-shear whole-space Kelvin propagator and its cubic-in-time optimized decay are classical; see the Kelvin formula in [Stability threshold of the two-dimensional Couette flow in the whole plane](https://www.sciencedirect.com/science/article/pii/S0022123625004537) and [Grebenkov–Helffer–Henry, arXiv:1603.06992, Appendix B](https://arxiv.org/pdf/1603.06992). The imaginary Airy operator has known wall eigenfunctions, spectrum and asymptotic decay estimates: [Stoller–Happer–Dyson, 1991](https://doi.org/10.1103/PhysRevA.44.7459), [Helffer–Sjöstrand, §3.1](https://arxiv.org/abs/1001.4171) and [Grebenkov–Helffer, 2018](https://doi.org/10.1137/16M1088387). The latter paper and [Averseng–Frantz–Hérau–Raymond, arXiv:2512.20202](https://arxiv.org/abs/2512.20202) study generator eigenfunctions, quasimodes and their spatial localization, rather than singular directions of a general modulated finite-time evolution. [Herberthson et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5367081/) directly compute restricted diffusion under time-varying gradients and include gradient reversal between reflected walls, with a fixed preparation/readout magnetization, not the present worst-case norm theorem. [Quantitative hydrodynamic stability for Couette flow on unbounded domains with Navier boundary conditions](https://arxiv.org/abs/2404.02412) proves enhanced-dissipation estimates on several geometries for a different hydrodynamic problem. We make no priority claim for the Kelvin transform, Airy modes, echo adjoint product, or general use of localization. The prospective addition is the explicit amplitude-uniform finite-channel norm comparison, the exact bulk essential-spectrum edge for arbitrary histories and finite-time singular-vector localization, rather than spectral results for a time-independent generator.

## 1. Path bounds and the finite-width inequality

Use the Feynman–Kac representation for the heat process with generator \(\partial_x^2\), killed at D endpoints and reflected at N endpoints. The phase multiplying each surviving path starting at \(x\) is
\[
 \exp\{-i\int_0^t a(t-r)X_r\,dr\}.
\]
The time reversal in this row-kernel convention follows by applying the time-ordered products from the row variable backwards. Its modulus is one. The formula follows for bounded spatial potentials by the heat product formula or Feynman–Kac; spatially truncate the linear potential for a line or half-line and pass to the limit by pathwise boundedness over finite time and heat domination. The row and column estimates below also hold for the adjoint, whose phase uses the reversed history with opposite sign. Every propagator is a contraction by domination by the substochastic symmetric heat kernel \(h_I\).

Two facts suffice.

1. **Second spatial moment.** The real heat process on a finite interval can be realized by folding \(x+\sqrt2 W_r\) using the triangular periodic fold of the real line, then killing on first hit of a D endpoint. On a half-line use the corresponding single absolute-value fold; on the line no fold is needed. Each fold is 1-Lipschitz and fixes the starting spatial point. Therefore \(\int |x-y|^2h_I(t;x,y)\,dy\le 2t\) uniformly in \(x\). By symmetry the identical column bound holds.
2. **Replacing a distant boundary.** Suppose two model domains have identical walls wherever the starting point is closer than \(d>0\) to any wall on which the domains or conditions disagree. Couple their Brownian paths by identical reflection/killing up to the first such boundary hit. On identical paths the phases coincide. Starting within the common domain, the total variation distance of the phase-weighted endpoint measures on any common output set is at most
\[
 2\,\mathbb P(\sup_{0\le r\le t}|\sqrt2W_r|\ge d)\le 8e^{-d^2/(4t)} .
\]
Before divergence the folded path travels at least \(d\) in physical displacement to hit a disagreement wall; the 1-Lipschitz fold gives the Brownian comparison even near an unchanged reflecting wall. Killing before divergence makes both measures zero for that path. The reflection principle and a Gaussian Chernoff bound give the last inequality. The same estimate holds for columns by reversing/conjugating time order.

Let real, nonnegative cutoffs \(\chi_j\) satisfy \(\sum_j\chi_j(x)^2=1\), at most two nonzero at each \(x\), with \(\|\chi(x)-\chi(y)\|_{\ell^2}\le K|x-y|/L\). On each support use a model domain \(I_j\) agreeing with all nearby walls, at distance at least \(d\) from every disagreement boundary. Operators with \(\chi_j\) on each side are understood via zero extension into the model domain and restriction back. The map \(Jf=(\chi_j f)_j\) is an isometry. Thus \(\|\sum_j\chi_j S_{I_j}^a\chi_j\|=\|J^*\operatorname{diag}(S_{I_j}^a)J\|\le\sup_j\|S_{I_j}^a\|\).

To bound the approximation error, write its integral kernel as
\[
 (1-\textstyle\sum_j\chi_j(x)\chi_j(y))K_I(x,y)
 +
 \textstyle\sum_j\chi_j(x)\chi_j(y)(K_I(x,y)-K_{I_j}(x,y)).
\]
The first prefactor is \(\|\chi(x)-\chi(y)\|_{\ell^2}^2/2\); using the second moment bounds and Schur's test its operator norm is \(\le K^2 t/L^2\). For the second term the row bound is \(2\cdot8e^{-d^2/(4t)}\) (at most two nonzero weights at the row point and each cutoff \(\le1\)); the same column bound follows by starting at the column point with the adjoint. Thus the localization error is at most \(K^2t/L^2+16e^{-d^2/(4t)}\). Countably many supports pose no convergence issue, by local finiteness and the same row/column estimates.

Set \(h=L/2-d\). Apply the estimate in both directions.

* **Finite channel to half-lines.** On \((0,L)\), use \((\chi_l,\chi_r)=(\cos\theta,\sin\theta)\) with \(\theta=0\) through \(d\), linear of slope \(\pi/[2(L-2d)]\) up to \(L-d\), and thereafter \(\theta=\pi/2\). Use half-lines to the right of 0 and to the left of \(L\), respectively. Each cutoff is at least \(d\) from the omitted wall. The localization bound gives \(N_L-m\le \pi^2 t/[4(L-2d)^2]+16e^{-d^2/(4t)}\).
* **Half-line to finite channels.** On the left-wall half-line \((0,\infty)\), use the model cells \((nh,nh+L)\) with centers \(c_n=L/2+nh\), \(n=0,1,\dots\). Set \(\chi_0=1\) through \(c_0\). On each \([c_n,c_{n+1}]\), use \(\chi_n=\cos[(\pi/2)(x-c_n)/h]\) and \(\chi_{n+1}\) the corresponding sine, with all other cutoffs zero. The first cell shares the real wall at 0; its support ends at \(c_1=L-d\). Every later support lies within \(h\) of its own model center and at least \(L/2\) from the real 0 wall, hence is at distance \(\ge d\) from all disagreement walls. Use the same two ordered endpoint conditions in every model cell. Their norms all equal \(N_L\) by translation covariance (the only extra factor is the scalar phase \(e^{-inh A(t)}\)). The cutoff-vector speed is \(\pi/(2h)=\pi/(L-2d)\), so \(N_{\rm left}-N_L\le \pi^2t/(L-2d)^2+16e^{-d^2/(4t)}\). For the right wall, reflect the entire problem, swap the wall type ordering, and restore the sign of \(a\) using complex conjugation; the same bound holds for \(N_{\rm right}\).

These give (1). Its stated sharp leading constant follows since \(d/L\to0\) and \(e^{-d^2/(4t)}=(L/\sqrt t)^{-3}=o(t/L^2)\) for the prescribed \(d\). The unavoidable lower coefficient follows from the first Dirichlet eigenvalue at zero shear. No comparison to whole-line norm is needed for this two-sided estimate.

## 2. Exact essential threshold

On the line, Fourier transform in \(x\) gives the multiplier and translation
\[
 \widehat{S_\mathbb R^a(t)f}(\xi)
 = e^{-\int_0^t(\xi+A(t)-A(r))^2\,dr}\,\widehat f(\xi+A(t)).
\]
Completing the square yields \(B\). Moreover \(S_\mathbb R^{a\,*}S_\mathbb R^a\) is the Fourier multiplier \(B^2 e^{-2t(\xi-\bar A)^2}\), so the line positive operator has spectrum \([0,B^2]\).

Write \(g_t(x)=(4\pi t)^{-1/2}e^{-x^2/(4t)}\). Let \(P\) denote spatial restriction to the positive half-line combined with the corresponding zero-extension in compressions, and \(C=P S_\mathbb R^a P\). For \(x,y>0\), the kernel of \(D=S_+^a-C\) satisfies
\[
 |D(x,y)|\le k\,g_t(x+y),\quad k=1\ (D),\quad k=3\ (N). \tag{3}
\]
Indeed for the absorbing wall subtract only paths that hit 0 and finish at positive \(y\), with unweighted density exactly \(g_t(x+y)\) by reflection. For the reflecting wall compare \(|x+\sqrt2W|\) with free Brownian motion; unhit paths agree. The surviving difference involves both free paths hitting 0 and finishing at positive \(y\), whose two phase factors differ by at most two, and free paths ending at negative \(-y\), re-entering the reflected measure with density \(g_t(x+y)\). More explicitly, conditional on a first hit of zero at time \(r\), the free endpoint is symmetric and the densities on either side are both \(g_{t-r}(y)\); on the + endpoint the two path phases need not agree. This proves \(k=3\) for N without assuming any phase cancellation. Equality of phases is asserted only on nonhitting paths. Either estimate makes \(D\) Hilbert–Schmidt, since \(\int_0^\infty\int_0^\infty g_t(x+y)^2\,dx\,dy=1/(4\pi)\).

Consequently \(S_+^{a\,*}S_+^a-C^*C\) is compact. The cross-operator \((1-P)S_\mathbb R^a P\) is also Hilbert–Schmidt by heat domination with kernel \(g_t(|x|+y)\) on \(x<0,y>0\). Thus \(C^*C\) differs compactly from the compression \(P S_\mathbb R^{a\,*}S_\mathbb R^a P\). The latter is bounded between 0 and \(B^2 I\). To see that every value of the line multiplier, and also 0 by closure, lies in its essential spectrum, take increasing-width smooth compactly supported near-monochromatic packets at a frequency attaining the value and translate each arbitrarily far inside the positive half-line. Before compression the normalized residuals tend to zero by the multiplier formula. The packets can be translated to tend weakly to zero; compression preserves the vanishing-residual bound because their inputs lie in the half-line. Weyl's criterion and compact invariance now prove the claimed essential spectrum. In particular \(N_\pm\ge B\) and strictly larger singular values produce discrete eigenvalues of the associated positive operator. This is a statement about singular modes, not the eigenfunctions of the imaginary-potential generator.

## 3. Exponential localization in the gapped regime

The following weighted argument covers both a half-line and a finite interval. For the domain \(I\), let \(r_I(x)\) be distance to its finite boundary, so \(r_+(x)=x\) and \(r_{I_L}(x)=\min(x,L-x)\). Let \(C_I\) be the line propagator compressed to \(I\), \(D_I=S_I^a-C_I\), and \(T_0=C_I^*C_I\). Then \(\|T_0\|\le B^2\) and its kernel is dominated by \(g_{2t}(x-y)\) by Gaussian convolution.

The half-line Gaussian bound for \(D_I\) is (3) with the safe constants. On a finite interval the same path argument gives
\[
 |D_{I_L}(x,y)|\le 3\sum_{n\in\mathbb Z}g_t(x+y+2nL).
\]
In detail: the free contribution from paths hitting one or both walls and ending in the interval is at most \(g_t(x+y)+g_t(2L-x-y)\), by taking the union and reflection at each wall. The bounded-domain surviving contribution is at most that of the fully reflecting process after it has hit some wall. That density is \(h^{NN}_{I_L}-h^{DD}_{I_L}=2\sum_n g_t(x+y+2nL)\) by folding and the classical image kernels. Homogeneous mixed killing only removes paths from this latter contribution. Thus both the row and column estimates follow, and there is no independence assumption about double hits.

For any fixed \(\alpha>0\), use the weights \(w_R(x)=\exp(\alpha\min(r_I(x),R))\). The Gaussian bounds show \(\|w_R D_I\|\) and \(\|w_R D_I^*\|\) are bounded in \(R\) and in \(L\): for the interval, \(r_I(x)\le |x+y+2nL|\) for all \(n\) and \(x,y\in I\). In Schur's test the integral over \(y\) and sum over \(n\) is at most a whole-line integral of \(g_t(v)e^{\alpha|v|}\), since these \(v\)-intervals of length \(L\) are disjoint. The column follows the same way, and the analogous half-line estimate has just one image term. The weight is multiplicatively Lipschitz, \(w_R(x)/w_R(y)\le e^{\alpha|x-y|}\). Consequently
\[
 \|w_RT_0w_R^{-1}-T_0\|
 \le G_t(\alpha):=\int_{\mathbb R}g_{2t}(v)(e^{\alpha|v|}-1)\,dv
 \longrightarrow0 \quad (\alpha\downarrow0).
\]
Likewise \(\|w_R C_I^*w_R^{-1}\|\) is uniformly bounded. With
\(\Delta_I=S_I^{a\,*}S_I^a-T_0=C_I^*D_I+D_I^*C_I+D_I^*D_I\),
these estimates show \(\|w_R\Delta_I\|\) uniformly bounded. If \(S_I^{a\,*}S_I^a f=\lambda f\), \(\lambda>B^2\), choose \(\alpha\) with \(G_t(\alpha)<\lambda-B^2\). Solve the weighted eigenvalue equation by the convergent resolvent series about this norm-smaller operator. It bounds \(\|w_R f\|\) uniformly, so monotone convergence yields \(\|e^{\alpha r_I}f\|<\infty\). At a fixed uniform spectral gap the bound is uniform in \(L\).

For (2), the \(\lambda_L=N_L^2\) maximizing finite eigenvalue is attained by compactness of the finite heat-dominated propagator, and (1) gives for all large \(L\) a uniform gap above \(B^2\) if \(m>B\). Its unit input \(f_L\) therefore loses at most \(C e^{-\alpha L/3}\) in norm when its mass on \((L/3,2L/3)\) is dropped. Each remaining input piece is supported in the left or right outer third. Split the output into left and right halves. The blocks that cross from one input outer third to the opposite output half have operator norm \(O(e^{-L^2/(144t)})\) by the heat-domination displacement bound and Schur's test (distance at least \(L/6\)). The two diagonal blocks differ from their associated one-wall restrictions by \(O(e^{-L^2/(16t)})\) by the distant boundary comparison: both the row and column supports are at least \(L/2\) from the disagreement wall. Orthogonal halves then give \(N_L\le m+O(e^{-\alpha L/3}+e^{-L^2/(144t)})\).

Conversely choose a normalized maximizing singular vector of the half-line attaining \(m>B\). Truncate it at \(L/3\), compare its finite-channel output on the near half by the same distant-boundary bound, and discard the half-line output beyond \(L/2\). The lost input norm is exponentially small by localization and the lost output of the truncated input is Gaussian small by the distance \(L/6\) and heat domination. This proves \(N_L\ge m-O(e^{-\alpha L/3}+e^{-L^2/(144t)})\), completing (2).

## 4. Exact wall-dominance control and reversal

At constant \(a=1\), the whole-line Kelvin norm is \(B=e^{-t^3/12}\). A negative real zero \(z=-\beta\) of \(\operatorname{Ai}\) for a D wall (or of \(\operatorname{Ai}'\) for an N wall) gives the square-integrable eigenfunction \(\operatorname{Ai}(e^{i\pi/6}x-\beta)\) of the forward generator, with eigenvalue \(-e^{i\pi/3}\beta\). This is the known Airy eigensolution, verified by the Airy equation and the decay in the indicated complex sector. Thus \(\|S_+^1(t)\|\ge e^{-\beta t/2}>B\) whenever \(t^2>6\beta\). This is a sufficient gap condition, not a claim that it is necessary or that this eigenfunction is the top singular vector.

The offline rational replay certifies a Neumann channel with \(t=3,L=128\). The Airy derivative has a negative zero between \(-1.02\) and \(-1.00\), giving the half-line bound \(>e^{-1.53}>0.216\). With \(d=12\), the (1) error is \(<0.003\), so the finite norm is \(>0.213\), while \(B=e^{-9/4}<0.106\). The squared-norm statistic is still more sharply separated. In fact for a second half with sign-reversed, time-reversed shear, \(a(t+r)=-a(t-r)\), the second half is \(S_I^{a\,*}\), as is seen by adjointing the time-ordered evolution. The refocusing/echo evolution is therefore \(S_I^{a\,*}S_I^a\), with norm exactly \(N_I^2\). An integer number of identical echoes has norm \(N_I^{2j}\) by the spectral theorem. The example thus certifies a finite-channel echo norm more than \((0.213/0.106)^2>4.03\) times the infinite-line prediction, even though the inviscid flow reverses its displacement exactly. The horizontal diffusivity factor is common and cancels in this comparison. This ratio is for worst-case scalar survival, not a prediction of a scalar signal for a fixed preparation/measurement or a demonstration of increased kinetic energy.
