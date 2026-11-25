# What finite-frequency heat response determines about static restoring feedback

## Frozen question, assumptions and result

Suppose a stable, linear thermal network is driven and measured at the same surface temperature port. All heat capacities and reciprocal internal conductances are positive where present. Internal heat transfer conserves total heat, and only the surface loses heat to the exterior, through one unknown restoring coefficient \(\lambda>0\). The surface heat capacity \(C_0>0\) and total *connected* subsurface heat capacity \(C>0\) are independently known. At \(m\) distinct positive angular frequencies \(\omega_1,\ldots,\omega_m\) one observes exact complex forcing-to-temperature admittances \(Y_j=\widehat F(\omega_j)/\widehat T_0(\omega_j)\) in stationary harmonic response. There is no separate measurement of the radiative or total energy imbalance. The equilibrium (DC) temperature response per unit sustained forcing is \(1/\lambda\).

The **largest possible diameter** of the compatible restoring-coefficient set, over all feasible measured vectors in this model, is exactly
\[
 W_{\rm all}=\frac{C}{2\left(\sum_{j=1}^m\omega_j^{-1/2}\right)^2}. \tag{1}
\]
There is an explicit affine estimator with error at most \(W_{\rm all}/2\); no estimator driven by only these exact observations has a smaller worst-case absolute error. The supremum is approached by pairs of finite, strictly positive connected star networks which have the *same* values of every reported complex admittance. Arbitrarily slow/fast limiting rates are used for sharpness, not treated as physically present.

If positive thermally visible transport rates are known not to exceed \(b\), the exact worst-data diameter improves to
\[
 W_{[0,b]}=\frac{C b}{2\left(\sum_{j=1}^m\sqrt{1+\sqrt{1+(b/\omega_j)^2}}\right)^2}. \tag{2}
\]
For a **known strictly positive finite interval** \(0<a<b<\infty\), the analogous exact formula is \(W_{[a,b]}=2C/\kappa_{a,b}\), with the elementary rational function and derivative for \(\kappa_{a,b}\) below. There the equality pair uses only rates in the given interval, with no singular limits. Each of the two catalogs needs at most \(m+1\) reservoirs. These formulas concern rates of the interior temperature dynamics with the surface port fixed, not eigenrates of the full free network.

**Scope.** Nothing here infers an actual climate feedback or establishes a range for Earth. The transfer law does not encompass a spatially distributed forcing, noncollocated temperature, pattern-dependent or time-varying feedback, transport with advection-induced nonreciprocity, or unmeasured internal forcing. With forcing, radiative imbalance and surface temperature measured as exact same-model independent channels, the direct balance identifies \(\lambda\), and this ambiguity does not apply. Taking \(1/\lambda\) gives no uniform finite sensitivity interval near \(\lambda=0\) without an independent positive lower bound. If the model and data establish \(\lambda\ge\lambda_*>0\), then the width in \(1/\lambda\) is at most \(W/\lambda_*^2\).

**Primary baselines.** [Geoffroy et al. 2013](https://doi.org/10.1175/JCLI-D-12-00195.1) give the standard calibrated two-layer response; [Fredriksen–Rypdal 2017](https://doi.org/10.1175/JCLI-D-16-0877.1) discuss multiple relaxation scales in multibox EBMs. [Lovejoy 2021, §3.3](https://doi.org/10.5194/esd-12-469-2021) explicitly treats complex thermal admittance and periodic surface forcing for a half-order heat equation. [Bastiaansen–Ashwin–von der Heydt, climate response and late tipping points](https://arxiv.org/abs/2207.06110) and [Finite-Response Complementarity, 2026](https://arxiv.org/abs/2608.07356) distinguish finite responses from DC sensitivity in other contracts. [Armour 2017](https://doi.org/10.1038/nclimate3278) explains why time-varying climate feedbacks and ocean heat-uptake efficacy take real climates outside this simplified reciprocal model. [Ziegeler, Nolte and Schweizer 2021](https://doi.org/10.3390/en14217068) describe established Foster/Cauer thermal-network response and reconstruction. Rational Chebyshev/Markov methods are classical; see [Borwein–Erdélyi–Zhang 1994](https://www.cecm.sfu.ca/~pborwein/PAPERS/P90.pdf) and [Kalmykov–Nagy–Totik](https://arxiv.org/abs/1610.06706). The claim is this explicit sharp identification width and a physically realizable indistinguishability construction, not the invention of thermal transfer functions, the general inverse problem, or rational alternation.

## 1. Exact physical reduction to positive transport-rate capacity

Label surface temperature \(T_0\) and hidden temperatures \(u\). Write the diagonal hidden capacity matrix \(D_c\), the vector \(g\) of direct surface conductances, and \(K=L_{\mathrm{inner}}+\mathrm{diag}(g)\), with \(L_{\mathrm{inner}}\) the symmetric internal graph Laplacian. Connectedness to the surface makes \(K\) positive definite even if some \(g_i=0\); \(K\boldsymbol 1=g\). The surface equation reads \(C_0\dot T_0=F-\lambda T_0-\boldsymbol1^\top g\,T_0+g^\top u\), and \(D_c\dot u=g T_0-Ku\). Hence
\[
 Y(s)=\lambda+C_0s+\boldsymbol1^\top g-g^\top(sD_c+K)^{-1}g.
\]
Diagonalize \(M=D_c^{-1/2}KD_c^{-1/2}\) with rates \(r_\ell>0\), and \(v=D_c^{1/2}\boldsymbol1\). Then \(D_c^{-1/2}g=Mv\), \(\boldsymbol1^\top g=v^\top Mv\), and
\[
 Y(s)=\lambda+C_0s+\int \frac{r s}{r+s}\,d\mu(r),
 \qquad \mu=\sum_\ell |\langle v,e_\ell\rangle|^2\delta_{r_\ell},\qquad
 \mu((0,\infty))=\|v\|^2=C. \tag{3}
\]
Invisible eigenvectors have zero mass and impose no restriction. Conversely a finite positive rate measure of mass \(C\) gives a connected star: an atom of mass \(c_\ell\) at \(r_\ell>0\) is a capacity \(c_\ell\) connected only to surface by conductance \(c_\ell r_\ell\). Every extremal construction below is therefore a real passive thermal network, not a signed or independently tunable complex response.

The **same upper bounds also apply to a finite continuous diffusion column** of heat capacity density \(c(y)>0\), positive conductivity, fixed driven surface temperature and insulating bottom. For uniformly positive bounded coefficients on finite depth, let \(A=-c(y)^{-1}\partial_y(k(y)\partial_y)\) be the positive diffusion operator with homogeneous top Dirichlet and bottom Neumann boundary conditions, self-adjoint on \(L^2(c\,dy)\). For the harmonic interior state \(u\), write \(w=T_0-u\), which has homogeneous boundary data and satisfies \(\dot w+Aw=\dot T_0\,\boldsymbol1\). The flux into the column equals the time derivative of its total heat by conservation. Thus its term in \(Y(s)\) is \(s\langle\boldsymbol1,A(s+A)^{-1}\boldsymbol1\rangle\), with (3) by the spectral theorem and \(C=\int c(y)\,dy\). The rate support is ordinarily unbounded above, so (1) is the applicable bound absent further restrictions. Sharpness over all thermal networks does not assert sharpness within this narrower diffusion-column family.

Define real data
\[
 R_j=\Re Y(i\omega_j)=\lambda+\int h_j\,d\mu,\qquad
 Q_j=\frac{\Im Y(i\omega_j)}{\omega_j}-C_0=\int f_j\,d\mu,
\quad
 h_j(r)=\frac{\omega_j^2r}{r^2+\omega_j^2},\quad f_j(r)=\frac{r^2}{r^2+\omega_j^2}.
\]
Here \(\omega_j\) is angular frequency, not cycles per unit time. Only \(R_1\) contains a new free common level: differences \(R_j-R_1\) determine \(\int (h_j-h_1)\,d\mu\).

## 2. What the observations freeze: one missing derivative functional

Let \(D(r)=\prod_{j=1}^m(r^2+\omega_j^2)\) and \(\mathcal V=\{p(r)/D(r):\deg p\le 2m\}\) over the reals. Let \(L(F)=F'(0)\). The known-moment space is exactly
\[
 \mathcal U=\mathrm{span}\{1,f_1,\ldots,f_m,h_2-h_1,\ldots,h_m-h_1\}=\ker L. \tag{4}
\]
Indeed the first \(m+1\) functions give all even numerators: multiplying by \(D\) and evaluating at \(r^2=-\omega_j^2\) first isolates each \(f_j\), then the constant. The \(m-1\) odd differences are independent by evaluation at their distinct poles (or at the corresponding complex roots of \(D\)); each has derivative zero at 0 because \(h_j'(0)=1\). Their dimension is the dimension of odd numerators with no linear term. Since \(D'(0)=0\), these two subspaces fill exactly the stated kernel.

Let the closed allowed rate support (including endpoint limits if needed) be \(K\). Put
\[
 \kappa_K=\max_{\|F\|_{\infty,K}\le1,\ F\in\mathcal V}|F'(0)|.
\]
For an extremizer \(F_*\) with \(F_*'(0)=-\kappa_K\), \(g_*=h_1+F_*/\kappa_K\in\mathcal U\). Then the observable estimator
\[
 \widehat \lambda=R_1-\int g_*\,d\mu
\]
has absolute error at most \(C/\kappa_K\). Two networks with the same data must therefore differ by at most \(2C/\kappa_K\). This also gives a specific sufficient condition for target precision without assuming a number of boxes or selecting an arbitrary branch of a nonlinear fit.

## 3. Exact rational extremizer and physical minimax pair

### Unrestricted positive rates

For \(u=\sqrt r\), let
\[
 Q(u)=\prod_{j=1}^m\big(u^2-\omega_j-i\sqrt{2\omega_j}\,u\big),\qquad
 F_*(r)=\Re\{Q(u)/\overline{Q(u)}\}.
\]
Each squared factor modulus is \(r^2+\omega_j^2\). Thus \(F_*=\Re Q^2/D\), whose numerator is a real even polynomial in \(u\) (a polynomial in \(r\) of degree at most \(2m\)). As \(u\) rises from zero to infinity, each factor has its argument increasing strictly from \(-\pi\) to zero: its derivative is \(\sqrt{2\omega_j}(u^2+\omega_j)/|q_j(u)|^2>0\). Consequently \(F_*\) takes alternating values \(+1,-1,\ldots,+1\) at exactly \(2m+1\) ordered points \(r_0=0<r_1<\cdots<r_{2m-1}<r_{2m}=\infty\). Its expansion at zero is
\[
 F_*(r)=1-4\left(\sum_j\omega_j^{-1/2}\right)^2 r+O(r^2).
\]
This follows by taking the cosine of twice the accumulated phase, whose \(u\)-slope at zero is \(\sum_j\sqrt{2/\omega_j}\).

To prove optimality for the specified derivative, interpolate \(F\in\mathcal V\) at these nodes. Write \(L(F)=\sum_{j=0}^{2m}w_j F(r_j)\). The rational cardinal function for a finite node \(r_j\) is its usual Lagrange polynomial in the \(2m\) finite nodes, times \(D(r_j)/D(r)\); the infinity node is \(\prod_{j=0}^{2m-1}(r-r_j)/D(r)\), with monic \(D\). The derivative weights have signs \(\mathrm{sgn}(w_j)=(-1)^{j+1}\): for \(j>0\) differentiate the numerator factor \(r-r_0=r\) at zero and count its remaining factors, for \(j=0\) the derivative is the negative sum of reciprocals of all finite nodes other than zero, and for infinity the derivative is negative by the odd number of nonzero finite nodes. The factors \(D(r_j)/D(0)>0\) do not change signs, and \(D'(0)=0\). All \(w_jF_*(r_j)\) have the same negative sign. Thus
\[
 |L(F)|\le\big(\textstyle\sum_j|w_j|\big)\|F\|_{\infty,K},\qquad
 \sum_j|w_j|=-L(F_*)=\kappa_K,
\]
proving the extremizer without assuming an alternation theorem for the missing derivative.

**Sharpness of the thermal inference.** The signed measure \(\sum_j w_j\delta_{r_j}\) integrates any \(F\in\mathcal V\) to \(L(F)\). Its total mass is zero and variation \(\kappa_K\). Scale its positive/negative parts separately by \(2C/\kappa_K\), giving two positive measures each of mass \(C\), agreeing on \(\mathcal U\), and differing in their \(h_1\) moment by \(2C/\kappa_K\). Their \(h_j\) moment difference is the same for every \(j\), so use opposite restoring-coefficient difference of that amount to obtain identical complex admittance at every frequency, with both coefficients strictly positive by choosing a positive offset. To avoid endpoint rates 0 and \(\infty\), shift them slightly into positive finite domain and recompute the interpolation weights. These weights vary continuously under this change: \(\mathcal V\) has a continuous basis in the compactified coordinate \(r/(1+r)\), the node-value matrix at the limiting nodes is invertible, and its derivative functional is fixed. Thus total variation converges to \(\kappa_K\), while moment equality and positive-star realizability hold exactly at every finite stage. The limiting pair proves both exact worst-data diameter (1) as a supremum and the minimax estimator error of half that width.

**A single fixed dataset can express the full supremal width.** Fix the positive star on the \(m\) odd nodes and its positive restoring coefficient (chosen sufficiently large). In the limiting even-node star, keep the endpoints at prescribed perturbed positive finite rates, and solve for the \(m+1\) masses and \(m-1\) interior even positions so that all \(2m\) moments in \(\mathcal U\) match the fixed odd star. The Jacobian in these \(2m\) variables is nonsingular at the limiting star. A vector in the left nullspace would be a function in \(\mathcal U\) vanishing at all even nodes and with zero derivative at the interior even nodes (the corresponding masses are strictly positive). For \([0,\infty]\), its numerator has degree at most \(2m-1\) because it vanishes at infinity, and has the factor \(r\prod_{j=1}^{m-1}(r-r_{2j})^2\), already of degree \(2m-1\). The only possible numerator is a multiple of that factor, whose derivative at zero is nonzero unless the multiplier is zero; membership in \(\mathcal U\) forces it to be zero. For \([0,b]\) the analogous product has the additional factor \(r-b\), the full degree \(2m\), and the same argument applies. The implicit function theorem (using inverse rate as a coordinate at infinity) supplies positive finite even catalogs after moving all forbidden endpoints inward. Adjust only their restoring coefficients to match the full common real levels. Thus one finite positive reference star's exact measured vector has the stated diameter as a supremum, rather than merely a limit of changing data vectors. The finite star on the odd nodes is fixed throughout.

### Bounded positive rates

For the closed interval \([a,b]\), \(0\le a<b<\infty\), reparameterize with \(r=(a+b u^2)/(1+u^2)\), \(u\ge0\). Write
\[
 \alpha_j=\sqrt{a^2+\omega_j^2},\quad \beta_j=\sqrt{b^2+\omega_j^2},\quad
 \zeta_j=\sqrt{2(ab+\omega_j^2+\alpha_j\beta_j)},\quad
 Q_{a,b}(u)=\prod_j(\beta_j u^2-\alpha_j-i\zeta_j u).
\]
Now \(|Q_{a,b}(u)|^2=(1+u^2)^{2m}D(r)\). The same expression \(\Re(Q_{a,b}/\overline{Q_{a,b}})\), continued as a rational function of \(r\), is in \(\mathcal V\): its even numerator in \(u\) is transformed by \(u^2=(r-a)/(b-r)\) and the factor \((b-r)^{2m}\) cancels the substitution denominator. The phase argument gives \(2m+1\) alternating extrema including both finite endpoints. If \(a>0\), ordinary Lagrange differentiation at the external point 0 still gives \(\mathrm{sgn}(w_j)=(-1)^{j+1}\): the sign of the cardinal numerator at zero is \((-1)^j\) and its logarithmic derivative at zero is the strictly negative sum of reciprocals of all *other* nodes. When \(a=0\), differentiate the factor \(r\) as before. Thus the explicit rate-interval constant is
\[
 \kappa_{a,b}=-F_{a,b}'(0)>0.
\]
For \(a>0\) all alternating nodes are finite positive and the signed-interpolation construction attains the exact diameter with two finite stars. When \(a=0\), move the zero endpoint slightly inwards and take a limit as above. Setting \(a=0\), expansion at \(u=0\) uses \(r=b u^2+O(u^4)\) and each phase slope \(\zeta_j/\omega_j\), giving
\[
 \kappa_{0,b}=\frac{4}{b}
       \left(\sum_j\sqrt{1+\sqrt{1+(b/\omega_j)^2}}\right)^2.
\]
This is (2). For \(b\ll \min_j\omega_j\), it gives \(W_{[0,b]}\sim Cb/(4m^2)\); as \(b\to\infty\) it approaches (1). No assumption of small \(b\) is needed for the exact bound.

## 4. Checks, use, and what noise changes

* **One frequency.** With \(q=Q_1/C\), the real transport moment is \(\omega_1 C\,\mathbb E\sqrt{x(1-x)}\) under \(x=f_1(r)\) and \(\mathbb E x=q\). Its closed range is \(0\) to \(\omega_1 C\sqrt{q(1-q)}\), by the concave upper bound and endpoint mixtures approaching the lower bound. The largest width \(C\omega_1/2\) agrees with (1). This simple single-frequency control is not the new result. If frequencies duplicate, count each distinct frequency just once; an exact duplicate adds no information.
* **Two-frequency physical indistinguishability.** In consistent normalized units use \((\omega_1,\omega_2)=(1,4)\), \(C=1\), \(C_0=1/2\) and \(r_j=(1/1000,469/1000,2,8531/1000,1000)\). Compute exact Lagrange derivative weights at these nodes, split by sign, normalize each star to \(C=1\), and offset their feedbacks about the positive reference \(\lambda=2\) so both real levels coincide. The exact replay checks both complex admittances match and \(\Delta\lambda>199/900\), more than \(99.5\%\) of \(2/9\), the theorem's sharp supremum on all positive rates. No zeros or infinite conductances occur in these two stars.
* **Affine robust certificate.** Expand \(g_*=a_0+\sum_j a_j f_j+\sum_{j\ge2}b_j(h_j-h_1)\) by a finite linear system. Then
  \[
  \widehat\lambda=(1+\textstyle\sum_{j\ge2} b_j)R_1-\sum_{j\ge2}b_j R_j
      -a_0 C-\sum_j a_j(\Im Y_j/\omega_j-C_0).
  \]
  This is exact-data minimax; with bounded calibration errors one adds the ordinary signed-coefficient norm for the measured error set to \(C/\kappa_K\), and also capacity uncertainties (use an upper bound for true \(C\) in the structural term). For unbounded rates and \((1,4)\), the exact coefficients are \(a_0=1/9\), \(a_1=16/15\), \(a_2=-16/15\), \(b_2=-1/5\); for componentwise admittance errors \(|\delta\Re Y_j|,|\delta\Im Y_j|\le\eta\), with exact capacities, the total deterministic error radius is at most \(1/9+(7/3)\eta\) when \(C=1\). Nearly repeated frequencies can make the coefficient solve ill-conditioned; no claim that clustering frequencies offers the noiseless gain under fixed instrument precision.
* **A deterministic spacing limit under calibration error.** Suppose every tone lies within \(\eta/C\) of one reference tone \(\omega_0>0\), one observes that tone, and each complex admittance has absolute error at most \(\eta\). Take an indistinguishable single-tone pair approaching the width \(C\omega_0/2\). For a positive branch kernel, \(\left|\partial_\omega[r i\omega/(r+i\omega)]\right|=r^2/(r^2+\omega^2)\le1\). With the same \(C_0\) and \(C\), the pair's admittances therefore have separation at most \(2C|\omega-\omega_0|\), after their exact matching at the reference tone. Taking their midpoint independently at each additional tone yields a single allowed noisy record from both. No estimator on such an adversarially calibrated cluster beats the single-tone structural minimax lower bound \(C\omega_0/4\), regardless of the number of accurate-in-form-but-imprecise neighboring tones. This is a lower bound for the stated worst-case calibration model, not for independent unbiased observation noise.
* **Interpretation.** With exact harmonic pairs at \(\omega_j=j\omega\), (1) gives \(W\sim C\omega/(8m)\), not zero feedback uncertainty from any finite number of harmonics. With one known rate ceiling, (2) quantifies precisely why independent transport constraints can narrow the uncertainty. Neither frequency idealization produces an unconditional bound from weather variability or a measurement whose forcing input is unknown.

Run python verify.py for the exact offline controls. The main proof is analytic and neither fits a climate record nor requires data downloads.
