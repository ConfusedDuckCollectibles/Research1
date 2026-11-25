# A source track does not identify a thin lens in wave optics

## Setup and conclusions

Consider scalar coherent wave propagation in the usual weak-field, static, single-plane thin-lens diffraction approximation. Fix the angular and frequency-to-delay scales. The lens-plane position is \(\boldsymbol x=(u,v)\), source-plane offset is \(\boldsymbol y\), projected convergence is \(\kappa=\tfrac12\Delta\psi\), and dimensionless delays are
\[
 S(\boldsymbol x)=\tfrac12|\boldsymbol x|^2-\psi(\boldsymbol x),\qquad
 T_{\boldsymbol y}(\boldsymbol x)=S(\boldsymbol x)-\boldsymbol x\cdot\boldsymbol y+\tfrac12|\boldsymbol y|^2.
\]
In a common delay gauge the amplification is the usual Fresnel integral
\[
 F(w,\boldsymbol y)=\frac{w}{2\pi i}\int_{\mathbb R^2}e^{iwT_{\boldsymbol y}(\boldsymbol x)}\,d^2\boldsymbol x. \tag{1}
\]
Integrals can be defined with \(\Im w>0\) and boundary values on real \(w\). For any differences with compact spatial support appearing below, one directly subtracts the matching tails and obtains ordinary integrals for real \(w\).

1. **Exact straight-track nonidentifiability.** Infinitely many distinct smooth, nonnegative finite-mass lenses can give *exactly* the same calibrated complex amplification for every frequency and every source position on an entire known straight line. They can also agree on the full local arrival-time surface around every geometric image from every source on that line: not just image positions, magnification ratios, or time-delay ratios. A concrete positive-density control below has a strictly local same-point density change exceeding \(1/400\) (more than 1% at that point). No mass sheet, source rescaling, time shift or magnification factor is used.
2. **What the line probes, and a contrasting sufficient aperture.** For compact changes to a potential, the full-frequency response on any open segment of a fixed straight line is exactly a constraint on one joint position–delay distribution; it need not fix the spatial potential. Conversely equality of absolutely calibrated complex data on an open **two-dimensional** source-offset patch at even one nonzero real frequency rules out a nonzero compact potential difference. This is uniqueness, not stable recovery from a small or noisy aperture, and does not assert two-dimensional coverage is minimal.
3. **Circular prior changes the answer.** In circular lenses with a single radial delay well, the full spectrum at exact alignment alone still has infinitely many indistinguishable circular positive finite-mass profiles (even positive decreasing spherical three-dimensional density is possible). The aligned spectrum **plus the second derivative with respect to calibrated transverse source offset**, as a function of frequency, uniquely determines the radial delay profile up to a constant time gauge. This statement survives an unknown per-offset smooth frequency-independent amplitude and time origin, but not general unknown intrinsic waveform distortions. A fixed nonzero small source displacement separates the explicit circular family at real frequencies with a quantitative lower bound.

**Known baseline and limits.** The Fresnel/coarea formula and the sensitivity of diffraction to the full distribution of arrival delays are classical: [Ulmer–Goodman 1995](https://doi.org/10.1086/175422), [Takahashi–Nakamura 2003](https://doi.org/10.1086/377430), [Tambalo et al. 2023](https://doi.org/10.1103/PhysRevD.108.043527), and [Meena 2024 §§2–4](https://arxiv.org/html/2409.00145v2). Meena's explicit example distinguishes a point mass from its particular centered redistribution; nothing here disputes that comparison or the paper's global degeneracy analysis. Section 5 argues more broadly that any local delay-surface degeneracy changes the amplification. The present compact local counterexample identifies a precise exception: differently shaped contours can have identical joint area distributions. [Liesenborgs–De Rijcke 2012](https://doi.org/10.1111/j.1365-2966.2012.21751.x) treats nonuniqueness in geometric imaging; [Villarrubia-Rojo et al., GLoW](https://arxiv.org/abs/2409.04606) and [GWLENS, Nov 2026](https://arxiv.org/html/2611.24621v1) compute wave-optical responses including moving-source settings, not the exact line-track null construction. The equal inverse-branch shear itself is the classical Abel construction for isoperiodic potentials, explicitly presented by [Asorey–Cariñena–Marmo–Perelomov, §2](https://arxiv.org/abs/0707.4465); the physical diffraction null and its complementary lens data are the claims here. No new astronomical observation or claim that a calibrated offset derivative is available for a current event. A source-line orbit approximation, unknown unlensed signal, finite detector band, near-field propagation or time-dependent lens each demands separate analysis.

## 1. An exact null family for an entire straight track

Choose axes so the source track is \((s,0)\), \(s\in\mathbb R\). For an affine line \((s,c)\) with \(c\ne0\), replace \(S\) in this construction by \(S-cv\) (or translate the origin in the normal coordinate). Choose a real compactly supported function \(H(u,t)\). Suppose the preimage of its support under \((u,v)\mapsto(u,S(u,v))\) lies within a compact set where baseline \(\kappa\) has a strictly positive margin, and avoids all stationary transverse fibers \(S_v=0\) and their neighborhoods. These conditions hold in the explicit example below. Set
\[
 \chi_\epsilon(u,v)=(u,v+\epsilon H(u,S(u,v))),\qquad
 S_\epsilon(\chi_\epsilon(u,v))=S(u,v). \tag{2}
\]
For sufficiently small \(|\epsilon|\), \(\partial_v\chi_{\epsilon,v}=1+\epsilon H_t S_v>0\), and the map is a diffeomorphism because each fiber is strictly increasing and equals the identity outside a compact set. Here the diffeomorphism need not preserve lens-plane area. Rather, it preserves the relevant *joint distribution* as follows.

For any smooth test \(f(u,t)\), after the change of variable the difference of the two integrals of \(f(u,S)\) is exactly
\[
 \epsilon\int du\int_{-\infty}^{\infty} f(u,S(u,v))H_t(u,S(u,v))S_v(u,v)\,dv=0. \tag{3}
\]
For each fixed \(u\), the inner integral is the difference of a primitive of \(f(u,t)H_t(u,t)\) at the two fiber ends; \(S\to+\infty\) at both. No branch count is needed, no approximation in \(\epsilon\) is made, and orientation signs in the primitive are essential. The compact difference permits oscillatory tests in (3). With \(f(u,t)=\exp[iw(t-su)]\) and the common \(e^{iws^2/2}\), (3) proves exact equality of (1) for every \(s,w\).

The transformed potential is \(\psi_\epsilon=|\boldsymbol x|^2/2-S_\epsilon\). It equals \(\psi\) outside a compact set, so total excess mass is zero by integrating \(\Delta(\psi_\epsilon-\psi)\) and the exterior deflection is unchanged. Values of \(\Delta S_\epsilon\) vary continuously in the \(C^2\) norm under a small fixed smooth compact perturbation; hence the new convergence stays positive on the changed region and is exactly unchanged outside. Conversely this density produces the prescribed potential: solving the two-dimensional Poisson equation with the same exterior convention differs only by a globally harmonic function vanishing outside, thus by zero. The lens remains finite mass. Distinct shifts of a nontrivial regular contour give distinct profiles; the choices of \(H\) form an infinite-dimensional family.

For every \(s\), \(T_{(s,0),\epsilon}=T_{(s,0)}\circ\chi_\epsilon^{-1}\) because the map fixes \(u\). The critical sets are consequently carried by \(\chi_\epsilon\). Every original critical point has \(S_v=0\), where the map is locally the identity by hypothesis. There are no new critical points, and all old images retain the same neighborhoods of the lens potential, all differential invariants and their absolute delays in this gauge. This uses neither a mass-sheet transformation nor a rescaling of source coordinates. The same construction applies for an affine track with the transverse-critical condition on \(S-cv\). Moreover, when \(|S_v|\) is bounded below on the compact changed region, \(|\partial_v S_\epsilon|\) is likewise bounded below there. Thus both lenses actually agree on all geometric image neighborhoods in an open **strip** of source offsets around the track. Fresnel equality is only asserted on the line.

### What line data identify; why a two-dimensional patch is different

For a compact change, write \(\nu\) for the finite signed difference of the joint pushforwards of \(du\,dv\) by \((u,S_1)\) and \((u,S_0)\), computed over a common compact region containing the change. At real frequency \(w\ne0\), the difference of (1) with common prefactors removed is \(\widehat\nu(ws,-w)\) in the convention \(\widehat\nu(k,\ell)=\int e^{-i(ku+\ell t)}\,d\nu\). Compact support makes this transform real analytic jointly in its two arguments. Any open rectangle in \((s,w)\) with \(w\ne0\) maps to an open two-dimensional patch under \((s,w)\mapsto(ws,-w)\); equality there implies \(\nu=0\) by analyticity and Fourier uniqueness. Thus one need not assume data at all offsets or frequencies for this exact comparison (stability is a separate issue). The fiber construction changes the potential without changing precisely this joint measure.

For different lenses at **one fixed** nonzero frequency \(w\), (1), with common phase and prefactor removed, has difference \(\widehat{e^{iwS_1}-e^{iwS_0}}(w\boldsymbol y)\). The screen difference is compactly supported. Its Fourier transform is real analytic on \(\mathbb R^2\); vanishing on an open 2D patch implies zero everywhere. Hence \(S_1-S_0\) is an integer multiple of \(2\pi/w\) pointwise, and continuity plus equality outside the compact support forces zero globally. This last result assumes absolute complex calibration at that frequency across the entire patch; arbitrary per-offset gain/phase changes would erase the assertion. It yields no quantitative stability of analytic continuation.

**First transverse derivative under a one-interval fibre prior.** With absolutely calibrated data, taking the first derivative in the source coordinate normal to the line, after removing the common prefactors, inserts \(-iwv\) into the integrand. Hence on the same open \((s,w)\) rectangle it determines equality of the second joint measure weighted by \(v\,du\,dv\) for compact differences. Disintegrate the two measures in \(u\); for almost every \(u\) the unweighted sublevel length and first \(v\) moment of every fixed rational delay level agree (take the simultaneous countable intersection), so right continuity extends equality to all levels. If each fibre sublevel is a single bounded interval, its length \(A_0\) and first moment \(A_1\) determine the endpoints \(A_1/A_0\pm A_0/2\) where it has positive length. This determines the continuous delay field (up to changes on null sets before imposing continuity). The derivative is a sufficient extra observable for this class, not a claim of its availability or of unique recovery without a fibre-shape prior.

## 2. Circular alignment and the offset-curvature remedy

For a circular lens write \(x=r^2/2\), so \(r\,dr=dx\), and let \(\Phi(x)=x-\psi(\sqrt{2x})\). The convergence in this coordinate is
\[
 \kappa(x)=1-\Phi'(x)-x\Phi''(x).
\]
Assume a smooth finite-core finite-mass lens with one strict delay minimum and compact interval sublevel sets; positive decreasing projected convergence with \(\kappa(0)>1\) is a sufficient physical condition. Indeed \(\Phi'=1-\bar\kappa\), where \(\bar\kappa=x^{-1}\int_0^x\kappa(z)\,dz\) strictly decreases from \(\kappa(0)>1\) to zero. A baseline is \(\kappa=2/(1+x)^2\), \(\Phi=x-2\log(1+x)\) up to a constant.

Let \(x_-(t),x_+(t)\) be the endpoints of the delay sublevel interval. Choose a nonzero real \(h\in C_c^\infty((\min\Phi,\Phi(0)))\) so both branches exist and are regular on its support. With \(x\mapsto x+\epsilon h(\Phi(x))\), define the new radial \(\Phi_\epsilon\) by composition with the inverse, as in (2). Both endpoints shift by exactly \(\epsilon h(t)\), so the area measure \(A_0(t)=x_+-x_-\) and aligned transfer are unchanged. The center, the Einstein ring and the entire exterior retain untouched neighborhoods. The new projected convergence stays strictly positive and decreasing: both margins on the compact region persist for a sufficiently small \(C^3\) perturbation. The above baseline is the projected form of a three-dimensional Plummer-type \(\rho_0(r)=C(1+r^2/2)^{-5/2}\). For completeness the circular family can keep \(\rho>0\) and decreasing too. In Abel deprojection,
\[
 \delta\rho(r)=-\frac1\pi\int_0^\infty f(\sqrt{r^2+z^2})\,dz,\qquad
 f(R)=\frac{\delta\kappa'(R)}{R},
\]
where the prime on \(\delta\kappa(R)\) here uses physical radius \(R\). The support of \(f\) is in a fixed \([R_0,R_1]\), \(R_0>0\). Thus \(|\delta\rho|\le R_1\|f\|_\infty/\pi\) and \(|\delta\rho'(r)|\le r R_1\|f'\|_\infty/(\pi R_0)\), with no perturbation beyond \(R_1\). On the compact range the baseline \(\rho_0\) and \(-\rho_0'(r)/r\) have positive lower bounds. These margins persist as \(\epsilon\to0\). This stronger spherical condition is asserted for the circular family, not for the nonsymmetric track example.

For angular offset \(y\) along a known signed axis, remove common factors from (1) to obtain
\[
 U(w,y)=\frac{i}{w}e^{-iwy^2/2}F(w,y)=\int_0^\infty e^{iw\Phi(x)}J_0(wy\sqrt{2x})\,dx.
\]
For \(\Im w>0\) this is convergent and differentiable, since \(\Phi(x)/x\to1\) at large \(x\); the exponential from \(J_0\) grows only as \(e^{O(\sqrt x)}\) for bounded \(y\). The known Bessel expansion gives
\[
 U(w,0)=\int e^{iw\Phi(x)}dx,\qquad
 \partial_y^2U(w,0)=-w^2\int x\,e^{iw\Phi(x)}dx. \tag{4}
\]
These are the transforms of the pushforward measures for \(dx\) and \(x\,dx\); define their cumulative functions \(A_0,A_1\). Real coherent full-frequency data with the causal prescription give the same uniqueness as Laplace data. More explicitly a tempered Fourier boundary ambiguity supported only at \(w=0\) transforms to a polynomial distribution on the full delay line; none can have bounded-below delay support unless zero. Thus the missing value at zero frequency does not introduce a new physical measure freedom. For a nonempty sublevel interval,
\[
 A_0=x_+-x_-,\quad A_1=(x_+^2-x_-^2)/2,\qquad
 x_\pm=A_1/A_0\ \pm\ A_0/2. \tag{5}
\]
The endpoints determine all sublevel sets and hence \(\Phi\). Its additive gauge does not affect \(\kappa\).

### Smooth relative gain and time origin

Suppose instead the measured germ is \(c(y)e^{iw b(y)}U(w,y)\), with \(c(y)>0\) and \(b(y)\) real and both twice differentiable near zero; the same unknown factors apply across frequencies. The frequency scale and geometric offset coordinate are still fixed. This only allows scalar gain/time ambiguity, not an arbitrary frequency-dependent unlensed waveform. To compare two lenses, normalize \(c(0)\) using \(\lim_{p\downarrow0}pU(ip,0)=1\) (from \(\Phi(x)/x\to1\)), and align the additive delay gauge using the common lower support of the aligned measures. Then \(\mu_0\) is common. Equality of first measured derivatives fixes the same \(c'(0)/c(0)\) and \(b'(0)\), since \(U_y(w,0)=0\). Equality of second derivatives leaves
\[
 \Delta L_1(w)=\frac{A+iw B}{w^2}U(w,0),\qquad
 L_1(w)=\int x\,e^{iw\Phi(x)}dx
\]
for real constants \(A,B\) from second-order gain and time shifts. On \(w=ip\), the right side has the form \(U(ip,0)(-A/p^2+B/p)\). Its inverse Laplace measure therefore has density \(-A\int_{-\infty}^t A_0(s)\,ds+B A_0(t)\): dividing a causal transform by \(p\) integrates once. But for all sufficiently high \(t\), both radial sublevels start at \(0\), and \(A_1(t)=A_0(t)^2/2\) in both lenses. The measure difference there is zero. Since \(A_0(t)\sim t\), the density terms grow as \(t^2/2\) and \(t\), so \(A=B=0\). Equation (5) then applies as before. This does not treat nonsmooth per-offset phase normalizations (e.g. tracking a minimum through a caustic without unwarping its phase convention).

### Quantified separation for one displaced circular source

In the circular shear, assume \(h\ge0\) with support in a delay interval of length \(T\), both branch endpoints and their shifted positions \(\le X\) there. Write \(M=\int h(t) A_0(t)\,dt>0\) and \(a=wy\). The weighted sublevel integral for \(U(w,y)\) changes by
\[
 \Delta A_f(t)=\int_0^{\epsilon h(t)}[f(x_+(t)+s)-f(x_-(t)+s)]\,ds,\quad f(x)=J_0(a\sqrt{2x}).
\]
The elementary angular-average formula for \(J_0\) and \(|\sin v-v|\le |v|^3/6\) give \(|f'(x)+a^2/2|\le a^4x/8\) for real \(a\); use \(\langle\cos^2\theta\rangle=1/2\), \(\langle\cos^4\theta\rangle=3/8\). Thus \(\Delta A_f=-\epsilon h A_0 a^2/2\) up to a pointwise error \(\le|\epsilon|h A_0 a^4 X/8\). Integrating the compact Stieltjes difference by parts yields
\[
 |\Delta U(w,y)|
 \ \ge\ \frac{|\epsilon||w|^3 y^2}{2}\,M
  \left[\cos(|w| T/2)-\frac{w^2y^2 X}{4}\right] \tag{6}
\]
for \(0<|w|T<\pi\); to lower bound the leading transform, rotate by the phase at the midpoint of the support. Thus any known sufficiently small but nonzero offset distinguishes the constructed family at a real frequency satisfying the displayed bound. This absolute distinction is conditional on scale/gauge, and the small offset produces a quadratic signal; it is not an assertion of detectability from finite noisy data. At the level of full spectra, these compact changes also cannot be made to agree by a constant gain and time shift at that fixed offset once the calibrated spectra differ: outside the affected delay range their cumulative sublevel areas are identical, strictly increasing for large \(t\) and asymptotic to the same leading \(t\) after angular normalization. Matching by a gain/shift would force gain one by the leading asymptotic and shift zero by strict increase.

## 3. Explicit noncircular track witness with certified positive convergence

Choose a baseline circular \(\kappa_0(x)=3/2-x\) for \(x=r^2/2\le7/5\); extend smoothly and strictly decreasing with a positive finite-mass tail. One such extension is obtained by letting \(k(s)\) be a nonnegative smooth transition from \(1/(3/2-s)\) near \(7/5\) to 1 before \(29/20\), and taking \(\kappa_0(x)=(1/10)\exp[-\int_{7/5}^x k(s)\,ds]\) thereafter. A smooth convex combination in that transition interval is positive, and the piece before the cut matches all derivatives by using \(k(s)=1/(3/2-s)\) initially. In the inner region \(S=\Phi=(x-1)^2/4\) after an additive shift.

Use \(B(z)=256z^4(1-z)^4\) on \([0,1]\), zero outside, \(C^3\). Let
\[
 H(u,t)=B\!\left(\frac{u-3/2}{1/16}\right)
        B\!\left(\frac{t-1/64}{5/256}\right),\qquad
 \epsilon=1/200000.
\]
On its preimage \(3/2\le u\le25/16\), \(5/4\le x\le11/8\); on these fibers \(x\ge u^2/2>1\), the delay rises strictly with \(|v|\), and at \(v=0\) one has \(u^2/2<5/4\). Thus all modified fiber portions avoid every fiber critical point, consistent with the no-new-images hypothesis. Baseline convergence is at least \(1/8\) there, \(\|DS\|\le3/8\), \(\|D^2S\|\le3\). The replay verifies \(\|B'\|_\infty\le4,\|B''\|_\infty\le32\) by checking elementary extrema in \(q=2z-1\). Consequently for \(K(u,v)=H(u,S(u,v))\), the chain rule gives \(\|DK\|\le141\), \(\|D^2K\|\le31000\), with operator norms. Put \(\delta=141|\epsilon|\), \(L=(1-\delta)^{-1}\). At transformed points the inverse-map chain rule bounds the Hessian change by
\[
 \|D^2S_\epsilon(\chi_\epsilon(\boldsymbol x))-D^2S(\boldsymbol x)\|
 \le 3 (L+1)\delta L+\frac38 L^3\,31000|\epsilon|<1/8.
\]
Here \(\|D^2\chi^{-1}\|\le L^3\|D^2\chi\|\); half the trace difference is bounded by the same operator-norm difference. Hence \(\kappa_\epsilon>0\) everywhere.

At the midpoint of both bumps with \(v>0\), \(u=49/32\), \(x=1+\sqrt{13/128}\), and \(H=1\), \(DK=0\). Exactly there \(\Delta K=-32(16^2+(256/5)^2\|DS\|^2)\) and the transformed convergence is \(\kappa_\epsilon(\chi(q))=\kappa_0(q)+\epsilon S_v(q)\Delta K/2\). Correcting the comparison to the exact same point as in the *unperturbed* model adds \(\epsilon v+\epsilon^2/2\), since \(\kappa_0\) is linear in \(x\). The exact replay bounds the resulting density difference away from zero by \(1/400\), with the unperturbed same-point convergence below \(19/100\). The whole-line wave and geometric degeneracies follow exactly from (2)-(3), not from small estimated numerical differences. The explicit polynomial witness is \(C^3\); the general theorem has \(C^\infty\) witnesses, and this positive-margin witness can also be smoothed with arbitrarily small \(C^2\) changes to \(H\).

Run the offline verify.py for exact rational geometric, positivity, and circular sublevel moment checks. No numerical Fresnel quadrature or data fit is used to certify equality.
