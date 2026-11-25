# Independent audit: lens-delay rearrangements and offset information

**Verified in the ordinary single-plane scalar thin-lens model:** (i) perfectly aligned broadband complex data do not identify a circular finite-mass lens, even among projections of smooth nonnegative decreasing spherical three-dimensional densities; (ii) the second offset coefficient identifies the circular potential up to an additive delay when all delay sublevels are single intervals in area coordinates; (iii) every fixed affine source-offset line admits nonidentifying compact, physically positive lens perturbations, with all geometrical image positions and jets identical on the track; (iv) a calibrated open two-dimensional offset patch at just one nonzero real frequency distinguishes compact potential differences. Single-frequency absolute data in (iv) cannot be replaced by arbitrary per-offset complex gauges. An additional first-transverse-derivative sufficiency rule for interval fibres is given below.

## 1. Circular reduction, physically admissible shear, and two inverse branches

Use dimensionless lens-plane coordinates and fix the source and frequency scales. Set $x=r^2/2$, so $r\,dr=dx$, convergence $\kappa=\Sigma/\Sigma_{\rm crit}$, $m(r)=2\int_0^r\kappa(q)q\,dq$, and $\psi_r=m(r)/r$. Thus in area coordinates
\[
 \tau(x)=x-\psi(\sqrt{2x}),\qquad
 \kappa(x)=1-\tau_x(x)-x\tau_{xx}(x).
\]
The derivative identity follows by differentiating $x(1-\tau_x)=\int_0^x\kappa(q)dq$. Start with a Plummer projected profile $\kappa(x)=k(1+x)^{-2}$, $k>1$, so $\tau(x)=x-k\log(1+x)$ (a time constant is fixed here). The unique minimum occurs at $x_E=k-1$, the central delay is $\tau(0)=0$, and strict monotonicity holds on each branch. For delays $s$ strictly between that minimum and 0 let $x_-(s),x_+(s)$ be the decreasing and increasing inverse branches. Replace both by $x_\pm(s)+h(s)$ with smooth $h$ of compact support in that open delay interval. A sufficiently small $C^4$ perturbation preserves positivity, ordering and branch monotonicity, because the inverse-branch derivatives have strict sign and the support stays away from both endpoints. Inverting produces a smooth $\tau_h$ near $\tau$ in $C^4$, exactly unchanged near the center, Einstein ring and outside a compact annulus or union of annuli.

The reference $\kappa,\ -\kappa',\ \kappa''$ have positive minima on the compact altered region. The formula for $\kappa$ then guarantees those strict signs survive for small $C^4$ shifts. At infinity the potential and enclosed mass are exactly unchanged. For the spherical deprojection write $\Sigma_v(v)=\Sigma(\sqrt v)$, $v$ projected squared radius. The Abel inverse is
\[
 R(v)=\rho(\sqrt v)= -\pi^{-1}\int_0^\infty \Sigma_v'(v+s)\,s^{-1/2}\,ds .
\]
Indeed the forward formula is $\Sigma_v(v)=\int_0^\infty R(v+s)s^{-1/2}ds$ and the two kernels convolve to $\pi$ (differentiate then integrate using exterior decay). The same inversion differentiated once gives $R'\le0$ from $\Sigma_v''\ge0$. Reference decay and compact perturbations justify both operations. Hence $\rho$ is nonnegative and strictly decreasing in radius away from the center; its total mass is finite by positive projection/Fubini. This argument uses convexity of projected $\Sigma$ in squared radius, not simply monotonicity of $\Sigma$.

In these coordinates the aligned complex diffraction amplification is $(w/i)\int_0^\infty e^{iw\tau(x)}dx$. Its datum is the pushforward of $dx$ by $\tau$, with cumulative sublevel width $W(s)=x_+(s)-x_-(s)$ (put $x_-=0$ after the central delay). The common shear preserves $W$ for every delay, including the entire large-delay tail; hence equality is exact, at absolute complex amplification and at all frequencies, not just leading stationary phase or magnification. Distinct shifts give distinct potentials; the freedom contains an infinite-dimensional small ball of shifts supported inside any fixed compact subinterval of the two-branch delay range.

## 2. Abel boundary and offset-squared inversion

The exact angular integral for radial offset $y$ is
\[
 F(w,y)=\frac{w}{i}e^{iwy^2/2}
         \int_0^\infty e^{iw\tau(x)}J_0(wy\sqrt{2x})\,dx.
\]
For $\Im w>0$ and $\tau(x)/x\to1$ at infinity, this integral and its offset derivatives converge absolutely, since the possible Bessel growth is exponential in $\sqrt x$, not $x$. The physical oscillatory values are the Abel limit from $\Im w>0$. For the aligned integrals the corresponding delay measures are supported on a half-line and have polynomial cumulative growth ($W=O(s)$ and the area moment $O(s^2)$). Their Fourier transforms can be taken as tempered distributions. Using $\int e^{iws}$ as transform convention, let $\widehat A,\widehat B$ be the transforms of $\tau_\#dx$ and $\tau_\#(x\,dx)$. Then
\[
 F(w,0)=\tfrac wi\widehat A(w),\qquad
 [y^2]\left(e^{-iwy^2/2}F(w,y)\right)=-\tfrac{w^2}{2}\tfrac wi\widehat B(w).
\]
The factor $1/2$ comes from $J_0(z)=1-z^2/4+\cdots$. Do not confuse the $y^2$ coefficient with the second derivative in $y$ (twice the coefficient). Knowledge away from DC at both frequency signs suffices in this distributional convention: any missing distribution supported at $w=0$ in a transform would yield a polynomial in delay, which cannot vanish on an entire interval below the delay support except when zero. For real potentials, negative frequencies follow by complex conjugation.

Suppose the area-coordinate sublevel sets are intervals (in particular a circular single well, possibly with one-branch levels). The two cumulative observables are $W$ and $M=\int_{\tau\le s}x\,dx=(x_+^2-x_-^2)/2$. For all $W(s)>0$ one obtains
\[
 x_\pm(s)=M(s)/W(s)\ \pm W(s)/2.
\]
These determine the delay function and $\kappa$. Where the sublevel touches the origin the formula recovers $x_-=0$. The minimum endpoints are determined by limits and the class's continuity, not by a division by zero. For the sheared family $M_h-M=hW$, so every nonzero shear is detected by the second coefficient. If a lens admits multi-interval sublevels, two moments need not determine all endpoints; no such generalization is claimed.

A common absolute time shift translates $A,B$ together and does not change the reconstructed $\kappa$. As to per-offset gauges, if one compares two physical models only up to a frequency-independent multiplicative amplitude $a(y)$ and a delay $c(y)$ (factor $e^{iw c(y)}$), the raw second coefficient gains $(a_2/a_0+iw c_2)F_0$ after the common gauge at zero offset is removed. This algebra shows the gauge must be handled, not ignored. In the exact full-band model, the prescribed low-frequency normalization removes both of these relative gauges: with $w=i\sigma$,
\[
 F_0(i\sigma)\to1,\qquad
 \frac{[y^2]F(i\sigma,y)}{F_0(i\sigma)}
 =-\frac{\sigma}{2}+\frac{\sigma^2}{2}
 \frac{\int x e^{-\sigma\tau(x)}dx}{\int e^{-\sigma\tau(x)}dx}
 =o(\sigma).
\]
Rescale $\sigma x$, use $\tau(x)/x\to1$ and an exponential dominating function, then both moments of the unit exponential are 1. The full small-$y$ ratio also satisfies $F(i\sigma,y)/F(i\sigma,0)=1+o(\sigma)$ at fixed $y$, by the analogous $I_0$ expansion. These limits can be extracted from exact Abel data (not from intensity or from noisy finite-band samples). For the finite-center isodelay class, an independent cross-check is that equal aligned $W$ forces the same outer inverse beyond both central delays, so the first-moment-density difference is compact and has zero integral; the difference of raw offset coefficients is then $O(w^4)$ at zero and cannot be disguised by a constant plus linear-in-frequency gauge. Uncalibrated frequency-dependent per-offset source waveforms or frequency/angle scale changes are a different problem.

## 3. Entire straight-track data can agree for a nonaxisymmetric positive lens family

The two-dimensional dimensionless arrival delay is
\[
 T_\psi(u,v;y)=\tfrac12|(u,v)|^2-\psi(u,v)-y_u u-y_v v+\tfrac12|y|^2.
\]
Rotate the track to $y=(t,c)$; **include the normal offset** in $\sigma_0(u,v)=\tfrac12(u^2+v^2)-\psi_0(u,v)-cv$. Choose a smooth $H(u,s)$ nonzero on realized levels and supported in a compact rectangle whose entire preimage under $(u,\sigma_0)$ is away from $\partial_v\sigma_0=0$ and lies in a region where $\kappa_0>0$. Such rectangles exist for the logarithmic Plummer potential and every $c$: with $u$ in a bounded interval, large delays force large $|v|$ and all $\partial_v\sigma_0$ there are bounded away from zero. Define
\[
 \Phi(u,v)=(u,v+\epsilon H(u,\sigma_0(u,v))),\quad
 \sigma_\epsilon=\sigma_0\circ\Phi^{-1},\quad
 \psi_\epsilon=\tfrac12(u^2+v^2)-cv-\sigma_\epsilon.
\]
Because $\sigma_0$ is proper, the preimage support is compact. The warp is identity outside that compact, and $\partial_v\Phi_v=1+\epsilon H_s\sigma_{0,v}>0$ for small $|\epsilon|$; hence $\Phi$ is a diffeomorphism. There is a compact margin around all moved points where the reference convergence has positive minimum. Therefore the $C^2$-small compact $\delta\psi$ preserves positive convergence $\kappa_\epsilon=\Delta\psi_\epsilon/2$. Its total projected mass is the same since $\int\Delta\delta\psi=0$. The compact perturbation also equals the physical logarithmic potential of its signed convergence difference: the log convolution tends to zero at infinity due to zero total perturbation mass; their difference is harmonic globally and tends to zero. No extraneous harmonic or external shear has been added.

For every compact test function $f(u,s)$ the Jacobian change of the joint $(u,\sigma)$ pushforward is
\[
 \epsilon\int du\,dv\, f(u,\sigma_0)H_s(u,\sigma_0)\partial_v\sigma_0=0.
\]
Indeed for each fixed $u$ the inner integral is the difference of a primitive in $s=\sigma_0$ at both ends of the fibre, which both tend to $+\infty$. The equal joint pushforwards imply equal diffraction for every $t$ and every $w$ by the same phase $w(\sigma-tu)$; Abel convergence, or compact subtraction, justifies the noncompact test. This proof never assumed the track goes through the lens center. The result can be nontrivial and nonaxisymmetric because generic $H$ moves two differently placed fibre points through a non-axis-preserving map, and where $\sigma_{0,v}\ne0$ a sufficiently small nonzero move changes $\sigma_\epsilon$. It is not only a reflection ambiguity.

All stationary points at each track offset correspond exactly because $\sigma_\epsilon-tu=(\sigma_0-tu)\circ\Phi^{-1}$. They have $\sigma_{0,v}=0$, where the map is identity in a neighborhood, so positions, delays, and all spatial derivative jets are unchanged, and no new stationary points arise. In fact the deformation can leave stationary data identical in a narrow source strip: on the moved region both $|\partial_v\sigma_0|$ and $|\partial_v\sigma_\epsilon|$ have a common positive lower bound, so small transverse track shifts still place no stationary points there. Wave-optical equality itself was proved only on the track.

## 4. How more offset information can distinguish the lenses

Suppose two smooth delays differ by a compact perturbation. At fixed real $w\ne0$, after removing the common geometric source phase, their amplification difference is the two-dimensional Fourier transform at $w y$ of the compact function $e^{iw\tau_1}-e^{iw\tau_0}$ with $\tau_i=\tfrac12|\theta|^2-\psi_i$. If calibrated absolute complex values agree on an open offset patch, this Fourier transform vanishes on an open real set and hence identically, by real analyticity. Fourier uniqueness and continuity give $\tau_1-\tau_0\in (2\pi/w)\mathbb Z$ pointwise; on the connected plane it is constant, and zero on the exterior. This is exact injectivity by analytic continuation, not a stability claim. Per-offset unconstrained complex amplitude/time calibration would make one-frequency equality meaningless. Also a common global time convention must be fixed for the condition that the potential difference vanishes outside a compact set.

For an additional sufficient measurement more closely tied to the track, suppose for each $u$ the sublevel sets $\{v:\sigma(u,v)\le s\}$ are intervals. Along $y=(t,c)$, remove $e^{iw|y|^2/2}$ before differentiating. The track itself gives the joint transform of the $(u,\sigma)$ pushforward of $du\,dv$, evaluated at $(\xi,w)=(-wt,w)$; the first **transverse** derivative gives $-iw$ times the same transform for the pushforward of $v\,du\,dv$. For two potentials differing compactly those differences of signed measures are compact in $(u,s)$. Thus values for $t$ and $w$ in open intervals (with $w\ne0$) determine equality of the measures by Fourier analyticity, because $(t,w)\mapsto(-wt,w)$ maps open sets to open sets. Inverting and taking cumulative measures in delay gives fibre widths $W(u,s)$ and moments $M(u,s)$ for almost every $u$; regularity extends to ordinary nonsingular levels. The same interval formula $M/W\pm W/2$ recovers both transverse endpoints. This extra sufficiency criterion needs the stated interval-fibre hypothesis and absolute transverse derivative, unlike the open-two-dimensional-patch result.

## Primary comparators and limits of novelty

* [Ulmer–Goodman (1995)](https://doi.org/10.1086/175422) introduced the time-delay contour diffraction calculation. [Diego et al., §4 and Eq. 13](https://doi.org/10.1051/0004-6361/201935490) explicitly Fourier-transform equal-delay-bin areas, and [Tambalo et al., §IV.A, Eqs. 31–32](https://arxiv.org/abs/2210.05658) give the pushforward and contour integral in the modern waveform notation. The forward delay-area identity is not new.
* [Asorey–Cariñena–Marmo–Perelomov, §2](https://arxiv.org/abs/0707.4465) explicitly give the classical Abel isoperiodic inverse-branch shear, including $x\mapsto x+g(U(x))$. Neither the abstract redistribution argument nor the optical Fourier transform is claimed new. The checked contribution is their realization as exact thin-lens degeneracies under positivity/finite mass (and radial spherical deprojection), together with the specific offset-moment diagnostics.
* [Meena 2024, §§4.1 and 5](https://arxiv.org/html/2409.00145v2#S5) distinguishes a point/sphere monopole pair that **changes delay-bin area**; §5 argues broad local degeneracies can break. The present family shows that changing a Fermat surface away from the geometrical images does not by itself force a change in diffraction: the area pushforward can stay exact, even on an entire offset track. This is a limitation to the broad reading, not a rebuttal to that paper's example. Its §4.2 and [Chen et al.](https://arxiv.org/abs/2406.06463) concern global mass-sheet and scaling/source-distance ambiguities, which differ from these equal-calibration compact perturbations.
* [Maretzke's inverse Fresnel work](https://arxiv.org/abs/1805.06185) treats finite-detector inverse Fresnel regimes; the open-patch analytic argument above is a straightforward compact-support Fourier uniqueness consequence. The track is a collection of ideally calibrated stationary measurements, not a claim that one moving astrophysical source supplies those data. Source and coordinate scaling, noisy/stochastically varying data, multiple curved/nonparallel tracks, and a general nonaxisymmetric finite-data stability theory are out of scope. Targeted primary search found no matching physical construction with the same line-scan and local-jet scope; no exhaustive priority claim.
