# How much atmospheric information can arrive before opaque contact?

## Scope and decision

A heterogeneous but static stellar photosphere cannot cause an **occultation** when an exactly located opaque planet is entirely outside its emitting disk. A planetary absorbing annulus can. This suggests a geometric check on transit-light-source contamination, distinct from fitting the overall ingress/egress durations. Its usefulness is limited by how little of the annulus is backlit and how briefly. We give the exact upper-envelope annular profile and a photon-information ceiling, then evaluate the privately frozen 100 ppm benchmark in [FROZEN.md](FROZEN.md). No science light-curve outcome is fit.

**Result for the frozen benchmark.** Under the calibrated conditions below, ten ideal complete transits cannot distinguish the stipulated e-like atmospheric feature using only the two exterior-contact windows with both size 5% and power 80%. Certified power is less than 59.1%, and at least 18 ideal visits would be necessary, not sufficient. This is a result about one geometric discriminator and a defined absorber class; it does not apply to full transit fits or to all atmospheric scenarios.

## 1. Geometric model

Distances are in stellar radii. The stellar image is the closed unit disk. At each wavelength (a subscript \(\lambda\) is suppressed in this section), its nonnegative intensity may otherwise vary arbitrarily, but satisfies \(I(x)\le L\bar I\), where \(\bar I\) is the actual disk mean. The core has independently known opaque radius \(p>0\), and the outside-core radial attenuation \(a(r)\) is nonnegative and at most one, vanishes outside \(p<r\le p+H\), is centered on the same core, and acts purely by removing light. Its area as a fraction of a uniform stellar disk is
\[
 \Delta=\int_p^{p+H}2r a(r)\,dr .
\]
This is **not** the observed interior depth in a nonuniform stellar image. A nonnegative difference in opacity between independently fixed absorption and reference models has the same formula. If a reference model contributes tails too, the proof compares suitable known-base rate models and the rate floor must be propagated; the numerical example is opaque-only versus extra attenuation.

Write \(d\) for star-to-planet projected separation and \(s=d-1-p\). The pre/post opaque-contact windows have \(s\ge0\); the compact model gives zero signal once \(s\ge H\). For \(d>1\) and the radial ranges of interest, the fraction of the circumference at radius \(r\) around the planet that lies in the stellar disk is
\[
 f(d,r)=
 \begin{cases}
 \pi^{-1}\arccos\!\left(\dfrac{d^2+r^2-1}{2dr}\right),&d-r<1<d+r,\ r<d+1,\\
 0,&\text{otherwise}.
 \end{cases}                                                     \tag{1}
\]
(The condition \(r<d+1\) explicitly covers rings beyond the far side of the star.) For all \(r\), the fractional atmospheric signal \(m(s)\) relative to the true unocculted flux satisfies
\[
 0\le m(s)\le L \int_p^{p+H} f(1+p+s,r)\,2r a(r)\,dr.              \tag{2}
\]
This is equality for \(L=1\) with a uniform image. The **absence** of an opaque-only decrement at \(s\ge0\) does not use the brightness cap; the bound on the alternative does use it. Stellar total-flux changes or instrumental drifts can still mimic a time-domain decrement and are not ruled out by contact geometry.

### Exact profile extremizer

Assume \(p+H<\sqrt{(1+p)^2-1}\) as holds in the benchmark. For every fixed exterior separation \(d\), the function \(f(d,r)\) is nondecreasing in \(r\) over the allowed shell. To see this where it is nonzero, differentiate its arccos argument, obtaining \((r^2-d^2+1)/(2dr^2)<0\); the fraction is zero until the ring first touches. If
\[
 r_0=\sqrt{(p+H)^2-\Delta},\qquad a_*(r)=\boldsymbol{1}_{[r_0,p+H]}(r),
\]
then \(a_*\) simultaneously maximizes the RHS of (2) at **every** exterior separation. Its total weighted area equals \(\Delta\). Subtract the value \(f(d,r_0)\) inside the integral of \(a_*-a\); the two factors now have the same sign on both sides of \(r_0\), proving the claim. Thus
\[
 m(s)\le L\,[\mathcal A(1,p+H,d)-\mathcal A(1,r_0,d)]/\pi,         \tag{3}
\]
where \(\mathcal A\) is the exact common area of two circular disks. No small-planet, flat-limb or delta-function reduction is used in (3). This profile is optimal for the geometric upper-envelope problem, not a claim that the most likely planetary opacity is concentrated at the top. It attains the full pointwise envelope for a uniform star \(L=1\). For \(L>1\) a real fixed brightness map need not saturate it; the bound remains valid. For models restricting opacity to be monotone decreasing (rather than permitting a chromatic excess shell), the supremum cannot be higher.

## 2. Explicit photon-information ceiling

The outside-contact geometry gives a simpler bound valid without the additional thinness condition used in (3), and retains the star's finite limb curvature. Write \(z=r-p\) and suppose the ring is partially covered. If \(\theta=\pi f(d,r)\) is its covered half-angle, then
\[
 \tan^2(\theta/2)=\frac{1-(d-r)^2}{(d+r)^2-1}.
\]
Writing \(q=z-s\), the numerator is \(q(2-q)\le 2q\), and the denominator is at least \(4p(1+p)\), because \(d\ge1+p\) and \(r\ge p\). Applying \(\arctan v\le v\) gives
\[
 f(1+p+s,p+z)\le \sqrt{\frac{2(z-s)_+}{\pi^2p(1+p)}} .
\]
The same bound applies in the uncovered case. Therefore
\[
 m(s)^2\le\frac{2L^2\Delta^2}{\pi^2p(1+p)}(H-s)_+,\qquad
 m(s)\le L\Delta .                                               \tag{4}
\]
An optional contact timing guard \(s\ge s_0\) makes the integrated information factor \(H^2\) below become \((H-s_0)_+^2\). For a single wavelength, Minkowski's integral inequality also permits replacing \(\Delta H\) by \(\int_p^\infty (r-p)2r a(r)\,dr\), without compact support, provided the same radial speed and photon rate are bounded over the selected exterior windows. Indeed the geometric kernel \(\sqrt{(z-s)_+}\) has \(L^2(ds)\) norm \(z/\sqrt{2}\). This is only an optional first-moment generalization; the data benchmark uses the compact support hypothesis as frozen.

Consider the optimistic experiment where every source photon in the band has an exact detection time and the assumed wavelength mark, the unocculted photon rate \(\phi_\lambda(t)\) is known, the core contacts (which must be independently constrained) and alternative template may be given to the test, and visits are independent conditional on the fixed models. Under \(P_0\) the rate on these windows is \(\phi_\lambda(t)\); under \(P_1\) it is \(\phi_\lambda(t)(1-m_\lambda(t))\). The weights in the quadratic cap are spectral *photon-rate weights*, not unweighted wavelength averages. At every time assume
\[
 \int \phi_\lambda(t)\,d\lambda\le\Phi,\qquad
 \int \phi_\lambda(t)\Delta_\lambda^2 d\lambda\le\Phi\,\Delta_{\rm rms}^2 ,
 \qquad 0\le m_\lambda(t)\le\epsilon<1 .
\]
A uniform per-channel cap \(\Delta_\lambda\le\Delta_{\rm rms}\) suffices. Ordinary coarse-bin eligibility is not enough to bound an ideal experiment retaining larger unresolved line-core information unless the quadratic condition actually holds at that resolution.

The Poisson-process likelihood gives the relative entropy, in natural-log units,
\[
 D(P_1\Vert P_0)=
 \int dt \int d\lambda\,\phi_\lambda(t)\big[(1-m_\lambda(t))
                         \log(1-m_\lambda(t))+m_\lambda(t)\big].
\]
The bracket is at most \(m_\lambda^2/[2(1-\epsilon)]\), since it vanishes to second order at zero and its second derivative is \(1/(1-m_\lambda)\). If the monotone radial speed \(|ds/dt|\) is bounded below by \(u>0\) on each of the two windows in each of \(n\) visits, then integrating (4) and adding visits gives
\[
 \boxed{\quad
 D(P_1^{(n)}\Vert P_0^{(n)})
 \le \frac{n\Phi L^2\Delta_{\rm rms}^2 H^2}
            {\pi^2 p(1+p) u (1-\epsilon)} .
 \quad}                                                         \tag{5}
\]
Speed and time here use stellar radii and seconds. The bound includes *both* exterior windows. Additional read/background noise, less than unit throughput, integration over finite exposures, missed windows, and data processing do not increase the information under the stated hypothesis-independent observation/noise model. Absolute stellar calibration must give a flux cap at transit times; absolute calibration does not make the occulted annular response independent of the surface map.

Formula (5) is not exact as a finite information maximization for an arbitrary stellar map. Its leading geometric coefficient is nonetheless sharp in a thin-shell limit. For a uniform star, constant photon rate and radial speed, and a family with \(\Delta(H)>0\) realized by outer opaque shells of thickness \(o(H)\), take \(H/p\to0\) at fixed \(p\). The exact ring fraction (1), at \(s/H\in[0,1)\) and \(r=p+H+o(H)\), divided by \(\sqrt H\), tends to \(\sqrt{2(1-s/H)/[(1+p)p]}/\pi\). The curved-limb bound (4) dominates the square over the full scaled interval. With \(\epsilon\to0\) the true two-window KL thus equals
\[
 \frac{\Phi\Delta^2 H^2}{\pi^2 (1+p) p u}(1+o(1))
\]
per visit. The factor \(1/(1+p)\) is the curved stellar limb correction in this control. An optional substitute for the right hand side of (5) is direct integral of (3) and the exact Poisson bracket; it is no larger when (3) is evaluated using the exact flux, geometry and speed.

### Testing gate

Let the reported decision be binary, including an auxiliary independent RNG if desired. By data processing,
\[
 D(P_1^{(n)}\Vert P_0^{(n)})\ge
  d_{\mathrm{Bern}}\!\left(P_1(\text{reject})\Vert P_0(\text{reject})\right).
\]
At size at most 0.05, power at least 0.80 requires KL at least \(d_{\mathrm{Bern}}(0.80\Vert0.05)>1.9064\); if power is at least 0.591 the requirement exceeds 1.1149. These are limitations for any decision rule on the specified experiment, not asymptotic SNR heuristics. No agglomeration across unmodeled shared calibration sessions or post-hoc windows can be counted as extra independent information.

## 3. Independent benchmark inputs and checks

- Geometry scale: [Agol et al. 2021, Tables 2, 5–7](https://doi.org/10.3847/PSJ/abd022) report \(R_p/R_\star=0.07079\pm0.00055\), \(a=0.02925^{+0.00025}_{-0.00025}\) au, \(P\simeq6.1010\) d, stellar radius \(0.1192\pm0.0013\,R_\odot\), and ingress scale 3.825 min under their photodynamic transit model. The implied circular central speed is about \(6.3\times10^{-4}\,R_\star/{\rm s}\); these values motivate the looser frozen \(p\ge0.068\) and \(u\ge4.5\times10^{-4}\,R_\star/{\rm s}\) regime, not an independently exact opaque-contact timestamp for any JWST observation. \(H=0.01R_\star\) is about 830 km at that scale and the geometric exterior window then lasts at most 22.3 s on each side under the speed cap. \(L\le2\) is a declared independent brightness cap, not established by a disk-integrated spectrum. A 100 ppm uniform-illumination area feature corresponds to about 59 km of thin opaque shell at \(p\approx0.071\); \(H\) is a generous permitted support relative to that *equivalent opaque thickness*, and it has not been measured.
- Stellar photon rate: The [IRSA AllWISE native source table and response](SOURCES.md) give \(W2=9.807\pm0.020\) for the correct source. This is independent of the target transit spectra. Applying the conservative response-integral and brightness conditions in [SOURCES.md](SOURCES.md), and using 26 m² (above JWST's 25.4 m² clear collecting area), gives a **quiet-source** upper value below \(4.302\times10^6\) photons/s over 4.10–4.50 μm at unit throughput. The proof explicitly allows double this upper value as a conditional cross-epoch safety margin; WISE does not certify a bound on all JWST epochs or flares. No intrinsic stellar brightness evolution is allowed *as an unknown mimic*; the measured independent unocculted light baseline is an input to the ideal experiment.
- Gate: For \(\Phi<8.604\times10^6\ {\rm s}^{-1}\), \(p\ge0.068\), \(p\le1\), \(H\le0.01\), \(L\le2\), \(\Delta_{\rm rms}\le10^{-4}\), \(u\ge0.00045\ {\rm s}^{-1}\), \(m_\lambda\le L\Delta_\lambda\le L(2pH+H^2)\le0.0402<0.041\) even if individual channels exceed the RMS. Using \(\epsilon=0.041\), exact interval arithmetic gives \(D(P_1^{(10)}\Vert P_0^{(10)})<1.113\); the \(n=17\) ceiling is below 1.892. Thus the 80% gate fails for 17 or fewer perfect transits, and even 59.1% power is unavailable with ten at the prescribed size. These are lower requirements, not guaranteed achievements with more visits.

## 4. Evidence boundary

Morris and colleagues already used the time duration of ingress and egress to help recover radii around spotted stars. The result here is only the optimized true annular overlap and model-level ceiling for the **spatial-map-independent null zero** in the exterior wavelengths/times. It does not constrain how much usable information survives partial or full ingress when there are independent transit chords, calibrated local spectra or good stellar models, nor does it supersede a published multi-transit retrieval. It cannot decide the presence or absence of any planet's atmosphere. Nonaxisymmetric extended tails, atmospheric brightening/refraction, ring-like occultation, high narrow-line equivalent depth, time-dependent stellar activity, and a wavelength-dependent stellar radius or unknown core ephemeris require different calculations. A geometric extra-contact decrement of itself would not identify a specific absorber.

The independent calibration audit and exact numerical bounds are replayable. No native transit light curves were used to test an atmospheric hypothesis. The recommendation for the frozen comparison is **no outcome fit on the two clean exterior windows alone**; park this route as a structural/power limit, not an observational discovery or a new description of stellar contamination.
