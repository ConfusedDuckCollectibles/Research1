# A two-grating separator for passive positive-path phonon transport

## Frozen data contract and result

Consider the real, normalized energy modulation after a small, approximately impulsive spatial sinusoidal temperature perturbation in a homogeneous bulk sample. Compare two magnitudes of the *same directional* wavevector at the same true time after the impulse. Preparation initially weights phonon modes by their equilibrium heat capacities, and readout is the same energy-defined temperature, not a rate-weighted fitting temperature or an arbitrary optical mixture. Fix an independently calibrated heat-capacity-weighted marginal \(\pi\) of the projected group velocity \(v\), with \(|v|\le V\). In particular the directional speed bound and velocity law are not tuned from the two responses under test. Define at a chosen time \(t>0\)
\[
 c=q Vt,\qquad B(c)=\int\cos(c v/V)d\pi(v),\qquad D(c)=\Theta(q,t)-B(c).
\]
The ballistic transfer \(B\) is the collisionless thermal-preparation baseline, not an independently adjustable component for each grating.

**Null \(\mathcal K\).** Any energy transport described by a positive random passive-carrier trajectory with projected velocity process \(v_s\in[-V,V]\), with marginal \(\pi\) at every time and displacement \(X_t=\int_0^t v_s ds\). The real grating transfer is \(E\cos(qX_t)\). This includes stationary arbitrary positive quasiballistic velocity mixing, reversible elastic direction changes and energy-conserving mode-dependent single-mode relaxation with matched preparation/readout. It is broader than a fit to any prescribed relaxation times; no collision number, spectral gap or independence assumption is used. It excludes coherent mode contributions, signed sources/readouts, and collision operators that cannot be represented by positive paths of this type.

For any \(0<c_1<c_2<\pi/2\), the null obeys the **exact finite two-grating cone**
\[
 0\le D(c_1),\qquad
 \frac{\cos c_2}{\cos c_1}\frac{D(c_1)}{c_1^2}
 \ \le\ \frac{D(c_2)}{c_2^2}
 \ \le\ \frac{D(c_1)}{c_1^2}. \tag{1}
\]
Both coefficients in the ratio bounds are sharp as limits even for the fixed uniform projected velocity marginal, using nonzero reversible Markov collisions. A passive, explicitly energy- and momentum-conserving normal collision model with small positive resistive rate violates (1) with **both grating amplitudes positive and both ballistic-relative differences positive**. For the isotropic gray Debye–Callaway model below, at \(\Gamma_N t=8,\Gamma_R t=7/100\), \(c_1=3/5,c_2=6/5\), a rational independent-oracle replay certifies
\[
 0.001047<D(c_1)<0.001048,\quad
 0.000478<D(c_2)<0.000479,\quad
 4\frac{\cos(6/5)}{\cos(3/5)}D(c_1)-D(c_2)>0.0013 . \tag{2}
\]
The witness is a full bounded continuous-angle transport evolution, not a finite-speed discretization or two fit exponentials. This rules out the null for the synthetic pair even with, e.g., \(4\times10^{-4}\) absolute error per ballistic-relative difference (exact dimensionless wavevector ratio and conservative verified cosine coefficients). It does not assert achievable precision in a current experiment. A violation rules out this null, not uniquely selects collective phonon momentum: optical or preparation effects, boundary leakage, velocity bound errors, and coherent phonon mechanisms must be independently excluded.

**Prior anchors.** The heat-capacity-weighted ballistic grating function and relaxation-time controls are standard [Hua–Minnich, 2014](https://doi.org/10.1103/PhysRevB.89.094302), [Forghani–Hadjiconstantinou–Péraud, 2016, Appendix D](https://doi.org/10.1103/PhysRevB.94.155439), [Observation of second sound in graphite over 200 K, 2022, Eq. 9](https://doi.org/10.1038/s41467-021-27907-z). [On the emergence of heat waves in the transient thermal grating geometry](https://arxiv.org/abs/2111.01578) explicitly demonstrates ballistic/quasiballistic as well as hydrodynamic sign changes. [A transient heat conduction phenomenon to distinguish the hydrodynamic and (quasi)ballistic phonon transport](https://doi.org/10.1016/j.ijheatmasstransfer.2021.121847) studies a temperature undershoot relative to the minimum initial temperature in multidimensional geometries, not the positive-signal two-grating cone here. [Huberman et al. 2019](https://doi.org/10.1126/science.aav3548), [graphite over 200 K](https://doi.org/10.1038/s41467-021-27907-z) and [room-temperature isotope-pure graphite 2026](https://doi.org/10.1038/s41467-026-70807-3) use transient gratings, linewidths/phase or material-specific modeling; these studies are not independently reanalyzed here. [Green's function with full scattering matrix](https://doi.org/10.1103/PhysRevB.104.245424), [dual-relaxation Green's function 2025](https://doi.org/10.1103/PhysRevB.111.035406), and [Interference/thermal-response functions 2026](https://arxiv.org/abs/2607.07644) give broader forward mechanisms, including ones outside this positive-path null. [Callaway approximation can fail in the strongly hydrodynamic regime](https://doi.org/10.1088/1361-648X/ade219) is an explicit guard against treating the model witness as a quantitative material waveform. The prospective addition is the exact stationary-positive-path finite two-q constraint, with sharp constants and an admissible collective violation without relying on either sign change. Neither convex ordering, the Callaway model nor ambiguity of oscillations is claimed new.

## 1. Positive carrier paths imply the convex envelope

Draw a time \(S\) uniformly on \([0,t]\), independently of the whole path. Put \(A=v_S/V\), \(Z=X_t/(Vt)\). Then \(A\) has the calibrated marginal, \(|A|,|Z|\le1\), and \(E[A\mid(v_s)_{0\le s\le t}]=Z\). Conditional Jensen therefore gives a nonnegative functional
\[
 J[\phi]=E\phi(A)-E\phi(Z)\ge0
\]
on all convex \(\phi:[-1,1]\to\mathbb R\). This only needs the one-time marginals, not the Markov property or time-reversal symmetry. For \(\phi_c(x)=-\cos(cx)/c^2\), \(J[\phi_c]=D(c)/c^2\), and \(\phi_c''(x)=\cos(cx)>0\). Moreover, pointwise
\[
 \frac{\cos c_2}{\cos c_1}\le \frac{\cos(c_2 x)}{\cos(c_1 x)}\le1 \quad (|x|\le1);
\]
the ratio decreases for \(x\ge0\), since its log derivative is \(-c_2\tan(c_2x)+c_1\tan(c_1x)\le0\). Thus the two bounding differences of convex generators are themselves convex, and application of \(J\) proves (1).

To see sharpness even with the uniform marginal, take its conditional velocity dynamics to reset at any fixed positive rate only within small velocity intervals near \(+1,-1\), and to leave all remaining velocities invariant. The reset law in each active interval is the equilibrium law conditional on that interval, making the process stationary and reversible. All trajectories and averages in an active interval lie within its convex hull, and give strictly positive Jensen difference in the nonconstant case. Inactive intervals contribute exactly zero. The same curvature comparison restricted to the active intervals squeezes the ratio to the lower constant as those intervals shrink to the endpoint speeds. Resetting in a shrinking interval around zero instead approaches the upper constant. These limiting directions need not give a large signal at the boundary. If irreducibility is desired, add an arbitrarily small positive full-space reset to these positive models, taken sufficiently small relative to the localized reset signal.

**Why ordinary energy-conserving single-mode RTA is in the null.** Write linearized mode energies as heat capacity times the deviation \(g(v)\), including frequency/polarization labels if desired. Its collision term in these variables is
\[
 -L_{\rm RTA}g(i)=-\gamma_i\left[g(i)-\frac{\int \gamma_j g(j)\,d\pi(j)}
                                                   {\int\gamma_j\,d\pi(j)}\right].
\]
This follows by choosing the common RTA relaxation temperature to conserve total energy under the mode-dependent rates. It is exactly the generator of a positive jump process with rate \(\gamma_i\) and reset law proportional to \(\gamma_j\pi(dj)\): its balanced edge weights are symmetric and \(\pi\) is invariant. In the energy-defined readout with the stated initial preparation, Fourier propagation is the expectation over this stationary process (via the usual tilted Markov generator or conditioning on the velocity jumps). The RTA *pseudo-temperature* used inside the collision term is not generally the measured energy component; the inequality does not automatically transfer to that different weighting.

## 2. Physical collective comparator and the exact finite witness

Use an isotropic linear Debye dispersion with a common magnitude \(V\). The projected normalized speed is \(z=v/V\), uniform on \([-1,1]\) under the heat-capacity-weighted solid angle. With \(P_0 g=\langle g\rangle\) and \(P_1g=3z\langle zg\rangle\), take
\[
 L=\Gamma_N(I-P_0-P_1)+\Gamma_R(I-P_0),\qquad
 \Theta(q,t)=\langle 1,e^{-t(L+iq V z)}1\rangle. \tag{3}
\]
This is the gray longitudinal sector of thermalizing Callaway (normal and resistive) collisions, not a postulated damped-wave fitting response. In a constant-rate linear-dispersion Debye gas, the infinitesimal displaced Bose equilibrium in temperature-normalized deviation coordinates is a temperature shift plus a multiple of \(z\): momentum per unit modal energy in the grating direction is \(z/V\). Conservation of energy and directional momentum in the normal relaxation therefore fixes the orthogonal projection onto \(1,z\); the equilibrium resistive collision preserves energy only. Other momentum directions are unexcited spectators by axial symmetry. Each projector complement is positive semidefinite and self-adjoint in the energy/linearized entropy inner product, and the streaming term is skew-adjoint. Equivalently, the full local equilibria are constrained Bose entropy maximizers: separately relaxing occupations towards either of them with a positive scalar rate preserves positivity and increases the concave Bose entropy, by the tangent inequality and their matched invariants. This supplies a positive physical neighborhood for the states whose linearized signed deviations are used here. The linear model is passive and energy conserving, with ordinary resistive relaxation of drift. It is a linearization near positive equilibrium occupations (the displaced Bose state is well-defined for small drift \(<V\)); it does **not** treat a signed excess mode energy as a population probability. The distinction from a positive single-carrier Markov collision is precisely what is tested. Model limitations include exact linear dispersion, gray rates and a bulk coherent energy readout.

In units of the observation time at \(c=qVt\), write \(N=\Gamma_Nt\), \(R=\Gamma_Rt\), and
\[
 A g=-(N+R)g+(N+R)\langle g\rangle+3N z\langle z g\rangle-ic z g.
\]
The generator is dissipative and \(\|A\|\le N+R+c\). Taylor propagating \(A^j1/j!\) yields exact polynomials in \(z\) of degree at most \(j\), coefficients in rational multiples of \(i\) for rational parameters, and exact means from \(\langle z^{2k}\rangle=1/(2k+1)\), odd moments zero. The operator integral remainder for propagating to unit time is bounded by \((N+R+c)^{j+1}/(j+1)!\) because \(\|e^{sA}\|\le1\). Thus the verification evaluates the continuum solution, with explicit rational intervals, no angular discretization or rounding premise. In the same marginal \(B(c)=\sin c/c\), with independent alternating Taylor remainder. Cosines in the cone have similarly exact intervals. This proves (2), while the positivity of both differences rules out a mere comparison to the sign of a single ballistic residual.

For independently bounded errors \(\epsilon_1,\epsilon_2\) on \(D(c_1),D(c_2)\), the lower-cone violation certificate is \(K D(c_1)-D(c_2)>K\epsilon_1+\epsilon_2\), \(K=(c_2/c_1)^2\cos(c_2)/\cos(c_1)\). In the example \(K<9/5\); taking both errors at most \(1/2500\) costs \(<0.00112\), below the verified \(0.0013\) margin. Joint uncertainties in the speed scale, angular law, normalization, and timing need independent propagation in an experiment. Unknown changes of sample or heating/readout as \(q\) varies would not meet this contract.

## 3. Displacement-moment and collision controls

An additional null test can be stated using derivatives at \(q=0\); it is not needed for (2). Let \(M_2(t),M_4(t)\) be the coefficients defined by \(\Theta=1-q^2 M_2/2+q^4 M_4/24+\cdots\). Under the null these are moments of \(X_t\). Let \(\delta_k=t^k\langle v^k\rangle-M_k(t)\). Convexity of \(x^2,x^4,6V^2x^2-x^4\) on \([-V,V]\) gives
\[
 0\le\delta_2,\quad 0\le\delta_4\le6V^2t^2\delta_2. \tag{4}
\]
The coefficient 6 is sharp in the approach to a fixed endpoint velocity, by the same stationary reset control.

For a bounded self-adjoint nonnegative collision operator \(L\) with \(L1=0\) and a real inversion-symmetric normalized response as in (3), the leading regular-time coefficients are
\[
 \delta_2(t)=\tfrac13 t^3 J_2+O(t^4),\quad
 \delta_4(t)=\tfrac15 t^5 J_4+O(t^6),\quad
 J_2=\langle v,Lv\rangle,\quad
 J_4=2\langle v^3,Lv\rangle+\langle v^2,Lv^2\rangle .
\]
These follow by expanding noncommuting words in \(L+iqv\); the first and last factors in each nonzero term are both \(iqv\). The passive Markov bound \(J_4\le10V^2J_2\) also follows directly from the reversible jump Dirichlet form: the edge multiplier of \((v-v')^2\) in \(J_4\) is \(3v^2+4vv'+3v'^2\le10V^2\). The limiting coefficient is sharp for two nearby high same-sign velocities, plus their inverted counterparts. In (3), \(J_2=\Gamma_R V^2/3\) and \(J_4=V^4(4\Gamma_N+22\Gamma_R)/45\). Thus this **short-time moment criterion**, not a hydrodynamic onset threshold, rules out the null for (3) when \(\Gamma_N>32\Gamma_R\).

There is a useful exact limiting contrast at \(\Gamma_R=0,\Gamma_N>0\). Because \(Lv=0\), the \(q^2\) response is exactly ballistic for every time. In the \(q^4\) Duhamel term the outer time gaps are irrelevant (\(Lv=0\)) and the middle gap \(s>0\) contributes \(\langle v^2,e^{-sL}v^2\rangle=m_2^2+(m_4-m_2^2)e^{-\Gamma_N s}<m_4\), where \(m_j=\langle v^j\rangle\). Thus \(\delta_4(t)>0\) for any \(t>0\), even though \(\delta_2(t)=0\). In a positive velocity model, no loss of this second moment means the conditional velocity variance is zero, so **all** displacement moments and the entire grating transfer must be ballistic at that time. A strictly positive resistive rate is used in the finite two-grating witness to avoid relying on this ideal limit.

## 4. Evidence boundary and replay

Run python verify.py for the exact rational continuum propagator, ballistic and cone bounds. No data download or optimization is needed. No current graphite waveform is tested. In particular no compatible same-time calibrated real energy amplitudes at both locked dimensionless \(c\)'s, independently bounded directional group-velocity law and angular support, and preparation/readout-to-energy and optical calibration covariance were verified together. A published normalized temperature dip plot or two independent fits is not a substitute for those nuisance checks. The result is a falsifiable sufficient model rejection rule, not a claim that previously reported material attribution is right or wrong.
