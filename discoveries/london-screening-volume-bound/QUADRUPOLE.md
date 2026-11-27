# Complement: a held-out nonuniform static field separates the radial equal-dipole pair

The uniform-field alias in PROOF.md is not a complete exterior boundary-operator alias. Two known radial stiffness profiles with a single ordering change can be compared exactly using just one further multipole. This note concerns radial profiles in the same vortex-free local London model, not the general arbitrary-layout identification problem or an existing pressure-cell measurement.

Let q_A and q_B be nonnegative bounded radial profiles in radius R. Assume for some a∈(0,R) that δq=q_A−q_B is nonpositive below a and nonnegative above a, and is not zero a.e. For an applied curl-free regular magnetic solid harmonic of degree ℓ≥1 in the sample region, use the corresponding divergence-free toroidal vector potential angular factor and radial potential u/r. The regular positive radial solutions obey
\[
 u_{\ell,j}''(r)=[\ell(\ell+1)/r^2+q_j(r)]u_{\ell,j}(r), \quad
 u_{\ell,j}(r)\sim c_{\ell,j}r^{\ell+1}\ (r\to0),\quad u_{\ell,j}(R)=1 .
\]
This is the same separation as the ℓ=1 uniform calculation, with the angular vector Laplacian eigenvalue ℓ(ℓ+1). Bounded q near zero gives the regular solution by integrating (r^(2ℓ+2)w')'=r^(2ℓ+2)q w for w=u/r^(ℓ+1), starting at positive w(0), and then normalizing. It is positive everywhere. Radial q preserves this harmonic and φ=0. An applied coefficient is calibrated against the known no-specimen harmonic; near the exterior of the sphere the total radial potential has a regular part proportional to r^(ℓ+1) and an induced decaying part proportional to −S_ℓ R^(2ℓ+1) r^(−ℓ). Let p_{ℓ,j}=u_{ℓ,j}'(R). Exact matching gives S_ℓ=(R p_{ℓ,j}−(ℓ+1))/(R p_{ℓ,j}+ℓ), strictly increasing in p; this normalizes perfect exclusion to 1. ℓ=2 is an applied field with a calibrated linear spatial gradient, not more uniform applied directions.

For each ℓ a cross-profile Wronskian integration yields
\[
 p_{\ell,A}-p_{\ell,B}=\int_0^R \delta q(r)u_{\ell,A}(r)u_{\ell,B}(r)\,dr .
\]
The Wronskian vanishes at 0 by the regular behavior and is p_A−p_B at R. Now suppose the two specimens give the exact same S_k for some k≥1. Within each profile, the derivative of Wronskian between ℓ>k and k is [ℓ(ℓ+1)−k(k+1)]u_{\ell,j}u_{k,j}/r²>0, so u_{\ell,j}/u_{k,j} is strictly increasing in r. Multiply the two ratios, calling the product G(r). Since the k moment vanishes, the ℓ moment is equivalently
\[
 \int_0^R \delta q(r)u_{k,A}(r)u_{k,B}(r)\,[G(r)-G(a)]\,dr>0 .
\]
Both signs agree on each side; strictness follows from nonzero δq on positive measure. Thus A has strictly stronger S_ℓ than B for every ℓ>k. For profiles with the reverse crossing the inequality reverses. This does not require an infinitesimal perturbation and imposes no unattested relation between penetration depth and the measured dipole. The full cone of arbitrary radial or nonradial q is not identified by two harmonics.

For the equal-dipole x=100, s=98 shell/solid example, A=shell and B=solid. replay-quadrupole.cjs gives the finite exact bracket 3.5×10^(−5)<S_2(A)−S_2(B)<3.7×10^(−5), in the stated dimensionless harmonic normalization. A real application would have to independently calibrate the applied gradient, induced quadrupole mixing, position and cell response at sufficient accuracy; this is not what an NV map under uniform excitation does. The model is a possible additional design check, not a published observational comparison or a proof of actual sample fractions.

Separate this comparison from standard spherical multipole separation and Sturm comparison (those tools are classical). Its application explains why the exact uniform-field pair is not generally ambiguous under a calibrated independent static harmonic. The two-excitation upper/lower design is not a claim of optimal general experimental strategy or proof of completely new prior art.
