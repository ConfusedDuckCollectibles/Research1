# Exact coupled stiffness sensitivities at common phase shear

## 1. Model and result

Let a periodic cell in `d >= 2` dimensions contain two perfectly bonded, isotropic linear elastic phases with positive bulk moduli `K1,K2`, an identical baseline shear modulus `mu>0`, and phase-1 indicator `chi` and volume fraction `f in (0,1)`. The two bulk moduli and **all geometry stay fixed**. Only phase 1's shear changes, to `mu+t`; phase 2 has shear `mu` throughout. Units for `t` are those of shear stiffness. Brackets denote the cell mean. For a symmetric constant strain `E`, use

```
W_t(E) = E:C*(t):E
       = min_u < K(x) [tr(E+e(u))]^2
                  + 2(mu+t chi) |dev(E+e(u))|^2 >,
e(u) = (grad u + grad u^T)/2,
```

with mean-zero periodic displacements and no factor of 1/2 in `W`. The norm is Frobenius. Define `H=I/d`, `B(t)=W_t(H)`, `C_x(t)=dev(C*(t)H)` and `S(t)=(1/(2m)) sum_alpha W_t(E_alpha)`, where `m=d(d+1)/2-1` and the `E_alpha` are any orthonormal basis of traceless symmetric strains. These are the strain-controlled hydro stiffness, a hydro-to-shear coupling tensor, and the orientation-invariant mean deviatoric stiffness. In the isotropic case, `B,S` are the usual bulk/shear moduli and `C_x=0`. No effective isotropy is assumed below. In particular `B` is **not** the inverse hydro compliance when there is coupling.

Put
```
v = 1-1/d;  a = 2 mu v;  wi = (Ki+a)^(-1);
h = f w1 + (1-f) w2;   delta = w1-w2;   beta=delta/h;
T = v f(1-f);   Q = -S''(0);
q0 = [2(d-1)f(1-f)/mu + 4 w2 T]/m.
```
If the bulks differ, then the full same-geometry tangent locus, allowing any periodic lattice, is exactly the following set for `b=B'(0)`, `c=C_x'(0)` and `Q`:
```
M = I/d + c/[2 beta f(1-f)] is positive semidefinite,
|c|^2/(2 f) <= b <= 2 beta^2 T - |c|^2/[2(1-f)],
m Q = m q0 + (2 h^2/delta) b - h |c|^2.                (1)
```
Here `c` is symmetric trace-free, so `tr M=1`. If a particular lattice is fixed, the set in (1) is its exact **closure**, and all points in (1) with `M` positive definite are already realized on that lattice. The singular part is realizable on a lattice compatible with the corresponding subspace; exact attainment on every predetermined lattice need not hold. Also
`B(0)=1/h-a, S(0)=mu, C_x(0)=0, S'(0)=f`.
The construction in §3 can use arbitrarily small inclusions and is an exact continuum result; finite-scale approximations are possible, not claims of a manufacturing specification. For identical bulks use (3) below without dividing by `delta`; the formula then gives both slopes zero for `B,C_x` and universal `Q=q0`.

For example, when the lower-bulk phase is the one stiffened, `delta>0`. At zero first-order coupling, the larger the first-order gain in bulk stiffness, the larger the second-order deficit from the mean-shear linear tangent. This is not a claim that shear stiffness decreases.

## 2. Derivation

Set `D = dev Hess Delta^(-1)(chi-f)` on the torus (the inverse fixes zero mean), `J=<chi |D|^2>` and `R=<chi D>`. Its nonzero Fourier multiplier is `n n^T-I/d` with unit reciprocal-lattice direction `n`. Parseval gives `<|D|^2>=T` and
```
I/d + R/[f(1-f)] = sum_{k ne 0} |chi_hat(k)|^2 n n^T/[f(1-f)].
```
This is a convex average of projectors. Since `<D>=0`, Cauchy–Schwarz on the two phases gives
```
|R|^2/f <= J <= T-|R|^2/(1-f).                          (2)
```

At baseline, the periodic-gradient longitudinal and divergence-free transverse parts of displacements are energy orthogonal. The trace term acts only on the longitudinal displacement. The cross term of the integrated shear vanishes by Fourier orthogonality because `mu` is constant. Write a longitudinal strain as `Hess Delta^(-1) omega` with `<omega>=0`; its baseline quadratic energy is `<(K+a)omega^2>`, since the integrated squared deviatoric part is `v<omega^2>`. The transverse energy is `mu <|grad u_trans|^2>`. The hydro load's optimal total trace is `w/h` (`w=w2+delta chi`), by a scalar weighted minimization with trace mean 1. Thus the optimal hydro total strain has deviatoric part `beta D` and no transverse displacement. The optimized hydro energy is `1/h-a`. The optimal deviatoric load at baseline is just the constant strain because its stress `2mu E` is constant. The hydro/shear effective cross term therefore vanishes.

Let the quadratic forms on the corrector-strain Hilbert space be `a_0` for baseline and `a_1(U,V)=2<chi dev U:dev V>` for the change. Use `a_0` as inner product and let `A` represent `a_1` on that space. Then `0 <= A <= I/mu`. If `eps_E^0` is the baseline minimizing total strain, define `g_E` as the Riesz representative of `a_1(eps_E^0,.)`. Baseline orthogonality and completing the square give the **exact** identity (for example at all `t>=0) and in a neighborhood of zero)
```
W_t(E) = W_0(E) + t a_1(eps_E^0,eps_E^0)
                    - t^2 <g_E,(I+t A)^(-1)g_E>_{a_0}.
```
Consequently
```
S'(0)=f,    b=2 beta^2 J,    c=2 beta R,
m Q = 2(d-1)f(1-f)/mu + 4 w2 T + 4 delta J
                                - 4 delta^2 |R|^2/h.  (3)
```
To check the last line explicitly, on the longitudinal subspace the functional for `E_alpha` is `2<F_alpha omega>`, where `F_alpha=D:E_alpha` by self-adjointness of the Hessian multiplier. The weighted constrained squared Riesz norm is `4[<w F_alpha^2>-<w F_alpha>^2/h]`. Sum over alpha, use `sum F_alpha^2=|D|^2` and `<w F_alpha>=delta R:E_alpha`. On the transverse modes the squared Riesz norm is `(4/mu) sum_{k ne 0} |chi_hat(k)|^2 |P_{n^perp}E_alpha n|^2` by integration by parts. For each `n`, the sum over alpha is `(d-1)/2` (rotate `n` to the first coordinate and use the normalized off-diagonal strains). The sum of these longitudinal and transverse norms squared is `mQ`, with the stated `W,S` normalization. Equations (2)–(3) now imply (1) as necessary.

## 3. Sharpness, including anisotropy

It suffices to realize (2) for each prescribed positive-semidefinite trace-one `M`, with `R=f(1-f)(M-I/d)`. First suppose `M` is positive definite, rotate to its eigenbasis, and consider a confocal ellipsoidal core and envelope. Write the core squared semiaxes as `a_i^2>0` and the outer squared semiaxes as `a_i^2+s`, so the inner/outer volume ratio is
`g=prod_i(1+s/a_i^2)^(-1/2)`. The construction below works for any `g in (0,1)`.

**Compact scalar potential.** Put `P(u)=prod_i sqrt(a_i^2+u)`. Inside the outer ellipsoid let `lambda=0` in the core and otherwise solve `sum_i x_i^2/(a_i^2+lambda)=1` for `lambda in (0,s)`. Define
```
Phi(x) = -P(0)/4 * int_lambda^s
                 [1-sum_i x_i^2/(a_i^2+u)]/P(u) du,
```
and set it to zero outside. This is the usual signed confocal Newtonian potential difference, here written as a finite integral valid also in dimension 2. Its gradient components are `[P(0)x_i/2] int_lambda^s 1/[(a_i^2+u)P(u)] du`: the derivative of the endpoint contributes zero by the equation for `lambda` in the shell. The value and gradient are continuous at both faces and vanish at the outer face; it has piecewise bounded second derivatives and belongs to `H^2`. By the derivative of `1/P`, the Laplacian in the core is `1-g`. In the shell the same integral gives `P(0)/P(lambda)-g` before differentiating the lower limit. The missing endpoint term is exactly `-P(0)/P(lambda)`, since implicit differentiation gives `sum_i [x_i/(a_i^2+lambda)] partial_i lambda=2`. Thus distributionally `Delta Phi=1_core-g 1_outer`.

On the core, the Hessian difference is diagonal with components
```
(1-g) M_i(p)
  = (1/2) int_0^s [p_i/(1+u p_i)] prod_j(1+u p_j)^(-1/2) du,
p_i = 1/a_i^2,    prod_j(1+s p_j)=g^(-2).                 (4)
```
Rescale axes to take `sum p_i=1`. For `p` in the closed simplex, the positive `s` is continuous, finite and uniformly bounded (the product is at least `1+s`). The map `p -> M(p)` is continuous; its coordinates are nonnegative, they sum to one by integrating the product derivative, and `M_i(p)=0` iff `p_i=0`. It maps onto the simplex. In detail, for a positive target `y` apply Brouwer's theorem to the map having coordinates `(p_i+y_i-M_i(p))_+` normalized to sum one. Its denominator is at least one. At any fixed point all coordinates `p_i` must be positive, as a zero coordinate would map to a positive one. Then no truncation occurs and its denominator is exactly one, so `M(p)=y`. It would suffice to work in a corresponding face for a boundary target.

**Both boundaries and filling their interval.** For given aspect ratio, orientation and `g`, fill a cell (or any open subregion) almost everywhere by a disjoint countable collection of scaled similar outer ellipsoids. Existence follows by Vitali covering (each fixed shape has bounded eccentricity), or by finitely many inscribed, uniformly proportioned ellipsoids in each dyadic cover stage and then exhaustion. Periodization can be performed with ellipsoids wholly inside the cell. Each potential has zero exterior value and gradient. Under scaling by `ell`, its Hessian is unchanged in shape and its squared `L^2` norm scales as volume. The summed potential is periodic `H^2` by disjoint supports and bounded cell size, and its Laplacian is `chi-f` when phase 1 is the union of cores and `g=f`. Its core Hessian is `(1-f)M`, so `R=f(1-f)(M-I/d)` and `J=|R|^2/f`. To realize the other boundary, choose `g=1-f` and let phase 2 be the union of cores. The potential difference changes sign, `D` in phase 2 is `-f(M-I/d)`, and now `J=T-|R|^2/(1-f)` with the **same** `R`. Since the compact potential works in any open subregion and both recipes have the same local fraction and `R`, one can partition a cell into two open pieces of prescribed volume up to a null interface and apply the two recipes independently. This gives every interpolation in `J` at fixed `R`. The resulting Hessian consists of the local fields without any limit argument for interaction between the two regions.

If `M` has rank `r<d`, carry out the corresponding construction on its range with an adapted lattice and extend the indicator and scalar potential constantly in the remaining dimensions. The Hessian in the varying coordinates then equals the required core Hessian, with zeros in the others. The trace-free operation is still for the full dimension `d`, and all averaging identities hold as above. When `r=1`, a simple two-phase layer suffices; the interval in (2) collapses at that `M`. On any *fixed* lattice, approximate a singular `M` by positive definite `M_e=(1-e)M+e I/d` and use the same interpolation parameter in (2). Necessity and the continuity of the locus give the closure assertion. In particular equality at the maximally anisotropic rank-one `M=n n^T` forces all nonzero Fourier modes to have directions parallel to `n`, which is not possible on a lattice with no such reciprocal direction.

The confocal potential and uniform inclusion fields are classical (see SOURCES.md). The use here is to preserve a **single geometry during shear change** and realize each joint sensitivity combination, rather than to solve for all static components of the stiffness of an ellipsoid packing.

## 4. Finite bounds and positive-phase check

For `t>0` define `b_t=(B(t)-B(0))/t` and `q_t=2(mu+f t-S(t))/t^2`. The resolvent identity gives
```
Q/(1+t/mu) <= q_t <= Q,
b/(1+t/mu) <= b_t <= b.                                  (5)
```
For the second line, the exact hydro increment is the minimum of `||eta||_0^2+t||eps_H^0+eta||_1^2`. By the triangle inequality for the second seminorm and `||eta||_1<=||eta||_0/sqrt(mu)`, it is bounded below by the scalar minimum of `x^2+t(sqrt(b)-x/sqrt(mu))_+^2` over `x>=0`, which is `t b/(1+t/mu)`. Setting `eta=0` gives the upper bound. These are valid envelopes, not claimed to be sharp finite-`t` static bounds.

In the predeclared dimension-2 case `f=1/2, mu=1, K1=1/10, K2=10, t=1/2`, `delta>0`. To bound `S(t)` from a measurement of `B(t)` without assuming symmetry, apply `|R|^2 <= (1-f)(T-J)` in (3), and then use `J >= b_t/(2 beta^2)`. The resulting lower bound on `Q` is increasing in `J`. Thus `B(1/2)>=13/10` implies `Q>=637/1452` and
`S(1/2)<=21143/17424 = 1.2134412305...`, excluding a simultaneous `S>=1.214`. The exclusion survives separate hard measurement errors of up to `1/4000` on each modulus (not errors in constituent parameters, geometry or intervention).

Each goal is individually achievable. Fill with coated disks having phase 1 as the core; `J=R=0` and `S(1/2)>=107/88>1.214` from (5); also `B=1` since `b=0`. In the reverse disk geometry phase 1 is the shell. Its exact hydro modulus is the common-shear Hill formula with shear set to the **shell** shear `3/2`, i.e. `343/262>1.30`. Indeed the common-shear hydro field is pure dilation on the core, so restoring the other shear there leaves the actual stress and the weak equilibrium unchanged.

These outcomes can be attained with finitely many interfaces. Keep disks covering all but a fraction `eta<=1/10000` of a square cell, then refill the remaining set to preserve the exact phase fraction. One finite choice is to take in each step a finite dyadic square cover of at least half the remaining open area (squares wholly inside the remaining region, disjoint interiors), then inscribe a disk in each; after 20 steps the residual is at most `(5/8)^20<1/10000`. Boundaries are null; the leftover phase fraction can be set by a straight cut in the residual. Compare with a continuation of the same disk packing to full measure. For the phase-1-core case the changed indicator lives on at most `eta`; the operator `dev Hess Delta^(-1)` has norm `sqrt(v)`, and in the ideal comparison `|D|<=sqrt(v)` pointwise (the radial shell factor is at most one). Hence `sqrt(J_new)<=2 sqrt(v eta)` and `J_new<=2 eta`. Apply (3)–(5) to get `S_new>=107/88-9 eta/22`, still above the target even with the hard reading error. For phase 1 as shell use the equilibrium ideal stress as a dual test. The exact radial hydro traces are `230/131` in the shell and `32/131` in the core. The ideal stress satisfies `|sigma|^2<12`: its hydro and deviatoric contributions in the shell sum to at most `2(23/131)^2+(9/2)(198/131)^2`, and in the core it is pure dilation with norm squared `2(320/131)^2`. The norm of the local compliance change is less than 5 by the hydro and deviatoric eigenvalues. Using `2<sigma>:H-<sigma:C_new^(-1):sigma>` gives `B_new>=343/262-60 eta>1.303`, including a margin for the reading error.

The finite test is an example, not an improved static bound: at the exact high bulk edge `B=343/262`, Bogdan's recent sharper **static** calculation for isotropic effective tensors caps `S` at `387/322`; his opposite bulk-edge lower bound is `61/50` (see SOURCES.md). This does not numerically evaluate the optimal interior static bound at `B=1.30`. No novelty is claimed for ruling out these illustrative finite goals by themselves.
