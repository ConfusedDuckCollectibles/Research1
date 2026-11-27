# Static London response: finite scalar aliases, with and without all-mode uniqueness

The finite scalar alias theorem below is complementary to an existing radial Schrödinger uniqueness theorem for the exact harmonic tail. The [tensor all-mode alias](TENSOR.md) does not contradict scalar uniqueness. These are exact conditional results for a known spherical sample boundary in the local, linear, static, vortex-free London model, not a statement about an actual specimen or an inverse estimate from noisy data.

## 1. Model, statement and normalization

Let the specimen be the ball of radius $R$ in vacuum with the same permeability on both sides and prescribed calibrated incident magnetic fields, regular and source-free near the closed ball. Let the unknown $q(r)$ be nonnegative, radial scalar stiffness in $L^\infty$, with units of inverse length squared. For spherical-harmonic degree $\ell\ge1$, a toroidal vector potential solves the problem with constant phase. Its radial factor is $u(r)/r$, where the normalized regular solution satisfies

$
u''=\big[\ell(\ell+1)/r^2+q(r)\big]u,\qquad
u(r)\sim c r^{\ell+1},\quad u(R)=1,\quad p_\ell(q)=u'(R).
$

Indeed $qA$ is distributionally solenoidal, including radial jumps, and the zero-phase toroidal solution equals the unique common forward solution. The incident exterior radial solution is proportional to $r^{\ell+1}-S_\ell R^{2\ell+1}r^{-\ell}$, with the second term carrying the induced field. Matching $u,u'$ gives

$
S_\ell(q)=\frac{Rp_\ell(q)-(\ell+1)}{Rp_\ell(q)+\ell}.
$

This $S_\ell$ is dimensionless and independent of the order within degree $\ell$, not a raw field amplitude or a constant promised for all source normalizations.

**Finite scalar alias.** Fix $0<a<R$, $0<b<Q$ and $q_0=b\,1_{(a,R)}$. For each finite set of distinct positive degrees, and even with fixed $\int_0^R r^2q(r)\,dr$, there are radial $q_\eta$ of strictly positive essential lower bound and at most $Q$, matching all stated data exactly for all sufficiently small $\eta>0$. One can take $q_\eta=\eta$ on $0<r<a$. For **any prescribed $a<r_0<R$** the changes can be confined to $0<r<r_0$, so the two profiles agree on an open outer annulus. The positive core stiffness can be arbitrarily small; the result does not uniformly lower-bound how much can be hidden as the selected list grows. The baseline has active volume $4\pi(R^3-a^3)/3$, versus the full ball for $q_\eta$, where active volume means Lebesgue volume of $\{q>0\}$. Equal integrated stiffness is not equal active volume. Angular decoupling and linearity give equality of the entire exterior field for each calibrated incident combination of the selected degrees.

## 2. Independent proof

Write a regular solution as $r^{\ell+1}w$. With center coefficient $w(0)=1$, integrate $(r^{2\ell+2}w')'=r^{2\ell+2}q w$ from zero. For nonnegative $q$, any first loss of positivity would contradict $w'\ge0$ on the preceding positive interval. After positive normalization this gives $0<u(r)\le(r/R)^{\ell+1}$, and $p_\ell\ge(\ell+1)/R$. For two profiles $q_A,q_B$, the Wronskian $u_A'u_B-u_Au_B'$ has derivative $(q_A-q_B)u_Au_B$ and vanishes at zero. Thus

$
p_\ell(q_A)-p_\ell(q_B)=\int_0^R(q_A-q_B)u_Au_B\,dr.
$

For finitely many amplitudes at fixed interfaces, center data are analytic in a constant core amplitude by the regular power series; on remaining intervals away from zero the usual parameter-dependent fundamental matrix is analytic. The regular value at $R$ remains nonzero in an open real neighborhood of the baseline, including small signed trial amplitudes. These signed values are only a mathematical extension of the data map for the implicit-function theorem. Consequently the derivative of $p_\ell$ is exactly $h\mapsto\int u_\ell(q_0)^2 h\,dr$, and the mass adds the kernel $r^2$.

Fix any open subinterval of $(a,R)$. The kernels $r^2$ and all selected $u_\ell(q_0)^2$ are linearly independent there. To see this, analytically continue each solution from the constant-$b$ shell as a solution of that same equation on $(0,R)$. This is not the physical solution inside the core. Let $y_+\sim r^{\ell+1}$ be the positive regular constant-$b$ solution. Its logarithmic derivative at $a$ is strictly greater than $(\ell+1)/a$, by the integrated equation just used and $b>0$. The hollow-shell solution has exactly $(\ell+1)/a$ at the interface (its physical core is free), so its continuation has nonzero coefficient of the second solution

$
y_+(r)\int_r^a y_+(s)^{-2}\,ds \sim \frac{r^{-\ell}}{2\ell+1}\quad(r\downarrow0).
$

After squaring, the leading singular powers are distinct, $r^{-2\ell}$, with nonzero coefficients. An identically zero kernel combination therefore loses its largest $\ell$ coefficient successively by multiplication by $r^{2\ell}$ and the limit $r\downarrow0$. The remaining $r^2$ coefficient is zero as well. Real analyticity justifies moving from any open shell subinterval to all of $(0,R)$.

There are thus as many distinct radial evaluation points as data, with nonsingular evaluation matrix, in $(a,r_0)$. Replace these by small pairwise disjoint intervals $I_j\Subset(a,r_0)$; their normalized integral columns converge to the evaluation columns, so the unnormalized matrix is also invertible. Set

$
q_{\eta,c}=q_0+\eta\,1_{(0,a)}+\sum_j c_j\,1_{I_j}.
$

The data map has exactly this invertible derivative in $c$. The implicit-function theorem gives $c(\eta)$, with $c(0)=0$ and $c(\eta)=O(\eta)$, preserving the selected $p_\ell$'s and mass exactly. For small positive $\eta$, strict headroom $0<b<Q$ puts every shell segment strictly between zero and $Q$; the core is filled. The same proof works without the mass row. If no degrees are requested, maintaining mass alone is immediate (or use the single $r^2$ row).

**Saturated shell is not a counterexample.** If $b=Q$, fixed mass and any one positive-degree response already distinguish this baseline from every other radial $0\le q\le Q$. Indeed $h=q-q_0$ is nonnegative below $a$ and nonpositive above. In the Wronskian integral use $u_Au_B/r^2=r^{2\ell}w_Aw_B$, strictly increasing for $\ell\ge1$, and subtract its value at $a$ using $\int h r^2=0$. The result is strictly negative for a nonidentical profile. This is a separate obstruction and explains why the finite alias does not extend to the saturated extremizer at the same mass.

## 3. A forward stability bound and what it does not say

Where $q_A-q_B$ is supported in $r\le r_0$, the Wronskian and $u_j(r)\le(r/R)^{\ell+1}$ give

$
|p_\ell(q_A)-p_\ell(q_B)|\le
\|q_A-q_B\|_\infty\,\frac{R}{2\ell+3}(r_0/R)^{2\ell+3}.
$

The derivative of $S_\ell$ with respect to $p_\ell$ is $R(2\ell+1)/(Rp_\ell+\ell)^2\le R/(2\ell+1)$. Therefore

$
|S_\ell(q_A)-S_\ell(q_B)|\le
\frac{R^2\|q_A-q_B\|_\infty}{(2\ell+1)(2\ell+3)}(r_0/R)^{2\ell+3}.
$

Take $r_0=R$ for arbitrary radial nonnegative bounded profiles. The bound applies to the exact aliases as well, with $\|q_\eta-q_0\|_\infty=O(\eta)$ for a fixed selected list. Equalities in those selected degrees are exact; unselected responses are forward-bounded, not asserted zero. This is not a uniform inverse accuracy guarantee, nor does the construction bound the compensator condition number uniformly over the list.

## 4. All sufficiently high scalar modes determine the radial profile

This is a transfer of prior theory, not a new inverse uniqueness theorem. On the **unit** ball put $s=|x|$, $V(x)=R^2q(Rs)$. The regular solution of the scalar zero-energy equation $(-\Delta+V)\Psi=0$ with boundary $Y_{\ell m}$ is

$
\Psi(s,\omega)=\frac{u_\ell(Rs)}s\,Y_{\ell m}(\omega).
$

This follows by calculating the scalar radial Laplacian, or putting $v(s)=u_\ell(Rs)$ so the centrifugal term is unchanged. The corresponding unit-ball Dirichlet-to-Neumann eigenvalue is $\lambda_\ell[V]=Rp_\ell(q)-1$. This encoding does not posit a physical monopole; it simply identifies the positive harmonic coefficients with known spectral samples of a scalar boundary map. Nonnegative $V\in L^\infty$ is safely in the radial class of the theorem below, and zero is not a Dirichlet eigenvalue (Poincaré and the nonnegative potential). No smoothness or agreement near either boundary or center is needed.

[Daudé, Macià, Meroño and Nicoleau, *Stable factorization of the Calderón problem via the Born approximation*, JLMS](https://doi.org/10.1112/jlms.70575) Theorem 1(i) and the paragraph immediately following give unique determination of their radial Born function $V^{\mathrm B}$ from any cofinite eigenvalue tail; Theorem 2 then gives $V$ on each annulus away from zero. In detail, choose a common large $L$ so both weighted Born functions are integrable against $|x|^{2L}dx$. By their formula (1.12), exact equality of $\lambda_\ell$ for $\ell\ge L$ gives zero moments of the weighted radial difference against all $(s^2)^j$, $j\ge0$. Push forward under $s^2$; polynomial density on $[0,1]$ forces this finite signed measure to vanish. Theorem 2 recovers equality of $V$ on every punctured annulus and hence a.e. on the ball. Since the measured $S_\ell$ uniquely gives $p_\ell$, knowing all $\ell\ge1$ is more than enough. Missing $\ell=0$ or any finite number of modes is immaterial for exact injectivity. We make no uniqueness claim for arbitrary infinite sparse degree sets here. [Macià, Meroño and Sánchez-Mendoza](https://doi.org/10.1017/S0308210525100929), especially their radial mode separation and integral identities, are related. [Ramm–Gutman](https://arxiv.org/abs/math-ph/0010033) give a related but different positive-energy scattering near-alias for practically indistinguishable phase shifts, not the stated exact finite calibration and same-mass shell theorem. A positive-energy scattering result is unnecessary for the zero-energy transfer.

## 5. Tensor contrast and controls

For the same ball with directional stiffness, all regular static magnetic excitations admit an exact null demonstrated in [Tensor stiffness](TENSOR.md). If $P_r,P_t$ are radial/tangential projectors, the profiles $K=q_t(r)P_t+q_r(x)P_r$ have the same full interior and exterior magnetic fields for arbitrary measurable bounded nonnegative $q_r$ at fixed radial $q_t$. The toroidal field is tangential and $KA=q_t A$; the weak phase formulation remains valid for singular $K$ without using its inverse. Strictly positive tensor spectral bounds and fixed integrated trace do not remove the all-harmonic ambiguity: sufficiently small smooth mean-zero angularly varying perturbations of the radial eigenvalue on an interior annulus preserve those constraints. This is not a shape or support change and does not imply failure of radial scalar uniqueness. The contrast is physical-model specific: finite calibrated modes cannot distinguish these positive-volume scalar aliases, perfect infinite-mode radial scalar data can separate every nonidentical profile a.e., and direction-dependent stiffness can hide exactly even with such data. No pressure-cell corrections, measured accuracy, finite-size penetration estimate or material-volume inference follows.

[The independent replay](replay-finite-scalar.cjs) reproduces the finite numerical example via independent elementary constant-layer solutions rather than the original Runge–Kutta integrations. It gives a degree-three slope difference about $-1.1034546473\times10^{-5}$ and response difference about $-1.4827851\times10^{-6}$; the decimals alone do not certify equality of the printed offsets. The script also certifies a strictly positive-shell Jacobian for a separate exact $R=1,a=1/2,b=1/25$ control at bands $[4/8,5/8],[5/8,6/8],[6/8,7/8]$, rows mass and degrees 1,2. From $0\le w'/w\le b r/(2\ell+3)$ at the hollow-shell baseline, the kernel $u_\ell^2$ differs from the free power downward by a relative at most $b/(2\ell+3)$ throughout. After clearing positive row denominators, every Leibniz product for this determinant changes by a fraction at most $b(1/5+1/7)=12/875$, so the exact determinant is positive from the free determinant $7238341440$ and the free permanent $510393577830$. This proof does not use a zero-stiffness determinant as a certificate at $b=2$. The analytic independence proof already applies to all strictly positive $b$.
