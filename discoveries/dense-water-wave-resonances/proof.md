# Dense active triad depths and a cubic small-divisor obstruction on the square water-wave torus
AI-assisted Research1 note, November 25, 2026. Not peer reviewed. No priority, experiment, energy-cascade or PDE-lifespan claim.

## Result

Fix any positive dimensionless capillarity alpha (not necessarily algebraic), horizontal periods 2 pi by 2 pi, positive unit gravity, a flat bottom at normalized depth h>0, and inviscid irrotational linear dispersion

    f(r)=sqrt(r(1+alpha r^2)),    omega_r(h)=f(r) sqrt(tanh(r h)).
    w(h)=omega_1(h).

Consider only the integer triads

    k1=(1,0), k2=(m,b), k3=(m+1,b)=k1+k2,  b>=1, m>=0.
    F_(m,b)(h)=omega_|k3|(h)-omega_|k2|(h)-w(h).

1. The positive depths with an exact triad of this form and a **nonzero canonical cubic coefficient** are countable dense. Every such triad has at most one positive root, and that root is simple. For any compact interval I=[h_-,h_+] in (0,infinity) of positive length, the number of these roots at a fixed sufficiently large b, counted as distinct m values, is

       (2/(3 sqrt(alpha))) [w(h_+)-w(h_-)] sqrt(b) + O_(I,alpha)(1).

   For fixed sufficiently large b, those roots on I are distinct. Summing b<=B gives (4/(9 sqrt(alpha)))[w(h_+)-w(h_-)] B^(3/2)+O_(I,alpha)(B) participating pairs (m,b), with no claim of distinctness across different b. One can reach any target depth along a sequence with errors O(1/b). At those roots and at the near-resonances below the canonical cubic coefficient is of order sqrt(b), not a null interaction.

2. A residual subset of (0,infinity) of **Hausdorff dimension zero** has no unbalanced exact frequency relation at any finite order, yet infinitely many triads above with 0<|F|<exp(-b). On this residual set no continuous homogeneous cubic scalar polynomial on a finite-order Fourier Sobolev phase space can solve the full cubic coefficient equation on the restrictions to finitely supported amplitudes. This is an obstruction to that scalar polynomial generator, not to general or partial/quasi-resonant transformations, nor to the PDE flow.

3. In contrast, for Lebesgue almost every fixed depth and every real p>3/2, there is c(h,p,alpha)>0 such that |F_(m,b)(h)|>=c b^(-p) for every pair in the displayed family. This is not a sharpness claim at p=3/2, nor a metric result for every possible wave triad.

If alpha is algebraic, all exact triad depths in (1) are transcendental by the existing [fixed algebraic-depth result](https://github.com/ConfusedDuckCollectibles/Research1/tree/79bb3017701a00fb83a68686a62bd39c9d4d16b3/discoveries/gravity-capillary-resonance); no bound at a specified algebraic depth follows here. “Mixed” means both gravity and capillarity are present. Only for 0<alpha<1 is the unit mode strictly gravity-dominated while its two large partners are capillary-dominated.

## Global uniqueness (known ratio mechanism)

For any fixed ordered radii r1<=r2<r3, real response factors f(r)>0 independent of h, consider S(h)=(omega_r1(h)+omega_r2(h))/omega_r3(h). The logarithmic h derivative of sqrt(tanh(r_i h)/tanh(r3 h)) is

    r_i/sinh(2 r_i h) - r3/sinh(2 r3 h) > 0             (i=1,2).

Indeed r/sinh(2rh) is strictly decreasing at r>0, since x cosh(x)>sinh(x) for x>0. Thus S'>0, a root S=1 is unique, and there F'=-omega_r3 S'<0. The limiting values are respectively (f(r1)sqrt(r1)+f(r2)sqrt(r2))/(f(r3)sqrt(r3)) and (f(r1)+f(r2))/f(r3). A finite positive root exists iff these strictly straddle 1; equality at 0 or infinity alone is not an interior root. This is the same highest-frequency normalization as [Durey–Milewski, Appendix A](https://arxiv.org/html/2301.02163v2) for pure-gravity cylinders, and the monotonicity extends immediately to positive fixed amplitude factors. No priority claim on it. Our family has |k3|>|k2|>=1.

## Signed high-frequency difference: uniform estimates, roots and exact count

Write M=m+1/2, r_j=|k_j|, kappa=(3/2)sqrt(alpha). Then

    r3-r2=2M/(r2+r3),
    f'(r)= (1+3 alpha r^2)/(2 sqrt(r(1+alpha r^2)))
          = kappa sqrt(r) (1+O_alpha(r^-2)).

Fix compact I=[h_-,h_+] with h_->0. Since w>0 and w'(h)>0, both have positive lower bounds there. Corrections from finite depth in the two high modes and their first h derivatives are O_I,alpha(b^(5/2) exp(-c_I b)), **uniform in all m>=0**, because all radii here are at least b and sup_(r>=b) r^A exp(-2h_- r) admits that bound for large b. This follows directly from 1-tanh(x)<=2 exp(-2x) and the derivative of sqrt(tanh(x)). In particular F'=-w'(h)+o(1), so F is uniformly strictly decreasing on I for every m for all large b.

To justify excluding distant m, for b large r3/r2<=2; elementary mean value bounds with f'~kappa sqrt(r) show

    c_alpha (m+1)/sqrt(r2) <= f(r3)-f(r2) <= C_alpha (m+1)/sqrt(r2).

Choose a small fixed delta_I>0 less than w(h_-)/2. If |F|<delta_I anywhere on I, then the deep difference is bounded above and below by positive constants. For m>=b, r2 is comparable to m and the lower estimate grows like sqrt(m), impossible for b large. For m<b, r2 is comparable to b, and both estimates force c_I sqrt(b)<=m+1<=C_I sqrt(b). Thus all possible roots or sufficiently small phases are in this band. The same estimate shows that for any fixed b only finitely many m have a root or |F|<delta_I in I (the exponentially small correction vanishes as m->infinity). Hence using an optional cutoff m<=b changes only finitely many small b for a fixed I.

Uniformly in 0<=m<=C_I sqrt(b), r2,r3=b+O_I(1); substituting into the derivative and radial-difference formula gives the signed phase expansion

    F_(m,b)(h) = kappa (m+1/2)/sqrt(b) - w(h) + O_(I,alpha)(b^-1)       (1)

uniform for h∈I. The exact difference can alternatively be written as (r3-r2)[1+alpha(r3^2+r3*r2+r2^2)]/[f(r3)+f(r2)]. The errors in (1) do not hide a term of order b^-1/2.

**Constructive closeness.** For target h0>0, take compact I containing it in its interior and choose b_m the nearest integer to (kappa(m+1/2)/w(h0))^2. This is positive and grows like m^2. The rounding causes O(1/b_m) change in kappa M/sqrt(b_m). Thus F(h0)=O(1/b_m) while the two compact endpoints have opposite signs for large m, by w strictly increasing. Uniform derivative separation implies a root with |h_m-h0|=O(1/b_m). The amplitude result below shows it is active for all sufficiently large m.

**Counting including endpoint roots.** Let v_m=kappa(m+1/2)/sqrt(b), e_b=C/b bounding (1), and w_±=w(h_±). Any m with root on the *closed* I has v_m∈[w_-−e_b,w_++e_b]. Conversely, each m with v_m∈[w_-+e_b,w_+−e_b] has weakly opposite endpoint signs and hence a root on closed I. To apply (1) to the converse, such m are automatically in the same O(sqrt b) band (enlarge the constant if needed); they are also positive and <=b when b is large. Each endpoint strip of width 2e_b contains only O(1) lattice positions, since their separation is kappa/sqrt(b). Counting consecutive grid points proves the formula in (1) of the result, including accidental endpoint zeros without excluding them. For m values in this band, F_(m+1,b)-F_(m,b)=kappa/sqrt(b)+O(1/b)>0 uniformly; comparing any pair of ordered m's gives distinct roots, increasing with m. The sum follows from Σ_{b<=B}sqrt(b)=(2/3)B^(3/2)+O(sqrt(B)) by integral comparison. This is a count of this subfamily only, with no unsupported exclusivity across all momenta.

## Cubic coupling and its high-frequency asymptotic

The first Dirichlet–Neumann shape derivative gives H3=(1/2)∫eta(|grad psi|^2−(G0 psi)^2), using the flat-bottom Fourier symbol lambda_i=r_i tanh(r_i h), normalized Haar measure and canonical variables

    E_i=sqrt(lambda_i/(2 omega_i)), rho_i=omega_i/lambda_i,
    eta_k=E_i(a_k+conj(a_-k)),
    psi_k=−i E_i rho_i (a_k−conj(a_-k)).

The exact signed finite-depth coefficient is already in [Chabane–Choi (2019), Eqs. (12), (13), and (37)](https://web.njit.edu/~wychoi/pub/SAPM_19_CC.pdf); see the independent audit for the direct relabeling. The step here is its uniform asymptotic in this specific lattice family. Capillary potential energy is an even function of grad eta and has no cubic term. The three directed Fourier momenta for a1 a2 conj(a3) are k1,k2,−k3. Treat positive/negative psi phases as −i,+i. For each choice of the eta factor, the two psi permutations cancel 1/2; the undecorated kernel is −p·q−lambda_p lambda_q. Their sum, divided by the positive common E1 E2 E3, is

    C = rho1 rho2 (m+lambda1 lambda2)
      + rho1 rho3 (m+1−lambda1 lambda3)
      + rho2 rho3 (Q−lambda2 lambda3),
    Q=k2·k3=m(m+1)+b^2.

There is a controlled cancellation in each large term. Combining the two terms with rho1 yields the exact identity
C = rho1 [rho2 m+rho3(m+1)−lambda1(w+F)] + rho2 rho3(Q−lambda2 lambda3).
The large pair obeys Q−r2 r3=−b^2/(r2 r3+Q) by the two-dimensional cross-product identity. Consider a sequence in I with b->infinity and F->0. The band and (1) imply M/sqrt(b)=w(h)/kappa+o(1), uniformly as these two limits hold. Furthermore rho2,rho3 ~ sqrt(alpha b) and Q−lambda2 lambda3=−1/2+o(1) uniformly; finite-depth factors in this expression incur only O(b^2 exp(-c_I b)). Since rho1 w=1+alpha exactly, this gives

    C/b = (4/3)(1+alpha)−alpha/2+o(1) = 4/3+5alpha/6+o(1).

It applies at the exact roots and for near zero F, not only at a signed zero. Finally E2 E3 sqrt(b)->1/(2 sqrt(alpha)), and E1(h) is continuous and strictly positive, bounded away from zero on I. Therefore the actual canonical coefficient H_{123}=E1 E2 E3 C obeys

    H_{123}/sqrt(b) = E1(h)/(2 sqrt(alpha)) [4/3+5alpha/6] + o(1).

In particular the resonant monomials are genuine. This is an algebraic nonvanishing result, not a statement of inevitable energy exchange or a theorem about interacting amplitudes at fixed nonzero steepness.

## All-order deletion, metric contrast and precise obstruction

For a fixed finite unbalanced signed list, group by distinct positive radii r, so that some coefficients n_r in Σ n_r omega_r are nonzero. The result is real analytic on h>0. Its expansion at large h is
Σ n_r f(r) − n_{r_min} f(r_min) exp(−2 r_min h)+o(exp(−2 r_min h)),
where r_min is the smallest radius with nonzero grouped coefficient. This proves it is never identically zero, whether or not the constant sums to zero and regardless of alpha algebraicity. Thus its zero set is discrete in (0,infinity) and countable (a compact exhaustion suffices). There are countably many finite lattice lists and sign choices. Their union E is countable; it does not remove always balanced shell combinations. In particular outside E no exact odd-order resonance of nonzero modes is possible, because equal per-shell positive and negative counts would make total degree even. This is a statement about finite lists, not convergence of infinite sums.

For every threshold B the union over b>=B,m>=0 of the open sets |F_(m,b)(h)|<exp(−b) is dense: for any open interval, the construction placed a root at arbitrarily large b in its interior. The Baire category theorem and removal of countable E give a residual set S with infinitely many strictly positive magnitudes smaller than exp(−b). The same construction works for any chosen positive b-dependent tolerances; no uniform lower tolerance is needed for Baire density.

For the measure side, on I all sufficiently large b have O(sqrt b) relevant m by the band estimate. Uniform |F'|>=c_I gives each nonempty sublevel |F|<epsilon_b a diameter at most 2epsilon_b/c_I once epsilon_b->0. At epsilon_b=exp(−b), the sum of s-powers of all lengths converges for every s>0 and tends to zero on tail unions. Thus S intersect I has Hausdorff dimension zero, and so does S by compact exhaustion. For epsilon_b=b^(−p), sum sqrt(b) b^(−p) converges when p>3/2; Borel–Cantelli gives eventual lower bound b^(−p) almost surely on I. For each of the finitely many remaining b, only finitely many m have bounded phase on I; excluding E gives a positive minimum after multiplying by b^p. Compact exhaustion and a countable sequence p↓3/2 prove the assertion simultaneously for all p>3/2. This proves only a sufficient metric exponent bound in this one family.

**Formal normal-form obstruction made precise.** Let h in S. Write h^s for the canonical amplitudes with finite-order Fourier norm (Σ⟨k⟩^(2s)|a_k|^2)^(1/2). If a continuous homogeneous cubic scalar polynomial χ3 on this real phase space had the Taylor coefficients of a generator solving the full cubic homological equation, its coefficient on a1 a2 conj(a3) would be of magnitude |H_{123}|/|F|. This follows from the diagonal quadratic Hamiltonian multiplying that monomial by ±iF in the Poisson bracket. On our subsequence it is at least c sqrt(b) exp(b). But boundedness of a continuous homogeneous scalar cubic on norm-bounded sets follows by homogeneity from continuity at zero. Restrict to the three independent complex modes of radii 1, r2, r3, choose their amplitudes ⟨k_j⟩^(−s) and average over the three independent phases against exp(−i(θ1+θ2−θ3)). This isolates the coefficient times the amplitude product and is bounded by that bounded-set norm. Thus every such coefficient is O_s(b^(2s)), a contradiction for each finite s. Equivalently a finite bounded third derivative/trilinear form cannot generate the full elimination. This does not exclude retaining near-resonances, nonpolynomial methods with different regularity hypotheses, or PDE solutions.

## Exact interval checks and the zero-capillarity negative control

The verifier encloses every radical, low-mode tanh and sign with exact rational intervals. For high arguments >=200 use 1−2^(-399)<=tanh(x)<=1 from e>2; no floating cancellation is needed. For both fixed checks it independently enumerates all six directed momentum placements in H3 and compares the algebraically stabilized coefficient. The roots are unique by the ratio proof, not by a numerical stopping rule.

| alpha | m,b | target | certified depth interval | C/b approximately |
|---:|---:|---:|---:|---:|
| 1/2 | 100, 14105 | 0.6 | (0.5999911374537900, 0.5999911374537902) | 1.75000001381028 |
| 2 | 100, 101759 | 0.15 | (0.1499997535178385, 0.1499997535178388) | 2.99999804934961 |

The first bracket corrects the last decimals of the preliminary floating Dreamer note; the correction is well within the stated asymptotic and does not change the choice b. Under alpha=0 the known pure-gravity result is a clean negative control: omega_r/r=sqrt(tanh(rh)/r) is strictly decreasing in r>0 while omega_r is increasing. Strict subadditivity plus the vector triangle inequality rules out every nonzero pure-gravity momentum triad at any positive depth on this flat torus. The code separately tests its selected high-b sign.

See SOURCES.md for exact scope and prior-art comparison. Run:

    sha256sum -c protocol.sha256
    python3 -B verify.py > /tmp/dense-depth-replay.txt
    cmp /tmp/dense-depth-replay.txt replay.txt
    sha256sum -c manifest.sha256

No dataset or external action; all mathematical controls are private.
