# Fixed-depth rounding law for one family of gravity–capillary triads

This is an addendum to [the dense-depth calculation (#21)](https://github.com/ConfusedDuckCollectibles/Research1/blob/ccd13a0b5f0eeed7d0cfc7b1caa8480dafcc3bdb/discoveries/dense-water-wave-resonances/proof.md), not a statement about all water-wave triads. Fix normalized depth \(h>0\), unit gravity, capillarity \(\alpha>0\), and the square torus with periods \(2\pi,2\pi\). Write
\[
 \omega_r=\sqrt{r(1+\alpha r^2)\tanh(rh)},\quad
 w=\omega_1,\quad \kappa=\tfrac32\sqrt\alpha,\quad A=(\kappa/w)^2,\quad M=m+\tfrac12.
\]
For integers \(m\ge0,b\ge1\), set \(r_2=|(m,b)|\), \(r_3=|(m+1,b)|\), and \(F_{m,b}=\omega_{r_3}-\omega_{r_2}-w\). The vector momenta obey \((1,0)+(m,b)=(m+1,b)\). Define the corrected real rounding center
\[
 g(m)=A M^2-\frac1{2A}.
\]

**Fixed-depth results.** When \(b_m\) is a nearest integer to \(g(m)\) (positive for all large \(m\)), the signed scaled phase is
\[
 \frac{b_m F_{m,b_m}}{w}=\frac{g(m)-b_m}{2}+O_{\alpha,h}(b_m^{-1}).
 \tag{1}
\]
For irrational \(A\), these values converge in empirical distribution to uniform measure on \([-1/4,1/4]\), and \(\liminf_m b_m|F_{m,b_m}|=0\), \(\limsup_m b_m|F_{m,b_m}|=w/4\). For every \(\varepsilon>0\) infinitely many pairs in the same rounded family satisfy
\[
 0<F_{m,b_m}\ll_{\alpha,h,\varepsilon} b_m^{-6/5+\varepsilon}.
 \tag{2}
\]
For every fixed \(\alpha>0\) and Lebesgue-a.e. \(h>0\), infinitely many such pairs satisfy, for any \(\eta>0\),
\[
 0<F_{m,b_m}\ll_{\alpha,h,\eta} b_m^{-5/4}(\log b_m)^{3/2+\eta}.
 \tag{3}
\]
The almost-everywhere result is separate from the all-irrational result, not assumed at a specified algebraic pair.

For rational \(A=P/Q\) in lowest terms, \(g(m)\bmod1\) is periodic and avoids both integers and half-integers. The limiting law in (1) is the finitely supported law given by the \(Q\) residuals, with their multiplicities. In the **entire** displayed triad family,
\[
 \lim_{R\to\infty}\inf_{\substack{m\ge0,\,b\ge1\\m+b\ge R}} b|F_{m,b}|
 =\frac w2\min_{0\le j<Q}\operatorname{dist}(g(j),\mathbb Z)>0.
 \tag{4}
\]
Only finitely many exact zeros can occur in this family at such a depth, and after their exclusion a positive uniform lower bound \(b|F|\ge c_{\alpha,h}\) holds. Neither upper exponent is claimed optimal, and neither conflicts with the #21 eventual a.e. lower bound at exponents \(p>3/2\).

## The correction term and scope of the remainder

For \(r=(b^2+x^2)^{1/2}\), first take the unit difference in \(x=m,m+1\) of the deep capillary term
\[
 \sqrt\alpha\,r^{3/2}
 =\sqrt\alpha\,b^{3/2}(1+x^2/b^2)^{3/4}.
\]
In the range \(0\le m\le C\sqrt b\), the Taylor coefficients are \(1,3/4,-3/32\). Writing \(t=x^2/b^2\), the omitted Taylor term has \(t\)-derivative \(O(t^2)\), so its \(x\)-derivative after multiplication by \(\sqrt\alpha\,b^{3/2}\) is \(O_\alpha(x^5 b^{-9/2})\). It is important to bound a derivative, rather than subtract two unrelated \(O(t^3)\) estimates. The gravity correction \(\sqrt{r(1+\alpha r^2)}-\sqrt\alpha\,r^{3/2}\) has \(r\)-derivative \(O_\alpha(r^{-3/2})\), so its \(x\)-derivative is \(O_\alpha(x b^{-5/2})\) here. Each resulting unit-difference remainder is \(O_{\alpha,C}(M b^{-5/2})\). Finite depth changes the two high frequencies by at most \(O_{\alpha,h}(b^{3/2}e^{-2hb})\) on this band.

The exact differences of the quadratic and quartic terms in \(x\) are \(2M\) and \(4M^3+M\). Thus
\[
 \omega_{r_3}-\omega_{r_2}
 =\frac{\kappa M}{\sqrt b}\left(1-\frac{M^2}{4b^2}+O_{\alpha,h,C}(b^{-2})\right),
 \tag{5}
\]
uniformly on the indicated band; the \(+M\) part of the quartic difference contributes \(-1/(16b^2)\) to the parentheses. With \(d=b-AM^2\) bounded, \(M^2=(b-d)/A\) and \(\kappa M/(w\sqrt b)=\sqrt{1-d/b}\). Substitution gives
\[
 \frac{F_{m,b}}{w}=-\frac{d+1/(2A)}{2b}+O_{\alpha,h,\sup|d|}(b^{-2})
 =\frac{g(m)-b}{2b}+O_{\alpha,h,\sup|d|}(b^{-2}).
 \tag{6}
\]
The absolute \(O(b^{-2})\) remainder in (6) requires bounded \(d\), not simply comparability of \(b\) and \(Am^2\). Before centering, the relative remainder in (5) holds throughout \(m=O(\sqrt b)\). Both statements can be taken locally uniform on compact positive parameter ranges with the other indicated bounds fixed.

For (4), consider **any** sequence escaping all finite boxes along which \(b|F|\) remains bounded. If \(b\) is bounded on a subsequence, \(m\to\infty\) there and the deep high difference is \(\asymp_\alpha\sqrt m\) while the finite-depth corrections vanish, impossible. Thus \(b\to\infty\). More generally the deep difference is comparable to \((m+1)/\sqrt{m+b}\) for large radii: use \(f'\asymp_\alpha\sqrt r\) with \(f(r)=\sqrt{r(1+\alpha r^2)}\) and \(r_3-r_2=2M/(r_3+r_2)\). The high-mode exponential correction tends to zero uniformly in \(m\) as \(b\to\infty\). As \(F=O(1/b)\), the deep difference stays bounded above and away from zero. This forces \(m<b\) eventually (otherwise the lower bound grows with \(\sqrt m\)), then \(m+1\asymp\sqrt b\). In (5) it follows that \(\sqrt{AM^2/b}=1+O(1/b)\), so \(d=O(1)\) and (6) applies. Its leading absolute value for any integer \(b\) is at least \(\tfrac12\operatorname{dist}(g(m),\mathbb Z)\). Taking \(m\) through a minimizing residue and rounding nearest attains the right side in (4) asymptotically.

## The two limiting laws

We have \(g(m)=A m(m+1)+A/4-1/(2A)\). For irrational \(A\), polynomial equidistribution of [Weyl](https://doi.org/10.1007/BF01475864) gives a uniform fractional part. Nearest-integer signed residual sends that distribution to uniform on \([-1/2,1/2]\); the map's isolated discontinuity has zero limiting mass, and two exact half-integer ties would already imply rational \(A\) by subtraction. Taking half and using (1) gives the asserted law and extremal limits.

For rational \(A=P/Q\), \(g(m+Q)-g(m)=P(2m+Q+1)\) is integral. An integer or half-integer \(g(m)\) would require \(P^2(2m+1)^2-2Q^2\) divisible by \(2PQ\). Modulo four the numerator is respectively \(1,3,2\) in the coprime parity cases (odd, even), (odd, odd), (even, odd). If \(PQ\) is even the numerator would have to be \(0\bmod4\); otherwise it at least must be even. Both contradictions prove the claim. Example: for \(\alpha=1/2,h=(\log2)/2\), \(\tanh h=1/3,w=1/\sqrt2,A=9/4\). The odd square is \(1\) or \(9\bmod16\), each on exactly two of the four residue classes. Formula (1) gives respectively \(49/288\) and \(-23/288\) with mass \(1/2\) each, and (4) gives an exact tail floor \(23w/288\) over all \(m,b\). At every fixed \(\alpha\), rational-\(A\) depths are a countable dense subset, since \(A(h)\) is a smooth strictly decreasing bijection onto \((9\alpha/[4(1+\alpha)],\infty)\).

At algebraic \(\alpha,h>0\), \(\tanh h\) is transcendental by Hermite–Lindemann and \(e^{2h}=(1+\tanh h)/(1-\tanh h)\), so \(A\) is in fact transcendental. The [fixed algebraic-depth companion (#17)](https://github.com/ConfusedDuckCollectibles/Research1/tree/79bb3017701a00fb83a68686a62bd39c9d4d16b3/discoveries/gravity-capillary-resonance) already excludes any exact unbalanced relation there, including these distinct-shell cubic relations (the two lower radii can coincide but their signed shell multiplicity remains unbalanced). That result is not needed for strict positivity in (2) or (3). #21 gives nonzero cubic coefficients of size \(\asymp_{\alpha,h}\sqrt b\) along these subsequences; no new Hamiltonian formula is asserted. The same arithmetic applies with any fixed nonzero lattice vector \(k\) as the first mode, by using the integral sublattice with basis \(k\) and its quarter-turn. With \(s=|k|\), replace \(\alpha,h\) by \(s^2\alpha,sh\) and rescale frequencies by \(\sqrt{s}\).

## Quantitative hits: (2)

Write \(e(x)=e^{2\pi i x}\). For convergents \(P/q\) to \(A\), use \(|A-P/q|<q^{-2}\), \(\gcd(P,q)=1\), and \(N=\lfloor q^{5/6}\rfloor\). Work on \(N\le m<2N\), a block of \(N<q\) consecutive residues. Define \(D_N\) as the supremum over circle arcs of the absolute difference between the empirical fraction in that arc and its length.

The complete Gauss sum \(G(a,j;q)=\sum_{x\bmod q}e((a x^2+jx)/q)\) obeys
\[
 |G|^2=q\sum_{\substack{d\bmod q\\2ad\equiv0\bmod q}} e((a d^2+jd)/q),\qquad
 |G|\le\sqrt{q\,\gcd(2a,q)}\le\sqrt{2q\,\gcd(a,q)} .
\]
Square the modulus and replace \(x\) by \(y+d\); summing over \(y\) is exactly the orthogonality condition. This also covers even \(q\). Complete the Fourier sum of \(e(k Pm(m+1)/q)\) on our block with the normalized finite Fourier transform of its indicator. The transform has absolute-value sum \(O(\log(2q))\): each nonzero coefficient is at most \(1/[2\min(\ell,q-\ell)]\) by the finite geometric series, and zero frequency is \(N/q\). Each completed sum has \(a=kP\) and a possibly shifted linear coefficient \(j\), so \(\gcd(a,q)=\gcd(k,q)\). Inserting this into the [Erdős–Turán discrepancy inequality](https://www.renyi.hu/~p_erdos/1948-03.pdf), Part II, Theorem III, with Fourier cutoff \(q\), gives
\[
  D_N\ll q^{-1}+
     \frac{\sqrt q\,\log(2q)}{N}
       \sum_{1\le k\le q}\frac{\sqrt{\gcd(k,q)}}{k}
  \ll_\gamma q^{-1/3+\gamma}\quad(\gamma>0).
\]
For the last step, overcount via \(d\mid q,d\mid k\); the sum is at most \(\sum_{d\mid q}d^{-1/2}(1+\log q)\le\tau(q)(1+\log q)\), and the elementary divisor bound \(\tau(q)\ll_\delta q^\delta\) absorbs the two logarithms. Discrepancy for a circle arc crossing zero follows equally from the usual interval inequality, by splitting it if needed.

The target for \(A m(m+1)\bmod1\) is \(c(A)=1/(2A)-A/4\bmod1\). To prove **nonzero**, take \(\rho=q^{-1/3+\gamma}\) with \(0<\gamma<1/3\), and ask the rational block to hit the oriented circle arc \(c(A)+[\rho,2\rho]\). Use a smaller exponent, e.g. \(\gamma/2\), in its discrepancy bound so that \(D_N=o(\rho)\). The rational-to-true phase shift on the block is \(O(N^2/q^2)=O(q^{-1/3})=o(\rho)\). Therefore the true fractional part of \(g(m)\) is in \([\rho/2,3\rho]\). For large \(q\) this rounds down, so (6) gives positive \(F\), since \(b^{-1}\asymp q^{-5/3}=o(\rho)\) inside the scaled phase. Also \(b\asymp m^2\asymp q^{5/3}\); hence \(F\ll b^{-6/5+3\gamma/5}\). Choose \(\gamma\) for the desired \(\varepsilon\). This uses infinitely many convergent denominators with \(N\to\infty\), so the selected pairs cannot stay finite.

## Metric hits: (3)

[Baker's general almost-everywhere discrepancy theorem for sequences of distinct positive integers](https://academic.oup.com/jlms/article-abstract/s2-24/1/34/862375) gives \(D_{[1,N]}=O(N^{-1/2}(\log N)^{3/2+\eta})\) for each \(\eta>0\), as stated explicitly in [Aistleitner–Larcher, §1, Eq. (1)](https://arxiv.org/pdf/1507.06472), and would apply to \(m(m+1)\) for \(m\ge1\). Differences of prefix counts yield the same rate for dyadic blocks. The following direct proof on blocks suffices, without invoking Baker. For each \(k\ge1\),
\[
 \left\|\sum_{m=N}^{2N-1} e(k A m(m+1))\right\|_{L^2([0,1],dA)}=\sqrt{N},
\]
by distinct integer frequencies. With Fourier cutoff \(N\), Erdős–Turán and the \(L^2\) triangle inequality give \(\|D_N\|_2\ll N^{-1/2}\log(2N)\). Chebyshev and Borel–Cantelli over \(N=2^j\) imply \(D_N=O(N^{-1/2}(\log N)^{3/2+\eta})\) a.e., since the exceptional measure is \(O(j^{-1-2\eta})\). Periodicity in \(A\) of this discrepancy extends the a.e. statement to \(\mathbb R\). Use countably many \(\eta\downarrow0\) to get all \(\eta>0\) on one full-measure set.

For a fixed such \(A\), use again the one-sided arc \(c(A)+[\rho_N,2\rho_N]\), with \(\rho_N=N^{-1/2}(\log N)^{3/2+\eta}\) and a smaller exponent \(\eta/2\) in \(D_N\). Uniform discrepancy over arcs is crucial: the target depends on \(A\), and moves with \(h\), but for each individual \(A\) the bound applies to *all* arcs. Thus \(\rho_N\le\{g(m)\}\le2\rho_N\) for some \(m\) per sufficiently large dyadic block. Since \(b\asymp N^2\), this arc is well separated from the \(O(b^{-1})\) scaled Taylor error, giving strict positivity and (3). For fixed \(\alpha\), \(h\mapsto A\) has nonzero derivative and a locally Lipschitz inverse on compact positive intervals, so null exceptional sets pull back to null depth sets.

## Reproduction and attribution

    python3 -B exact_replay.py | diff -u replay.txt -
    sha256sum -c SHA256SUMS

The replay independently encloses the original phase formula using 220-bit dyadic integer arithmetic, square roots by integer square root, and exact alternating rational exponential bounds. At the rational two-point control it validates 12 sign brackets and finds six instances where changing from unshifted to corrected rounding reverses the phase's sign. It exactly reproduces the 20,000 irrational-control bin counts for \(h=3/5,\alpha=1/2\) and \(100\le m<20100\), and checks 1,170 squared Gauss identities in integer cyclotomic quotients, including even moduli. These controls are checks of the derivation, not empirical evidence for the theorem.

The dispersion family, leading phase and interaction size are from #21. The addendum's contribution is the corrected center for this fixed-depth family, its rational/irrational law (including the full-family rational floor), and the quantitative phase approximations obtained by classical Weyl, Gauss and discrepancy machinery. Baker's primary publisher record was located; the complete primary text was not needed or obtained. [Aistleitner–Larcher's work on polynomial discrepancy](https://arxiv.org/abs/1507.00207) is close in the arithmetic input only; a lower bound for uniform discrepancy is not a lower bound on hitting this particular target. [Concus's finite-depth corrigendum](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Fluid%20Mechanics/1964%20Volume%2019/S0022112064000702.pdf) and [Abassi–Wilkening §3.1](https://arxiv.org/html/2401.00844v3) involve temporal harmonics rather than a noncollinear vector triad. [Ionescu–Pusateri](https://arxiv.org/abs/1807.02932) and [Feola–Montalto–Murgante](https://arxiv.org/abs/2604.20968) concern different infinite-depth estimates. No matching fixed-depth rounded lattice calculation was located in this bounded comparison; no priority claim.
