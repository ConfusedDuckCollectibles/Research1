# Fixed support and amplitude prevent arbitrarily close calibrated full-wave matching

**Claim.** Fix positive \(\tau_0,\omega_{\mathrm{ref}}\), a source position \(\beta\), a compact screen domain \(D\subset\mathbb R^2\), a finite uniform bound \(M\), a nontrivial compact observing frequency interval \(J\subset(0,\infty)\), and one vacuum point-mass comparator \(I_g\) at finite, nonzero impact and nonzero lag. Allow any real measurable \(N\) supported in \(D\) with \(\|N\|_\infty\le M\); in particular allow smooth nonnegative screens with arbitrarily large derivatives. Write \(I_N=|F_N|^2\) for the normalized cold-plasma Fresnel intensity below. Then
\[
\inf_{\substack{\operatorname{supp}N\subset D,\ \|N\|_\infty\le M\\ C\ge0}}
\;\sup_{\omega\in J}|I_N(\omega)-C I_g(\omega)|>0.
\]
This is a qualitative uniform gap; no numeric lower bound or instrument sensitivity is asserted. The constant \(C\) cannot vary with frequency. The geometry, time scale, band, comparator, and the two screen bounds are fixed, rather than allowed to diverge through a sequence.

**Proof.** Suppose there are \(N_n,C_n\) with the uniform mismatch on \(J\) tending to zero. For each positive integer \(m\), the sequence \(N_n^m\) is bounded by \(M^m\) in \(L^\infty(D)\). Since \(L^1(D)\) is separable (Lebesgue measure on a compact Euclidean domain), the weak-* closed bounded balls of its dual are metrizable and compact; extract successive weak-* subsequences and diagonalize. Write \(L_m\in L^\infty(D)\) for the simultaneous limits, still with \(\|L_m\|_\infty\le M^m\). No assertion that \(L_m=(L_1)^m\) is made.

Put \(Q(x)=|x-\beta|^2/2\) and \(a=\tau_0\omega_{\mathrm{ref}}^2\). At complex \(w\ne0\) the exact compact subtraction is
\[
 F_{N_n}(w)=1+\frac{w\tau_0}{2\pi i}
 \sum_{m\ge1}\frac{(-ia/w)^m}{m!}\int_D e^{iw\tau_0 Q(x)}N_n(x)^m\,dx.
\]
For any compact \(E\subset\mathbb C\setminus\{0\}\), \(\{e^{iw\tau_0Q(\cdot)}:w\in E\}\) is compact in \(L^1(D)\) by continuity on \(E\) and bounded \(Q|_D\). Uniform \(L^\infty\) bounds and a finite \(L^1\)-net consequently upgrade each weak-* integral convergence to uniform over \(w\in E\). The series tails are uniformly dominated there by a constant times the exponential series in \(aM/\min_{w\in E}|w|\). Hence \(F_{N_n}\) converges locally uniformly to the holomorphic \(F_*\) obtained by replacing \(N_n^m\) by \(L_m\). On positive real frequency, the summed \(m\ge2\) remainder is uniformly \(O(1/\omega)\) by the elementary exponential remainder (e.g. bounded by \((\tau_0|D|a^2 M^2/(4\pi\omega))e^{aM/\omega}\)). Therefore
\[
 F_*(\omega)=1-\frac{\tau_0 a}{2\pi}\int_D
 e^{i\omega\tau_0Q(x)}L_1(x)\,dx+O(1/\omega).
\]
The integral tends to zero by the Riemann–Lebesgue lemma: extend \(L_1\) by zero and integrate the angles around \(\beta\), then substitute \(s=|x-\beta|^2/2\) to get the Fourier transform of a compactly supported \(L^1\) function. Thus \(F_*\to1\) on the positive high-frequency axis.

The intensities converge uniformly on \(J\) to \(I_*=|F_*|^2\), real analytic on the positive axis (the holomorphic extension is \(F_*(w)\overline{F_*(\bar w)}\)). At any one point in \(J\) with \(I_g>0\), local boundedness of these intensities bounds \(C_n\). Such a point exists by the real analyticity and nontriviality of \(I_g\). Extract \(C_n\to C\ge0\). Then \(I_*=CI_g\) on \(J\), hence everywhere on the positive axis. For \(C=0\) this contradicts \(I_*\to1\). For \(C>0\) it contradicts the standard point-mass high-frequency two-ray asymptotic with persistent nonzero oscillation between distinct delays, again because \(I_*\to1\). This proves the gap.

**Related average, not needed for the gap.** For all such screens, the leading integral of \(N\) is the Fourier transform (with scale \(\tau_0\)) of their radial average in \(s\). These averages have uniformly bounded \(L^2\) norm because \(|N|\le M\) and the common radial \(s\)-support is bounded. Plancherel and the uniform \(O(1/\omega)\) remainder give a uniform bound for \(\int_1^\infty|F_N(\omega)-1|^2\,d\omega\). Therefore the mean squared voltage discrepancy on \([W,2W]\) tends to zero uniformly at \(O(W^{-1})\); the same holds for intensity discrepancy by uniform boundedness of \(F_N\) there. This is not uniform pointwise high-frequency decay, a quantitative lower bound on the original observing band, or a statement about arbitrary unknown frequency-dependent amplitude.
