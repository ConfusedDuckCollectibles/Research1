# No 34-point spectral set in $\mathbb F_{17}^3$

This note gives an exact, computer-assisted proof that no subset of $\mathbb F_{17}^3$ with 34 elements is spectral. Equivalently, no two 34-element subsets of the paired three-dimensional vector spaces form a spectral pair.

Here a spectral pair means that the $34\times34$ matrix of characters $\exp(2\pi i\langle x,y\rangle/17)$, divided by $\sqrt{34}$, is unitary. Orthogonality first limits the configurations of difference directions. The proof combines exact finite checks in $\mathbb F_{17}$ with a sequence of incidence, moment, and curve arguments. The last obstruction is algebraic: it is a cubic forced to have three roots and to have zero leading coefficient.

- [Proof and dependency map](PROOF.md)
- [A. Affine geometry and the six-direction pencil](A-AFFINE.md)
- [B. Rich sections and polynomial moments](B-MOMENTS.md)
- [C. Cubic intersections and the septic reduction](C-CURVES.md)
- [D. The horizon schemes and final contradiction](D-HORIZONS.md)
- [Exact computation inventory and commands](REPRODUCE.md)

All finite checks use exact integer or finite-field arithmetic. The proof distinguishes complete gates from supplementary convention checks. The terminal obstruction does not require a finite enumeration.

A 34-element subset also cannot tile $\mathbb F_{17}^3$ by translations, since 34 does not divide $17^3$. This note settles that cardinality only. It does not address size 68 or claim to settle the full $\mathbb F_{17}^3$ Fuglede problem. This is an AI-assisted research note, not an external review or a claim of priority.
