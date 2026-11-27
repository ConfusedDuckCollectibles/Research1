# Quantitative stability of the outer-shell optimum

A quantitative bound follows directly from the proof’s radial trial field. Let $t=4\pi r^3/3$, $t_a=4\pi a^3/3$, $w=u/r^2$ and $W(t)=B_0^2 r^2 w(r)^2/2$. Then $W'(t)\ge \kappa=B_0^2w(0)^2/(4\pi R)>0$. The inside excess $M_{\rm in}=\int_{r<a}q$ and outside deficit $M_{\rm out}=\int_{r>a}(Q-q)$ have sum $\|q-q_*\|_1$ and satisfy $M_{\rm out}\ge M_{\rm in}$. Subtract $W(t_a)$ in the trial comparison; its budget term is nonnegative for the deficit. After integrating over angular directions each side has density bounded by $Q$ per unit $t$. The minimal first moment with given mass $M$ and this density cap is $M^2/(2Q)$, obtained by filling at cap starting at $t_a$. Therefore the weighted deficit is at least $\kappa(M_{\rm in}^2+M_{\rm out}^2)/(2Q)\ge\kappa\|q-q_*\|_1^2/(4Q)$. Restoring the normalization of the three-axis response gives

$$
\bar S(q_*)-\bar S(q)\;\ge\;\frac{\kappa\,\|q-q_*\|_1^2}{24\pi B_0^2 R^3 Q}.
$$

For $f=0$ the feasible set has only one profile. For $f=1$ the argument is one-sided with $t_a=0$. The constant can be very small; this is stability of the maximizing layout, not a robust general inversion of a measured magnetization.

## The quadratic exponent is sharp

For $0<f<1$, let $q_h$ be the radial layout obtained by moving a filled volume interval of length $h$ immediately outside $t_a$ to the equal-volume interval immediately inside. Then $\|q_h-q_*\|_1=2Qh$. Normalize the two radial solutions by $u(R)=1$. The cross-profile Wronskian (see [higher harmonics](QUADRUPOLE.md)) expresses the difference in boundary slopes as the signed stiffness integral with kernel $u_*u_h/(4\pi r^2)$ per unit $t$. These kernels are uniformly Lipschitz near the fixed $a>0$: the regular solutions satisfy $0<u\le1$, their slopes increase in $r$, and their boundary slopes are bounded by the all-$Q$ solution by the same Wronskian. The moved masses cancel the value of the kernel at $t_a$, leaving an integral of order $h^2$. The boundary-slope-to-response map has bounded derivative on this range. Thus $\bar S(q_*)-\bar S(q_h)=O(h^2)$, so a lower bound with any smaller exponent cannot hold uniformly near the optimum. This does not optimize the constant.
