# Optional refinement: two real controls rather than three

Write $G_{\mathbf v}(r)=J_{\mathbf v}(\sqrt r)$ for squared horizontal input tangent norm $r\ge0$. The same free nominal controls make $G_{\mathbf v}$ concave: given a convex combination $r=\sum_i p_i r_i$, isometrically map the normalized incoming state and horizontal derivative to
\[
 \sum_i\sqrt{p_i}\,|i\rangle|\psi_i\rangle,\qquad
 i\sum_i\sqrt{p_i r_i}\,|i\rangle|\phi_i\rangle,
\]
with $\psi_i,\phi_i$ orthonormal for each $i$ and the labels orthogonal. The two target vectors have the original norms $1,\sqrt r$ and zero mutual inner product. At $r=0$ all $p_i r_i$ vanish. Reading the label with fixed probabilities to first order provides conditional tangent norms $\sqrt{r_i}$, and the conditional policies give $\sum_i p_i G_{\mathbf v}(r_i)$. A phase-dependent buffer was not prepared; the tangent norm was allocated from the incoming derivative. Arbitrary finite ancillary space suffices for this inequality; the independent tight recursion already returns the same value with constant quantum space.

For $0<v<1$, $d>0$, fixed $x\in(0,1)$ and real $z\in[0,\sqrt{x(1-x)}]$, the remaining admissible interval is
\[
 L=\frac{z^2}{x}\le y\le U=1-\frac{z^2}{1-x}.
\]
The two conditional squared tangent norms $(D_0,D_1)=(d_0^2,d_1^2)$ vary affinely with $y$, with $D_1=d^2(y/x-z^2/x^2)$ increasing and $D_0$ decreasing. Their fixed mean is
\[
 D=q_0D_0+q_1D_1
   =d^2+\frac{v x(1-x)+2vdz-(1-v)d^2z^2/x}{1-(1-v)x}.
\]
Both scores are zero at real $z$. For fixed $D$, the function
\[
 q_0 G_{\mathbf u}((D-q_1 D_1)/q_0)+q_1 G_{\mathbf u}(D_1)
\]
is concave on its natural nonnegative domain and takes its global maximum at $D_1=D$ by Jensen. Along the interval $[L,U]$ choose the closest feasible $D_1$ to $D$; $D_1=0$ at $y=L$, so only an upper clamp is needed:
\[
 y_{\rm opt}(x,z)=\min\!\left\{U,\ \frac{z^2}{x}+\frac{xD}{d^2}\right\}.
\]
This choice remains optimal if the continuation is not strictly concave. At $d=0$ the $y$ coordinate is irrelevant. At $x=0$ use $M=0$, and at $x=1$ take $z=0,y=1$ to preserve the incoming tangent in both branches. At $v=0$ discard the slot's action and use a dummy; at $v=1$ the jump is absent and $J_{(1,\mathbf u)}(d)=J_{\mathbf u}(d+1/2)$ by maximizing the surviving horizontal tangent and using monotonicity. Thus local optimization is over just $(x,z)$; no assertion about the global computational complexity of evaluating the recursion.

As all optimizing effect choices can be taken real, the recorded history can have identically zero first derivative of its probability at the nominal phase for every positive history. Its usefulness is the conditional quantum state and next control, not a first-order classical phase score.
