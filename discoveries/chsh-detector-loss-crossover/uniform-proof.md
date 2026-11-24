# A uniform quantitative remainder at the CHSH detector-loss corner

This gives an alternative to §§3–5 of the companion proof. Keep its full Jordan angle square and the lower-compression positivity from §§1–2. Use the same notation: $s=2-a-b$, $t=a+s=2-b$, $h=s/t$, $L=2+a+b$, $A=4-a^2$, $B=4-b^2$, and $C=8-2s+ab$. The determinant and product are
$$
F_{u,v}(d)=P(d)+8uv\{ab\,d-sC+(A-2v)u+Bv\},\qquad
P(d)=d(d+4+2a)(d+4+2b)(d+2(a+b)).
$$
Write $P_0=P'(0)$ and $x_h=2/(3+\sqrt{9-8h})$, $g(h)=x_h^2(1-x_h)$. The conclusion is the explicit all-path estimate
$$
\boxed{\ \left|\frac{Q(a,b)-L}{s^2h}-g(h)\right|\le\frac t5
       \quad (0<t\le 1/100,\ \ 0<s\le t)\ .\ }
$$

### Exact optimization at the base score

First assume $a>0$. Set $J_0(u,v)=uv[sC-(A-2v)u-Bv]$. A block with positive gap has $F_{u,v}(0)=-8J_0(u,v)<0$, by the positive lower compression. On $[0,2]^2$, $A-2v\ge-a^2$ and $u\le2$. Hence any positive $J_0$ requires
$$
v<\frac{sC+2a^2}{B}<A/2.
$$
The strict second inequality holds **uniformly throughout this corner**, even when $a/s\to0$. Indeed $C=(2+a)(4-t)$ and $B=t(4-t)$, so
$$
\frac A2-\frac{sC+2a^2}{B}
=a\left[\frac{2+a-2a/(4-t)}t-1-\frac a2\right]>0
$$
since $t\le1/100$ and $a<t$. Consequently a positive value also requires $v<V:=sC/B=h(2+a)$. Complete the square in $u$, bounding its maximum by
$$
\frac{v(sC-Bv)^2}{4(A-2v)},\qquad 0\le v\le V .
$$
With $x=v/V$ and $\kappa=2V/A=2h/(2-a)$, the envelope is
$$
\frac{V(sC)^2}{2A}\cdot\frac{x(1-x)^2}{2(1-\kappa x)} .
$$
Here $\kappa<1$ because $h=1-a/t<1-a/2$. The derivative of the rightmost factor has sign $(1-x)(1-3x+2\kappa x^2)$ on $[0,1]$. The sole interior maximizer is $x_\kappa\in[1/3,1/2]$, and the value at it equals $g(\kappa)$ using $1-x_\kappa=2x_\kappa(1-\kappa x_\kappa)$. The maximizing $v_\star=Vx_\kappa$ is below 1, and
$$
u_\star=\frac{sC(1-x_\kappa)}{2A(1-\kappa x_\kappa)}\le\frac{sC}{A}<2.
$$
Thus the bound is actually attained by allowed angles; points with nonpositive $J_0$ do not compete. In particular
$$
J_{\max}=\frac{V(sC)^2}{2A}\,g(\kappa),\qquad
\frac{8J_{\max}}{P_0s^2h}=R\,g(\kappa),\qquad
R=\frac{4(2+a)C^2}{AP_0}=\frac{(2+a)(4-t)}{2(2-a)(2-s)}.
$$

### Uniformly compare the determinant deficit and the spectral gap

Each positive block gap $d_{u,v}$ solves $F_{u,v}(d_{u,v})=0$. Because
$$
F_{u,v}'(d)=P'(d)+8uvab\ge P_0\quad(d\ge0),
$$
the mean value theorem gives $d_{u,v}\le8J_{\max}/P_0$. In the maximizing block $(u_\star,v_\star)$, $F_{u_\star,v_\star}(0)<0$, so it has a positive gap and the same theorem gives an attainable lower bound with denominator $P_0$ plus the derivative change along the interval to its root.

For an explicit bound, note
$$
R-1=\frac{6a+2s-a^2-3as}{2(2-a)(2-s)}\in[0,0.76t],
\qquad 0\le \kappa-h=\frac{ha}{2-a}\le\frac{t}{1.99}.
$$
On $0\le h\le1$ the square root in $x_h$ is at least 1. Differentiating gives $0\le x_h'\le1/2$; because $x_h\in[1/3,1/2]$, $0\le g'(h)=x_h'(2x_h-3x_h^2)\le1/6$. Thus the profile is Lipschitz all the way to 1 and $g(h)\le1/8$. This gives
$$
0\le Rg(\kappa)-g(h)\le t(0.76/8+1/(6\cdot1.99))<0.18t .
$$
Also $d_{u,v}\le s^2h\,R g(\kappa)<t^2$ uniformly. In the product $P$, the three shifted constants are positive and at most 8. Consequently for $0\le d\le t^2$,
$$
0\le P'(d)-P_0\le384d+72d^2+4d^3\le385t^2 .
$$
With $8uvab\le64t$, the derivative change is less than $68t$; $P_0\ge4\cdot7.98\cdot3.98>127$. Therefore the maximizing block has normalized gap at least
$$
\frac{R g(\kappa)}{1+68t/127}\ \ge\ g(h)-\frac{68t}{8\cdot127}\ >\ g(h)-0.07t .
$$
Together these estimates imply the displayed error bound, with no division by a vanishing angle or by $1-h$.

At $a=0$, one may pass to the limit at fixed $t$ since $Q$ is continuous in the tilts (indeed, changing a coefficient changes any strategy's score by no more than its absolute change). Independently, the known singly tilted value gives normalized gap $1/[\sqrt{16-8t+2t^2}+4-t]$, between $1/8$ and $1/(8-2t)$, also satisfying the bound. In the independent-erasure/no-click-to-$+1$ model the recorded excess is $\eta_A\eta_B$ times the tilted gap; none of these arguments changes the local CHSH bound of 2.
