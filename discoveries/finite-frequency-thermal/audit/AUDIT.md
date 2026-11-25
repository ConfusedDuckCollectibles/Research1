# Independent audit: finite-frequency ambiguity of reciprocal thermal-box feedback

**Result.** For finite connected reciprocal linear thermal boxes, strictly positive box capacities, no hidden losses/forcing, and one collocated surface heat input/temperature response, suppose surface capacity $C_0$, total hidden capacity $C>0$, and exact complex impedances at distinct positive angular frequencies $\omega_j$ are held fixed. Then every compatible pair has
\[
 |\lambda_+-\lambda_-|\le \frac{C}{2(\sum_j\omega_j^{-1/2})^2}.
\]
The constant is a sharp worst-case exact-data diameter. The pairwise bound is strict for finite strictly positive rates in the unrestricted class, but the same diameter is realized as a supremum at **one fixed admissible datum**. With a strictly positive closed spectral window $[a,b]$, the improved bound below is attained by finite exact-data pairs. Both model-class and frequency units matter: no climate estimation or noise robustness independent of spacing is claimed.

## 1. Positive rate measure and its realization

Let $K$ be the hidden principal block of the interbox graph Laplacian, $D$ the positive diagonal hidden capacity matrix, $k=K\mathbf 1$ the vector of direct surface conductances, and $b=D^{1/2}\mathbf 1$. The graph is connected, so $K$ is positive definite, and $S=D^{-1/2}KD^{-1/2}$ is positive definite. Eliminating hidden temperatures gives the driving-point admittance
\[
 Y(s)=Z(s)^{-1}=\lambda+sC_0+\mathbf 1^T k-k^T(K+sD)^{-1}k
 =\lambda+sC_0+b^T sS(S+sI)^{-1}b .
\]
Here $D^{-1/2}k=Sb$. In the spectral basis this is
\[
 Y(s)=\lambda+sC_0+\int\frac{sx}{s+x}\,d\mu(x),\quad
 \mu\ge0,\quad \operatorname{supp}\mu\subset(0,\infty),\quad \mu((0,\infty))=\|b\|^2=C,
\]
where the mass at a rate is the squared projection of $b$ onto its rate eigenspace; inactive modes carry no mass. Conversely an atom of mass $c>0$ at $x>0$ is an independent hidden box with capacity $c$ and conductance $cx$ to the surface. Any finite positive measure on strictly positive finite rates is thus realized by a connected positive thermal star, with the same $C_0$. At $s=i\omega$, $Y$ cannot vanish: its real part is strictly positive for $C>0,\lambda\ge0$ (and its imaginary part is also positive). Hence exact equality of $Z$ and of $Y$ is interchangeable.

With $g_j(x)=(x^2+\omega_j^2)^{-1}$ and $f_j(x)=\omega_j^2x/(x^2+\omega_j^2)$,
\[
 \Re Y(i\omega_j)=\lambda+\int f_j\,d\mu,\qquad
 \Im Y(i\omega_j)=\omega_j(C_0+C)-\omega_j^3\int g_j\,d\mu.
\]
Write $P(x)=\prod_j(x^2+\omega_j^2)$ and $\mathcal R=\{p/P:\deg p\le n=2m,\ p\text{ real}\}$. Partial fractions at the distinct simple poles show $1,g_j,f_j$ are a basis. Of these only $f_j$ have nonzero derivative at zero, and all $f'_j(0)=1$. Exact observation equality for $\nu=\mu_+-\mu_-$ gives
\[
 \int R\,d\nu=-(\lambda_+-\lambda_-)R'(0)\quad(R\in\mathcal R).
\]
Conversely this equality implies every observation match. Thus the missing direction is precisely the derivative evaluated at zero; $\|\nu\|_{\rm TV}\le2C$ gives upper bound $2C/D_*$, where $D_*$ is the norm of that derivative on $\mathcal R$ equipped with the support-window sup norm.

## 2. Exact norm on the entire positive-rate axis

Write $x=y^2$ and let $B(y)$ be the product of $(y-\zeta)/(y-\overline\zeta)$ over the upper-half-plane roots of $P(y^2)$. For each $\omega$ the pair is $\sqrt\omega e^{i\pi/4},\,\sqrt\omega e^{3i\pi/4}$. The even rational function $(B+B^{-1})/2$ is $R_*(y^2)$ for $R_*\in\mathcal R$; common denominator is $P(y^2)$ up to a constant. On $y$ real, $B=e^{i\Phi}$ with $\Phi(0)=0$ and $\Phi'>0$. Its derivative is a sum of $2\Im\zeta/|y-\zeta|^2$, even in $y$. The full-line increase is $4m\pi$, so $\Phi$ on $[0,\infty]$ rises by $2m\pi$. The $n+1$ points $x_i$ at which phase is $i\pi$ alternate with values $(-1)^i$; $x_0=0,x_n=\infty$. In particular $\|R_*\|=1$, $\Phi'(0)=2\sqrt2\sum_j\omega_j^{-1/2}$, and $R_*'(0)=-\Phi'(0)^2/2$.

For completeness the derivative's norm equality follows directly from rational Lagrange cardinal functions $L_0,\ldots,L_n$ at these points, including the value at infinity as a leading-coefficient condition. Their derivatives at zero have signs $(-1)^{i+1}$. For $i=0$, $L_0=P(0)\prod_{k=1}^{n-1}(1-x/x_k)/P(x)$ has strictly negative derivative since $P'(0)=0$. For $0<i<n$, the numerator includes the factor $x$ and terms $(x-x_k)$ for $k\ne i,\,1\le k<n$, giving sign $(-1)^{i-1}$ after its cardinal normalization. The infinity numerator is a positive multiple of $x\prod_{k=1}^{n-1}(x-x_k)$, again with the correct sign. Consequently
\[
 |R'(0)|\le \sum_i|L_i'(0)|\,\|R\|,\qquad
 \sum_i |L_i'(0)|=|R_*'(0)|=4(\sum_j\omega_j^{-1/2})^2 .
\]
This bypasses any unverified exact Markov constant. If a finite pair at strictly positive finite rates with nonzero $\lambda$ difference attained the bound, equality of the total-variation inequality forces support on the alternating sites with the alternating signs. The unique derivative weights at these sites are nonzero at 0 and infinity, which are forbidden physical rates. Thus no such pair attains equality.

## 3. Sharpness by exact stars, including a fixed-data version

Derivatives of the cardinal functions represent $R'(0)$ exactly as a signed weighted sum of $R(x_i)$. Their total weight is zero (test $R=1$). Move the endpoints to positive finite rates and keep the interior nodes or perturb them consistently; weights vary continuously because the interpolation matrix in the basis $x^k/P(x)$ extends to infinity and is invertible at distinct limiting sites. The two Jordan halves normalized by $2C/\sum_i |L_i'(0)|$ each have mass $C$ and give networks by the star realization. The $\lambda$ difference is exactly this normalization factor up to sign, and tends to $2C/D_*$. A common $\lambda$ shift makes both strictly positive and does not affect compatibility.

A fixed common datum stronger than these varying-data pairs is also possible. Freeze the odd-node measure (all $m$ rates strictly interior) from the exact compactified extremizer. The opposite measure has both endpoints and $m-1$ interior nodes. Holding mass, all $g_j$ moments and all $f_j-f_1$ moments fixed amounts to $2m$ real constraints. Treat the $m+1$ even-node weights and the $m-1$ interior locations as unknowns; treat endpoint perturbations as inputs. The Jacobian is invertible: a row dependence would give $R\in\ker R'(0)$ vanishing at both endpoints and vanishing doubly at every interior even node (the weights are strictly positive). On $[0,\infty]$ its numerator is necessarily a scalar multiple of $x\prod_{\mathrm{even,\ interior}}(x-x_i)^2$ (vanishing at infinity lowers the maximum degree by one), whose derivative at zero is nonzero unless the scalar vanishes. For a finite ceiling $b$ with $a=0$ the same polynomial has additional factor $(x-b)$, with the same contradiction. The implicit function theorem, using reciprocal rate as local coordinate near infinity, keeps capacities positive and rates ordered while moving forbidden endpoints inward. Adjust one $\lambda$ to match the common real offset, first taking a fixed common shift so all nearby $\lambda$ remain strictly positive. This holds the complete datum of the odd-node physical star fixed. In unrestricted rates both endpoints need movement; when $a=0,b<\infty$, $b$ is already physical. Thus sharpness need not involve drifting data.

## 4. Finite spectral windows

Assume every active rate is in $[a,b]$, $0\le a<b<\infty$. Use $x=(a+b y^2)/(1+y^2)$ and replace the original upper poles by $\zeta_j=\alpha_j+i\beta_j$ and $-\overline\zeta_j$, with
\[
 \zeta_j^2=\frac{i\omega_j-a}{b-i\omega_j},\quad \alpha_j,\beta_j>0.
\]
The same product $B$ and its symmetric part again produce an element of $\mathcal R$: the denominator in $y$ is $(1+y^2)^nP((a+b y^2)/(1+y^2))$ up to constant, and an even numerator of degree $\le2n$ becomes a numerator in $x$ of degree $\le n$. Phase alternation holds unchanged, giving $a=x_0<\cdots<x_n=b$. With $a>0$, the point of derivative evaluation is outside the permitted rate interval, but the derivative-cardinal signs still alternate: $L_i(0)$ has sign $(-1)^i$ and $L'_i(0)/L_i(0)=-\sum_{k\ne i} x_k^{-1}$ because $P'(0)=0$. For $a=0$, endpoint cardinal signs are the same, by the factor-$x$ calculation above (no infinity site now). Thus the exact width is $2C/|R_*'(0)|$; when the interval is closed and $a>0$, each alternating site is allowed and the sharp pair truly exists.

For $a=0$, $x=b y^2+O(y^4)$ and
\[
 |\Im(1/\zeta_j)|^2=\tfrac12\big(1+\sqrt{1+(b/\omega_j)^2}\big),\quad
 |R_*'(0)|=\Phi'(0)^2/(2b),\quad
 \Phi'(0)=4\sum_j|\Im(1/\zeta_j)|.
\]
This proves the width $Cb/[2(\sum_j\sqrt{1+\sqrt{1+(b/\omega_j)^2}})^2]$, only as a supremum when zero rates remain forbidden.

For $a>0$, write $\eta=\sqrt{a/b}$ and
\[
 h(v)=\sum_j\log\frac{\alpha_j^2+(v+\beta_j)^2}
                            {\alpha_j^2+(v-\beta_j)^2}.
\]
Then $x=0$ corresponds to $y=i\eta$; $B(iv)=e^{-h(v)}$ is positive, and $R_*(x(iv))=\cosh h(v)$. Differentiating $x(iv)$ gives
\[
 |R_*'(0)|=\frac{(1-\eta^2)^2}{2\eta(b-a)}\sinh h(\eta)\,h'(\eta)>0.
\]
Each summand of $h'$ is positive at $\eta$, since
\[
 h'_j(v)=\frac{4\beta_j(|\zeta_j|^2-v^2)}
 {(\alpha_j^2+(v+\beta_j)^2)(\alpha_j^2+(v-\beta_j)^2)},\qquad
 |\zeta_j|^2=\sqrt{(a^2+\omega_j^2)/(b^2+\omega_j^2)}>a/b .
\]
A degenerate one-point allowed positive window has no ambiguity since its measure is determined by $C$; no division by $b-a$ there.

## 5. Affine estimate and what errors it does not remove

If $R_*=A_0+\sum A_jg_j+\sum B_jf_j$, put $D_s=\sum B_j=R_*'(0)$ and
\[
 \widehat\lambda=\frac{A_0C+\sum_j A_j[(C_0+C)/\omega_j^2-\Im Y_j/\omega_j^3]
                        +\sum_jB_j\Re Y_j}{D_s}.
\]
On exact data $\widehat\lambda-\lambda=\int R_*d\mu/D_s$, so $C/|D_s|$ is a worst-case absolute error bound, optimal because the exact data can have twice this diameter. For real/imaginary additive admittance errors $\epsilon^R_j,\epsilon^I_j$, the estimate changes by $\sum_j(B_j\epsilon^R_j-A_j\epsilon^I_j/\omega_j^3)/D_s$. Its error bound is not independent of frequency separation. If measurement error is instead specified on $Z$, inversion adds its own conditioning; for $|\delta Z|<|Z|$, $|\delta Y|\le|\delta Z|/[|Z|(|Z|-|\delta Z|)]$. Exact distinct frequencies can approach one another while retaining separate moment constraints; duplicated frequencies and finite-error nearby measurements are not equivalent to exact separated measurements. The result does not convert into a uniform bound on $1/\lambda$ without an independent strictly positive lower bound on $\lambda$.

## 6. Primary-source comparators and checking scope

* [Geoffroy et al., two-layer EBM, Part I](https://doi.org/10.1175/JCLI-D-12-00195.1) gives the one-hidden-box climate equations and calibrates against simulated step forcings; [Fredriksen–Rypdal](https://arxiv.org/abs/1701.07355) treats multibox linear response via weighted exponential modes; [Lovejoy, §3.3](https://esd.copernicus.org/articles/12/469/2021/#section3) explicitly uses complex climate thermal impedance/admittance and phase. Neither the EBM nor that measurement convention is new here. An EBM with uptake efficacy/pattern effects not modeled by a reciprocal graph is not automatically in scope.
* [Ziegeler–Nolte–Schweizer, thermal structure-function network identification](https://doi.org/10.3390/en14217068), especially §2, gives the standard Foster/Cauer positive thermal spectral/network picture; it reconstructs/fits from transient response and does not give a radiative-feedback finite-frequency diameter. Here the use of *admittance* and the physical hidden capacity mass follow by the explicit elimination above.
* [Borwein–Erdélyi–Zhang](https://doi.org/10.1112/jlms/50.3.501) and [Borwein–Erdélyi, especially Theorem 5](https://people.tamu.edu/~terdelyi/papers-online/mathematika.pdf) supply the standard rational Chebyshev/Blaschke alternation and sharp Bernstein–Szegő framework. [Min, Thm. 2.1, Eq. 2.3](https://doi.org/10.4153/CJM-1998-008-3) treats global derivative bounds and gives endpoint values for the normalized rational Chebyshev extremal. The present endpoint/exterior derivative minimax application should not claim invention of the extremal or the general approximation technique.
* [Grabovsky, Stieltjes reconstruction from approximate values](https://sites.temple.edu/yury/files/2023/06/Stieltjes.pdf), especially §§2.1–2.3, and [Dyukarev](https://www.mathnet.ru/links/6f6120257449c4ade0d44040bb752b04/sm1276_eng.pdf) cover Nevanlinna–Pick feasibility/recovery or all solution families for general Stieltjes data. Those are substantially broader interpolation theories, not the specific exact one-missing-derivative positive thermal mass diameter with an explicit sharp frequency law checked here. No claim of exhaustive priority.

The replay uses exact fractions for independent $m=1,2$ witnesses, a one-frequency finite ceiling and a strictly positive attained window, and a direct nonstar node-elimination control. It is a set of finite controls, not the continuum proof.
