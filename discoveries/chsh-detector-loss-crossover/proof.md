# Sharp all-path crossover of the CHSH margin at a detector-loss corner

**Scope:** proved for the explicitly stated independent-erasure model, with an independent determinant derivation and 17 finite rational controls. This is a narrow asymptotic result, not a new Bell threshold, self-test, or claim that no one has previously considered the crossover. The exact quantum maximum was already characterized algebraically.

## Statement and operational interpretation

Consider CHSH with two binary settings each. Use every preidentified trial; a lost click is recorded as +1, regardless of local setting. The recorded local bound is exactly 2, with no fair-sampling or postselection assumption. To predict a quantum strategy's score, suppose Alice's and Bob's click losses are independent of each other, the inputs, underlying true outcomes, and shared state, with fixed efficiencies eta_A and eta_B. This is a **model for the attainable quantum margin**, not a hypothesis imposed on hidden-variable competitors or a way to rescue a post-selected score. Let

- a = 2(1 - eta_B)/eta_B, b = 2(1 - eta_A)/eta_A;
- s = 2 - a - b > 0, t = a + s = 2 - b;
- h = s/t in (0,1], x_h = 2/(3 + sqrt(9 - 8h)), and g(h) = x_h²(1 - x_h).

Write Q(a,b) for the quantum maximum of `CHSH + a<A0> + b<B0>`. As t tends to zero,

```
Q(a,b) - (2 + a + b) = (s³/t) [g(s/t) + o(1)] .
```

The [quantitative companion proof](uniform-proof.md) bounds the normalized error by t/5 for t at most 0.01. The remainder here is uniform over **every path** 0 < s <= t to the corner (a,b) = (0,2). For h bounded a positive distance below 1, the bracket remainder is O(t), uniformly even as h tends to zero. Exchange the party labels for the other corner.

The actual optimal all-trial excess above 2 is eta_A eta_B times this expression. This follows by conditioning on which detectors fire: joint clicks carry the underlying CHSH; a single missed detector produces a contribution proportional to the opposite party's first marginal; double misses give CHSH=2. Explicitly,

```
S_recorded,max = eta_A eta_B Q + 2(1-eta_A)(1-eta_B),
S_recorded,max - 2 = eta_A eta_B [Q - (2+a+b)] .
```

For any actual full-data distribution, a valid statistically significant certified CHSH score above 2 certifies nonlocality and hence entanglement without invoking this loss model. The asymptotic margin is a bound on tolerable *mean-score reduction* in the model, not a finite-trial p-value, a guarantee of extractable entanglement, an experimental prescription, or a universal result for correlated loss, other Bell expressions, or treating no-click as a third outcome.

Two noncommuting limiting operations cause the practical difference:

- One detector perfect (a=0, h=1): the tilted gap is s²/8 + O(s³), agreeing with the exact single-tilt formula. In raw CHSH, the leading coefficient tends to 1/16.
- The better detector nonideal and near-perfect with s << a << 1: h tends to zero, so the tilted gap is asymptotic to 2s³/(27a). The fixed-interior cubic law cannot be applied uniformly up to perfect efficiency.
- For a/s → kappa finite, the coefficient of s² is `g(1/(1+kappa))/(1+kappa)`. For equal scales (kappa=1), this coefficient is 0.0450849718747371, only 36.1% of the perfect-detector 0.125 coefficient despite both sequences tending to the same efficiency pair.

## Proof

### 1. Full-angle reduction and determinant identity

For binary quantum measurements, dilate observables to projective reflections on one enlarged Hilbert space per party (a shared embedding for both), then apply Jordan's lemma to each pair. In nontrivial two-dimensional blocks local rotations bring A0=B0=Z, `A1=(1-u)Z+xX`, `B1=(1-v)Z+yX` with u,v in [0,2], x²=u(2-u), y²=v(2-v). Signs of x and y can be changed by rotations preserving Z. Commuting/single-dimensional blocks are local. The Bell expression is affine in the behavior, so block mixtures cannot exceed the largest block value. It suffices to prove the upper bound on the **entire** angle square, not to assume the known maximizing angle solution or the claimed positive-cosine restriction. The trivial all-positive strategy realizes the local bound L=2+a+b; all 16 vertices satisfy it.

Since CHSH = 2A0B0 - (A1-A0)(B1-B0), write M(d) = (L+d)I - BellOperator = D(d)+K, with

```
D(d) = diag(d, d+4+2b, d+4+2a, d+2(a+b)),
K = [[-u,x],[x,u]] tensor [[-v,y],[y,v]] .
```

An exact polynomial expansion gives
```
F(d) = det M(d) = P(d) + 8uv[ab d - s C + (A-2v)u + Bv],
P(d) = d(d+4+2a)(d+4+2b)(d+2(a+b)),
A = 4-a², B = 4-b² = 4t-t², C = 8-2s+ab.
```
This identity is checked symbolically after reduction by x² and y². It doesn't use the published sextic or its corrections.

### 2. PSD check, including the opposite-angle corner

Prove the lower 3x3 compression Q of M(0) (basis |01>, |10>, |11>) is positive definite whenever 0<a<t<=1/100. The diagonal D on that subspace is at least (4-2s)I and the norm of K is 2sqrt(uv). Thus Q is strictly positive if uv < (2-s)². This leaves only the *other* angle corner, which a local series around u=0 would miss.

Otherwise let e=2-u, f=2-v. Both lie in [0,.02], since u,v<=2 and uv>=(2-.01)². The 2x2 subblock R on |01>,|11> has diagonals at least 3.98 and 7.94 and off-diagonal magnitude at most .4. It is positive definite. Its inverse has diagonals ≤.30 and .13, with off-diagonal magnitude ≤.02. The coupling from |10> into this block is (xy,uy) up to signs. We have (xy)²<=.08 f, (uy)²<=8 f, giving `w^T R^(-1) w <= 1.1f` by those rational inverse bounds. The Schur complement at |10> is at least
```
4+2a-uv-1.1f = 2a+2e+.9f-ef > 0.
```
The inequalities behind the inverse bounds are checked as rationals. Hence M(d) has at most one nonpositive eigenvalue for d>=0. Positive determinant at a trial d certifies all quantum block values are below L+d, and negative determinant is a strict witness. The a=0 endpoint will be handled exactly rather than extending a strict inequality there.

### 3. Localization uniformly away from h=1

Fix any epsilon in (0,1). Work with 0<h<=1-epsilon, t small. For an angle pair with a positive gap, F(0)<0 and u,v>0, so
```
(A-2v)u + B v < s C.
```
On the full angle square, A-2v>=-a² and u<=2, so initially
```
v <= (s C+2a²)/B = 2h[1+O(t)] + O(t).
```
Thus v stays a positive distance below 2 and A-2v is bounded below by a positive constant depending on epsilon. Applying the inequality again removes the additive O(t) leakage: `v <= s C/B = 2h[1+O(t)]` and `u=O_epsilon(s)`. This last refinement is essential when h tends to zero much faster than t. It is a genuine global localization, not just a bound for a presumed stationary point.

### 4. Scaling profile and global optimum

Use u=sU and v=2h x. For angles with positive gap, U is bounded and x<=1+O(t). If d_uv>0 is the gap for this pair, then d_uv≤||K||→0 uniformly from the block norm bound. The mean value theorem gives `d_uv=-F(0)/F'(xi)`, where xi is between 0 and d_uv. By the determinant identity `F'(xi)=128+O(t+xi)`; near zero it is bounded away from zero. Since F(0)=O_epsilon(s² h), this gives d_uv=O_epsilon(s²h) and improves the derivative remainder to O_epsilon(t). Substitution now yields
```
d_uv/(s² h) = (U x/8)[8(1-x)-4(1-h x)U] + O_epsilon(t).
```
For x>1 in the small allowable extension, the limiting leading expression is nonpositive. For x in [0,1], maximize first over U>=0: choose U=(1-x)/(1-h x), giving `x(1-x)²/[2(1-hx)]`. Its derivative has numerator `(1-x)(1-3x+2h x²)`. The unique interior maximizing critical point is x_h, including the continuous h=0 limit. Substitution, using `1-x_h=2x_h(1-hx_h)`, reduces the value to g(h). Conversely, choose exactly these U and x in the valid finite block. Its determinant at d=0 is strictly negative uniformly after scaling, and the same expansion achieves g(h)+O_epsilon(t). This proves both global upper and achievable lower bounds. No search remains in the leading coefficient.

### 5. Uniform extension to h→1

At a=0 the known singly tilted bound is `sqrt(8+2(2-s)²)`, so the normalized gap at h=1 tends to 1/8=g(1). At fixed s, Q(a,2-s-a) is convex in a (supremum of affine functions) and symmetric about a=(2-s)/2 (party exchange), so it is nonincreasing as a moves away from zero toward the midpoint. Fix h0<1 near 1. For h>=h0, at fixed s the gap is squeezed by the gaps at h0 and 1. Dividing by s²h yields lower bound h0/h times [g(h0)+O(s/h0)] and upper bound 1/h times [1/8+O(s)]. As h0 tends to 1 both converge uniformly to g(1); g(h) is continuous there. This establishes the claimed all-path o(1), without a false assertion that the away-from-endpoint error constant is uniform at h=1.

## Exact finite controls

`verify.py` checks the determinant, the local vertices, the loss transfer, the exact one-perfect limit, and rational Schur bounds. It then brackets the true optimum at 17 rational (a,s) points, with no numerical SDP or floating eigensolver. For trial d let J=sC-abd. If a gain term is possible it forces v<=(sC+2a²)/B<A/2 on these grids. Maximize the term u v[J-(A-2v)u-Bv] over u first. The upper envelope is `v(J-Bv)²/[4(A-2v)]` on v∈[0,J/B]; its unique interior maximum solves `4Bv²-3ABv+AJ=0` at the smaller root. A rational square-root interval upper-bounds the global envelope and a rational (u,v) inside certifies the lower. Bisection returns intervals with relative width <1e-17. Rational endpoints and lower witness parameters are preserved in `checks.json`. These exact finite comparisons are audits of the proof, not the proof of the uniform limit.

| Regime | t | Enclosed gap divided by s²h (rounded) | Limiting g(h) |
|---|---:|---:|---:|
| h=0.9 | 1e-7 | 0.114450742496182 | 0.114450738657388 |
| h=0.5 | 1e-7 | 0.090169948790128 | 0.090169943749474 |
| h=1/3 | 1e-7 | 0.083785180923414 | 0.083785175652046 |
| h=0.001 | 1e-7 | 0.074098781966667 | 0.074098776411729 |
| s=1e-10, a=1e-4, h≈1e-6 | 0.0001000001 | 0.074079654511107 | 0.074074098765418 |
| s=1e-5, a=1e-12, h≈1 | 0.000010000001 | 0.125000300000715 | 0.124999987500004 |

