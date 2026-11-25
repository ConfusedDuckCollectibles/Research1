# Size 52: the four-full-layer branch is impossible

Step 2 of the size-52 proof, 25 November 2026. If a spectral pair in \(\mathbb F_{13}^3\) has size 52, neither member contains a full 13-point affine plane. Combined with the previously verified size-52 reductions, **both plane caps are at most 10**. This step removes the complete-layer branch; later steps handle all other cases.

## Inputs and prior art

This is a continuation of the independently audited size-52 theorem. Use the [initial reductions](../01-reductions/PROOF.md) for proofs and runnable verification of the following input facts:

1. Both members have affine line cap 4 and plane cap 13; a determined difference direction of a member is a uniform normal of its partner, where each parallel affine plane has four points.
2. Every 13-point planar set of line cap four determines at least 11 projective directions (exact planar check). Every such 11-uniform pencil gives an injective projection along its perpendicular direction, with image four full parallel lines. The latter is the Radon profile argument; the four-special-direction endpoint is already in Qihang Wang, [arXiv:2609.20495v1](https://arxiv.org/html/2609.20495v1), and his constant/affine cases use the framework of [Kiss and Somlai](https://doi.org/10.1007/s10623-024-01404-y). The layered proof below only needs the 11-uniform case, not the 10-uniform endpoint.
3. All mixed moments of total degree at most five vanish for a spectral pair of size 52. The audited argument adapts the technique of [Ferguson and Iosevich, §11](https://arxiv.org/html/2610.30712v1#S11) and closes the prime-13 cap-13 exceptional configuration with a line-factor/ray obstruction.

Ferguson and Iosevich exclude a different four-layer situation at prime 11 by analyzing fourth and fifth moments of individually normalized layer graphs. Here the finite check couples two graph functions using the small-image condition below and needs only their first moment and one mixed moment. This note does not claim publication priority.

## 1. Dual coordinates for the layered case

Suppose one member has at least one full 13-point plane. Its 11 determined tangent directions force the partner into four full parallel layers. The partner now has a full plane, so the first member is also in four full layers. Denote one layer normal of the first member by \(n_E\), in the dual space \(V^*\). Its full planar slice gives the uniform pencil of normals of \(A\) in \(\ker(n_E)\subset V\). The associated projection of \(A\) along \(n_E\) is injective and the layer normal \(n_A\) it gives is in \(\ker(n_E)\). The full layer of \(A\) having normal \(n_A\) similarly makes projection of \(E\) along \(n_A\) injective. It is not necessary that the latter application reproduce \(n_E\) as layer normal. Choose a basis with \(e_2=n_A\), \(e_1\in\ker(n_E)\) independent of \(e_2\), and \(n_E(e_3)=1\). In these coordinates and their dual one has
\[
 E=\{(x,f_s(x),s):x\in\mathbb F_{13},s\in S\},\quad
 A=\{(u,t,g_t(u)):u\in\mathbb F_{13},t\in T\},
 \qquad |S|=|T|=4.
\]
All graph slices are nonaffine because their ambient line caps are four. Multiplying a monomial in \(x,f_s(x)\) of degree at most two by the cubic Lagrange indicator of layer \(s\) isolates its slice moment, which vanishes by input 3. We only use
\[
 \sum_x f_s(x)=0,\qquad \sum_x x f_s(x)=0 .
\]

## 2. At least eleven small-image shifts

Let \(f=f_s\), \(g=f_{s'}\), \(d=s-s'\ne0\), and use \(\zeta=\exp(2\pi i/13)\). For a shift \(h\in\mathbb F_{13}\), define \(D_h(x)=f(x)-g(x-h)\). Comparing the two \(E\)-rows at first coordinates \(x,x-h\) gives
\[
 \sum_{t\in T}\zeta^{tv}\Phi_t(h)=0 \quad
  (v\in\operatorname{im}D_h),\qquad
 \Phi_t(h)=\sum_u\zeta^{hu+d\,g_t(u)}.
\]
If the image contains at least four values, choose four and invert the indicated minor of the \(13\times13\) Fourier matrix, so every \(\Phi_t(h)=0\). Nonvanishing of all such minors is [Chebotarev's theorem](https://arxiv.org/abs/math/0312398). The optional check_minors.py also verifies the exact minors needed here without numerical evaluation: there are seven affine orbits of four-subsets, and the determinant is invariant up to a root-of-unity factor and Galois automorphism under independent affine changes of the two index sets. For each of 49 representative pairs the script accumulates signed counts of the 24 exponent residues; they are not all equal, and the cyclotomic polynomial is \(1+X+\dots+X^{12}\).

For any fixed \(t\), the graph of \(g_t\) determines at least 11 directions. Each such direction is the kernel of a functional \((u,y)\mapsto hu+d\,y\) for exactly one \(h\): the only projective normal not of this form has kernel vertical, a missing direction of any graph. At those at least 11 shifts a level is repeated, so not all level counts are one, and \(\Phi_t(h)\ne0\) by cyclotomic irreducibility. Thus at least eleven of the thirteen shifts have image size at most three.

## 3. Exact two-array enumeration

**Pair gate.** Let \(f,g:\mathbb F_{13}\to\mathbb F_{13}\) be nonaffine, with
\[
 \sum f=\sum g=\sum xf=\sum xg=0.
\]
If \(|\operatorname{im}(f(x)-g(x-h))|\le3\) for at least eleven shifts, it is impossible for both graphs to have line cap at most eight. In particular cap four is impossible.

To reduce completely to finite arrays, first note that every \(D_h\) has first moment zero. A function of 13 arguments with exactly two distinct values has nonzero first moment: it takes one value \(n\) times for \(1\le n\le12\), and the sum is \(n\) times the nonzero difference of the two values. At most one shift has constant \(D_h\): two such shifts would give a constant finite difference for \(g\) along a nonzero step, and hence an affine \(g\) and \(f\) because that step generates the additive group. Thus at least ten shifts have exactly three values. With at most three other positions in the cycle of length 13 there are four consecutive three-value shifts. Translate the \(x\)-argument of one of the graphs to make these \(h=0,1,2,3\). Translating or shearing does not affect the graph line caps; first and mixed moments remain zero under the graph transformations here since the sums of \(1,x,x^2\) over the field vanish and \(\sum f=\sum g=0\).

Write \(a=D_0\) and \(b=D_1\) for the two length-13 arrays. Both \(\sum a=\sum b=\sum xa=\sum xb=0\). There is also the matching bucket condition \(\sum x^2 a=\sum x^2 b\), because
\[
 \sum_x x^2(b(x)-a(x))
 =\sum_x x^2(g(x)-g(x-1))=-\sum_x(2x+1)g(x)=0.
\]
Their exact recovery is \(g(0)=0\), \(g(x)=g(x-1)+b(x)-a(x)\) for \(1\le x\le12\), and \(f(x)=a(x)+g(x)\). The equation at \(x=0\) follows from \(\sum(b-a)=0\); changing the shared vertical constant does not affect any condition. In particular all pairs of original graphs are covered; there is no assumption of an independent choice of the second- or third-shift image supports.

The following normalizations are valid in order.

* Subtract different vertical constants from \(f,g\) and scale both heights by a common nonzero factor. This allows restricting the first image support to any one representative in each affine orbit of three-subsets (3 possibilities). Every other support changes the same way; it is not separately scaled.
* A common shear \(f(x)\mapsto f(x)+k x,\ g(x)\mapsto g(x)+k x\) leaves \(D_0\) unchanged and adds \(k h\) to \(D_h\). Hence the second support can be put in its translation orbit representative (22 possibilities), independently after the first step. All image cardinalities and cap values are preserved.
* A common \(x\)-translation rotates both arrays, so take the lexicographically least cyclic rotation of \(D_0\) without any condition on the rotation of \(D_1\). This does not alter either support, and preserves the first/mixed moment conditions and matching buckets.

The first implementation generates all choices of the first eleven entries for each normalized support and determines the last two entries uniquely by the two moment equations \(\sum a=\sum xa=0\). Both supports must be used exactly, and pairs are formed only within the same second-weighted-moment bucket. Exact D2 and D3 images are checked first using \(D_h(x)=a(x)+\sum_{r=0}^{h-1}(b(x-r)-a(x-r))\). The second implementation independently generates the three-value arrays from all disjoint value masks, then recovers the actual graphs before checking direct differences at all shifts. In each implementation the cap is computed by the histograms of \(f(x)-k x\) and \(g(x)-k x\) for all \(k\in\mathbb F_{13}\); vertical lines each contain exactly one graph point. The independent program also retains the exact endpoint histogram.

| Scope | Count |
|---|---:|
| Moment-matched normalized array pairs | 49,895,820 |
| Exactly three values at shift 2 | 70,806 |
| Exactly three values also at shift 3 | 11,479 |
| At least 11 shifts of image size at most 3 | 5,418 |
| Of the latter, both graphs of cap at most 4 | **0** |

The 5,418 rows are normalized tuples, not distinct affine classes. The endpoint histogram leaves 5,310 rows with one graph affine (maximum line intersection 13), and 108 rows with both graphs nonaffine. The 108 rows all have exactly 11 small-image shifts; 66 have cap pair (9,9) and 42 have cap pair (10,10). Thus even cap eight is ruled out for both nonaffine graphs. No square-sum condition, determinant assumption or floating-point calculation was used.

## 4. Conclusion and replay

The two slices \(f_s,f_{s'}\) from Section 1 satisfy input moments and cap four, while Section 2 gives the pairs of shifts forbidden by Section 3. This contradicts four full layers of a size-52 spectral pair. By input 2 it also contradicts any full 13-point plane. An affine plane with eleven points would determine at least ten directions by the prior exact (11,9,4) gate. That gives either the already excluded layered case or the impossible 10-uniform pencil; thus the remaining plane cap is ten. No nonlayer implication beyond these bounds is made here.

Replay (tested C++17 and Python standard libraries, no external data):

~~~sh
g++ -O3 -std=c++17 -Wall -Wextra layer_pair.cpp -o layer_pair
./layer_pair | diff - layer_pair.out
g++ -O3 -std=c++17 -Wall -Wextra independent.cpp -o independent
./independent | diff - independent.out
python3 check_minors.py | diff - check_minors.out
~~~
