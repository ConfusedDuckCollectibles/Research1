# Size 78 in \(\mathbb F_{13}^3\): exclusion when neither member has six full parallel planes

Step 1 of the size-78 exclusion. **In any spectral pair of size 78, at least one member is supported in six parallel affine planes, each with thirteen points.** This is an exclusion of the entire branch with neither member layered, not exclusion of size 78 or proof that both members must be layered. The case in which one or both members are layered remains.

This step also records the basic reductions used in the remaining steps. The basic direction, profile and moment identities follow [Ferguson–Iosevich §§8, 13](https://arxiv.org/html/2610.30712v1) and [Fallon–Mayeli–Villano](https://arxiv.org/abs/1902.02936). The weighted identity is [Ferguson–Iosevich Proposition 8.6](https://arxiv.org/html/2610.30712v1#S8.SS1), used here with a deficit for the largest pencils and a different prime. [Malikiosis](https://doi.org/10.1007/s43670-023-00078-2) provides the broader LP context; its published unconditional threshold does not exclude 78. No claim of priority or publication.

## 1. Basic reductions for an entirely nonlayered pair

Let \(U_E\) (respectively \(U_A\)) be the projective set of normals with all parallel affine plane occupancies six, and \(V_E,V_A\) their complements. Write \(D(E)\) for projective directions of nonzero differences. Cyclotomic character orthogonality gives \(D(E)\subset U_A\) and \(D(A)\subset U_E\). A pencil of the 14 affine planes through a line containing \(a\) selected points has total occupancy \(78+13a\). Thus one cannot have 14 uniform directions in any pencil. A 14-point plane would determine all 14 tangent directions by pigeonhole, forcing such a pencil on the partner; therefore the plane cap is 13. Every pencil has at least one uniform direction because the partner has more than 13 points, hence two points with the same projection onto the perpendicular parameter. In particular line cap is six.

The exact planar inputs in the replay are:
* A six-point planar set of line cap five needs at least five directions.
* A ten-point planar set of line cap five needs at least eight directions.
* An eleven- or twelve-point planar set of line cap five needs at least ten directions.

If a member has a six-point line, the pencil through it has at most three uniform directions, since \(156\le 6K+13(14-K)\). The partner's 13 parallel slices in the corresponding tangent directions therefore determine at most three directions in total. Each slice has at most six points: any seven-point slice under line cap six contains a noncollinear six-subset, which determines at least five directions by the first gate. Their average is six, so all are six-point collinear, with at most three slopes. Two of those lines are parallel and their common plane determines at least 12 directions (one common slope plus at least eleven cross-slopes by Cauchy–Davenport on two six-element sets). The original member has a pencil with at least 12 uniforms; by the Radon lemma below it is layered. The partner has actual six-lines, so applying this argument to it yields both layered. Hence in the entirely nonlayered branch both line caps are five.

We also need that **each member has a slice of size at least 12**, even before separating the pencil cases. If \(E\) had plane cap 11, its uniform pencils would have at most nine points: \(K\ge11\) is ruled out by the Radon affine-profile lemma below, and the exact \(K=10\), cap-11 check has no projected solution. Also its line cap is at most five by \(78+13a\le14\cdot11\). If its partner had a six-line we would get two parallel six-lines in \(E\) by the preceding argument, contrary to plane cap 11. Thus the partner also has line cap five. An eleven-point slice of the partner would give the original \(E\) a pencil with at least ten uniforms by the 11-point planar gate. So the partner's plane cap is at most ten, and both uniform-normal sets have pencil cap nine. Both meet every pencil in at least five uniforms: the partner has a parallel six-point slice by averaging, to which the first gate applies. Define for a normal \(v\) of \(E\)
\[
 k_E(v)=\frac{q_E(v)-195}{13}
       =\frac{1}{26}\sum_{t\in\mathbb F_{13}}(n_t-6)^2 ,
\]
where \(q_E(v)=\sum_t\binom{n_t}{2}\) counts equal-level pairs. Then \(k_E\ge0\), and its support is \(V_E\). As there are 14 orthogonal normals for every difference direction and exactly one in any fixed pencil unless that direction is the pencil's perpendicular, the pair counts give
\[
 \sum_v k_E(v)=489,\qquad
 \sum_{v\perp d} k_E(v)=21+t_E(d),
\]
where \(t_E(d)\) counts unordered pairs of direction \(d\). For every \(d\in V_A\), \(t_E(d)=0\). Consequently
\[
 21|V_A|=\sum_v k_E(v)|V_A\cap v^\perp|. \tag{1}
\]
Here pencil cap nine would give \(|V_A|\ge117\), so \(|U_A|\le66\). Counting the disjoint sections of \(U_A\) on 14 lines through an external point gives \(|U_A|\ge70\), contradiction. Exchange \(E,A\). Now every nonlayered member, having line cap five, has a slice of at least 12 points, whose 12-point subset's tangent directions force a pencil with at least ten uniforms on the partner.

### Exact projection lemma for large pencils

Project along the perpendicular direction \(d\) of a pencil with \(K\) uniforms; choose a uniform slope vertical and put \(L=14-K\). If \(g_i\) are the integer profiles in the remaining slope families, their sum at each quotient point satisfies
\[
 13 f(x,y)=\sum_{i=1}^{L}g_i(y-s_i x)-6(L-1).
\]
Modulo 13, compare homogeneous coefficients in the reduced profiles. The Vandermonde relations give for a positive degree \(j\) the coefficients \(w_i Q_j(s_i)\), \(w_i=\prod_{k\ne i}(s_i-s_k)^{-1}\), \(\deg Q_j\le L-j-2\). Thus the degrees are at most \(L-2\), as in the earlier checkpoint.

* If \(K\ge12\), all profiles have constant residues. Nonuniform constant-residue profiles can only take values 0 and 13 with exactly six entries equal to 13. The constant Radon congruence forces \(L=1\), and \(f\) is exactly the indicator of six full parallel lines. This is the excluded layered branch. Under cap 11 even this is impossible.
* If \(K=11\), there are three affine residue profiles. Either all slopes are zero, impossible by the constant equation \(0\not\equiv12\), or all slopes are nonzero by the one-dimensional Vandermonde kernel. Each nonconstant profile is a permutation of \(0,\ldots,12\); its least residues already sum to 78, so no value is elevated to 13. In coordinates \(r=g_1\bmod13,\ s=g_2\bmod13\), the third profile is the least residue of \(12-r-s\). The integer Radon identity gives the Boolean triangle \(T=\{(r,s)\in[0,12]^2:r+s\ge13\}\). Under cap 11 this also is impossible.
* If \(K=10\), the four profiles have degree at most two. The exact finite check enumerates all three projective orbits of four exceptional slopes. For nonzero quadratic terms, translate to center all simultaneously since the linear/leading ratios are affine in \(s_i\); enumerate the common quadratic multiplier, admissible profiles with all possible zero-residue elevations, and constants. In the affine branch, a nonzero linear numerator vanishes in at most one of the four positions, and two nonzero rows may be translated to zero constant; the rest are enumerated. All-zero residue profiles contradict the constant equation. Test all integer Radon recoveries with \(0\le f\le\lfloor4(13-6)/13\rfloor=2\), verifying all fourteen line families on acceptance. There are no cap-11 projections, and at cap 13 exactly 804 normalized projections, of which exactly 24 are Boolean. All Boolean projections have exceptional slope maxima \((12,10,10,10)\). The verifier follows the earlier audited atlas code and regenerates these graphs rather than taking the counts as input.

## 2. Exact labelled six-fibre gate and the K10/mixed exclusions

The executable normalization for the planar gates sends two missing directions to the coordinate axes and any two points to \((0,0),(1,1)\). Extend the allowable nonzero direction mask to the fixed bound, then DFS adds points in increasing \(x\), tracking \(y\)-injectivity, directions and line caps. The threshold of six points in at most five directions gives 840 normalized arrays, all with exactly five directions; their injective projections along every missing direction yield exactly five affine abscissa mask types:
\[
 111,\ 187,\ 215,\ 311,\ 427 .
\]
Here a subset of \(\mathbb F_{13}\) is represented by its 13-bit mask, canonically minimized over all affine changes of coordinate. Computing these masks from the point arrays is part of the replay.

If an injective 78-point shadow has a quotient slope whose size-six row gives a mask outside these five types, the actual lifted planar six-set on that row must determine at least six directions: it has line cap five, and the labelled projection is along a missing direction of its planar graph. Direction lines for distinct quotient slopes meet only in the projection direction \(d\), and Booleanity says \(d\notin D\), so these direction lower bounds add without overlap.

In a pair with maximum \(K\le10\) on both sides, choose a 12-point subset of an actual affine slice of \(A\) of size at least 12. Its tangent direction line has at least ten points of \(D(A)\), so \(U_E\) has \(K=10\) there. The perpendicular projection direction \(d\) is the normal of the actual plane with at least 12 selected points, hence not in \(U_A\). Then \(d\notin D(E)\) by spectrality, and the corresponding \(K=10\) shadow of \(E\) is Boolean. Each of the 24 exact candidate projections has some forbidden size-six abscissa mask in every one of its ten uniform directions (at least seven such rows per direction). The three size-10 peaks give at least eight directions each, and the size-12 peak at least ten. Thus \(|D(E)|\ge10\cdot6+3\cdot8+10=94\). But if \(K\le10\) on the partner, \(|V_A\cap v^\perp|\ge4\) throughout (1), so \(|V_A|\ge94\) and \(|U_A|\le89\), contradicting \(D(E)\subset U_A\).

If either has a pencil with \(K=11\), its triangular shadow is automatically Boolean. The replay checks that all eleven uniform slopes have a forbidden size-six abscissa mask. The three 12-point slope peaks give at least ten directions each. Consequently \(|D(E)|\ge11\cdot6+3\cdot10=96\), ruling out a partner with \(K\le10\) as just above. Since \(K\ge12\) implies layers, the only prospective nonlayered branch is \((K_{\max},K_{\max})=(11,11)\).

## 3. Both K11: weight deficit forces global affine profiles and a divisibility contradiction

First recall the global eighth-moment step. For either \(U\), every projective line contains between five and eleven uniform directions: the minimum uses a size-at-least-six slice of the line-cap-five partner. If a nonzero degree-eight homogeneous polynomial vanishes on \(U\), each nonfactor line contributes at most zero to \((n_L-5)(n_L-8)\), and each of at most eight line factors contributes at most \((11-5)(11-8)=18\). With \(N=|U|\), projective incidence counting instead gives \(N^2-169N+7320\ge719/4>8\cdot18\). Therefore all directional moment forms through degree eight that vanish on \(U\) are identically zero (multiply by a linear power to degree eight as needed). In fact, equidistribution forces their vanishing on \(U\). For any normal, its reduced residue profile now has degree at most three: multiplying by \(t^j\), \(0\le j\le8\), and summing over \(t\) kills successively each coefficient of degree \(12-j\), using \(\sum t^i=0\) except positive multiples of 12.

The exact profile gate exhausts all \(13^4\) reduced cubics (including lower degrees), all possible elevations of zero residues to 13 subject to integer mass 78, and no elevations of positive residues. For a fixed polynomial this only involves the number of elevated zeros, since replacing 0 by 13 adds 13 to \(\sum(n-6)^2\) and all such replacements have identical weight. The zero polynomial yields exactly \(\binom{13}{6}=1716\) six-full-layer profiles, inadmissible in this branch. Of the remaining 23,765 exact profiles including the uniform one, the maximum \(k(v)\) is 12. No nonconstant polygon geometry beyond the base occupancy cap is assumed in this numerical gate.

Both triangular direction lower bounds give \(|V_A|,|V_E|\le183-96=87\). In (1) compare the weighted line loads with \(4\sum k_E=1956\); the actual total is \(21|V_A|\le1827\). Only a line with \(|V_A\cap v^\perp|=3\) can supply a positive deficit relative to four (loads are at least three). Each such term saves at most \(k_E(v)\le12\). Therefore the deficit at least \(1956-1827=129\) forces at least eleven distinct 11-uniform pencils in \(U_A\). This count is over \(v\in V_E\), and polarity gives different projective lines for distinct \(v\).

For \(j=9,10\), the homogeneous moment polynomial \(M_j^A(v)=\sum_{a\in A}(v\cdot a)^j\) vanishes on \(U_A\). Each of its eleven 11-uniform pencils has 11 distinct projective roots on its line; a nonzero restriction of degree at most ten cannot do so. Thus each line must be a factor if \(M_j^A\ne0\), impossible for eleven distinct factors. These moments vanish everywhere. Combined with the eighth-moment step, each reduced occupancy profile of \(A\) is now affine. The constant-zero residue is again impossible because it represents six full layers; any other constant residue is constant six. Nonconstant affine profiles are precisely permutations of \(0,1,\ldots,12\) without elevations and have weight \(2\sum_{j=1}^6 j^2/26=7\). Thus \(\sum k_A=489\) would be a multiple of seven, contradiction. This excludes both K11 and completes the entirely nonlayered proof.

## Replay and residual

Run from this directory
~~~sh
node check.js | diff - expected.txt
~~~
Standard-library exact integer arithmetic only. Default V8 replay checked all four output lines, no network or solver. The exact gates are finite cap-conditioned necessary reductions and the positive arrays were not claimed spectral. Separate proof dependencies for the standard incidence/pencils and the source of the exact slope normalizations are above. Any subsequent argument for six parallel full layers is a separate case.
