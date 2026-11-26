# Excluding spectral sets of size 204 in F_17^3

Research note, 26 November 2026. **No 204-point spectral pair exists in F_17^3.** This removes 204 from the checked cardinality ledger, leaving 34,68,102,136,153,170. It does not settle the full field. This is exact proof with finite certificates, not externally reviewed and no assertion of priority.

Use p=17, m=12. Regard E as a subset of V=F_p^3 and A as one of V*. Write D for the projective directions of differences, and U(E) for normals whose p parallel planes all have m points. Orthogonality of the square character matrix and of its transpose gives D(A) subset U(E), and vice versa. For a nonzero normal the Fourier sum vanishes exactly at uniformity, since the integer polynomial of plane counts vanishes at a primitive pth root and must be a multiple of Phi_p=1+X+…+X^(p−1). Thus the nonzero Fourier support is a union of punctured lines, regarded projectively below.

## 1. Dense cardinality forces full parallel layers

We use the standard elementary facts of [Fallon–Mayeli–Villano](https://arxiv.org/pdf/1902.02936) and [Ferguson–Iosevich §§8–13](https://arxiv.org/html/2610.30712v1). For an affine line containing a points of E, the p+1 planes through it have total mass mp+pa. A full pencil of uniform normals is impossible since it would give pa=m. A plane with more than p points determines every tangent direction, by pigeonholing into its p affine lines of that direction; the partner would then have a full uniform pencil. Hence plane cap p on both sides. In the pencil through any line parallel to a direction d, there is a uniform normal: because |A|>p, two partner points have the same value at d, so their difference is such a normal. Thus mp+pa≤m+p², giving line cap five when p=17,m=12.

We spell out the conditional LP step in [Malikiosis §6 Eq. (6.1)](https://link.springer.com/article/10.1007/s43670-023-00078-2#Sec6). If the projective nonuniform set of E is 3-fold blocking, take an inclusion-minimal 3-fold blocking subset B, with b=|B|. [Bishnoi–Mattheus–Schillewaert, Theorem 1.1](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v25i4p66) gives b≤17 sqrt(46)+20<136, so b≤135. The strict numerical inequality is 46*17²=13294<116²=13456. Malikiosis's explanatory display of the size bound appears to omit the +3, while Eq. (6.1) agrees with the full theorem because the formula uses b−3 below. We use the full established bound.

On V* define h(0)=b−3p, h(v)=1 for v≠0 whose direction is in B, and zero otherwise. Its unnormalized Fourier transform is p*(|B∩y^perp|−3) for y≠0 and p*(b−3) at zero: a punctured projective line contributes p−1 if orthogonal to the character and −1 otherwise. The transform is nonnegative. Because D(A) avoids B, summing h(a−a') gives |A|h(0). Fourier inversion and retaining only the zero-frequency term give

    |A| h(0) >= |A|² hhat(0)/p³,
    |A| <= p² (b−3p)/(b−3) <= 289*84/132 < 204.

The ratio increases for b>3 and a 3-fold blocker has b>3. Thus there is a projective line containing at least 16 uniform normals for E. Along its annihilator direction the projection injects E: any line in that direction has mp+pa≤(p−1)m+2p, so a≤1. Its Boolean planar shadow has at most two exceptional Fourier directions, and so, after linear planar coordinates, the indicator has the form a(x)+b(y) by Fourier inversion. If a is nonconstant two rows have fixed nonzero difference and Boolean entries force it to be ±1; both rows are constant, forcing b to be constant. This yields m parallel full lines in the shadow and therefore m parallel full layers (size 17 each) in E. The same argument applies to A. This is the known parallel-plane criterion, not a new application of the unconditional Theorem 1.4 at 204.

## 2. Choosing mutually perpendicular normals

The first finite gate establishes that **every 17-point planar set in F_17² with line cap five determines at least 14 directions**. Select five missing directions if they exist; send an ordered three to 0,1,infinity and enumerate all normalized masks. There are 105 masks and four PGL orbits. The exact cap-bounded graph search on them has 17,316; 16,003; 14,868; 14,377 recursion nodes respectively and no full graph; see Section 5 for completeness. This finite bound is stronger in this parameter setting than the general Rédei–Megyesi bound (see [Somlai 2024](https://doi.org/10.1007/s00013-024-01979-x)).

Choose any layer normal n_A∈V for A. It is nonuniform and therefore absent from D(E); projection of E along n_A is injective. The shadow Q in W=V/<n_A> is Boolean of size 204. Its dual W* is exactly the tangent direction plane of A's layers in V*. The at least fourteen secant directions of any single 17-point layer of A are uniform normals on Q. Thus its complement Q' of size 85 has at least fourteen uniform normals as well.

The needed conclusion that Q' is a union of parallel lines has a simple exact univariate check. Send a uniform normal to vertical and let L≤4 be the number of actual nonuniform finite slopes s_i, with corresponding line occupancies g_i(t) in [0,17] of mass 85. The Radon identity is

    sum_i g_i(y−s_i*x)=5(L−1)+17*1_Q'(x,y).

It counts every other point once and the reference point 18 times, then removes the uniform contributions. Reduce the identity modulo 17 and write each profile as a polynomial of degree <=16. The homogeneous degree-j coefficient equations for j>0 are sum_i [t^j]g_i * s_i^h=0 (0≤h≤j), so Vandermonde kills all positive coefficients of degree at least L−1; the degree is at most max(0,L−2)≤2. Nonconstant affine polynomials have least nonnegative residues summing to 136. All 272 centered nonconstant quadratics a*t²+c have least-residue sums at least 102; translation centers any such quadratic. Actual counts are at least the least residues, so neither fits mass 85. Thus every exceptional polynomial is constant modulo 17. Nonzero constant residue would give unique equal counts in [0,17], necessarily five, making it uniform. There is at least one exceptional slope, as uniformity in all directions would give 1_Q'=5/17 by the full Radon identity. Hence one exceptional profile is zero modulo 17 and uses only counts 0 and 17: Q' is five full parallel lines. Q is its twelve-line complement. This same Boolean shadow implication also follows from [Ghidelli's rich/poor direction theorem](https://arxiv.org/abs/1903.03881) and [Wang's four-special classification, Table 3](https://arxiv.org/html/2609.20495v2#S7), but neither is needed.

Pulling Q back gives an E layer normal n_E in the annihilator of n_A, so n_E(n_A)=0. This applies whether the original normals are perpendicular or not, without a uniqueness assumption. Because n_E is nonuniform for E, A projects injectively along n_E. Choose a basis of V with second vector n_A, first vector independent of it in ker n_E, third vector evaluated to one by n_E, and take the paired dual basis on V*. These coordinate changes preserve spectral pairing. The two injectivities give

    E={(x,f_s(x),s): x∈F_17,s∈S},
    A={(u,t,g_t(u)): u∈F_17,t∈T},  |S|=|T|=12.

Each layer is a full graph over the field, and none of the graphs is affine by line cap five.

## 3. Cross-layer chart palettes and many active shears

Let zeta be a primitive 17th root and put H_s(r)=sum_x zeta^(r*x+f_s(x)). Call r active if at least one H_s(r) is nonzero. The 17-point planar bound from Section 2 gives at least 14 active r: a single E layer determines ≥14 finite slopes, and slope −r is missing exactly when f_s(x)+r*x is a permutation, exactly when H_s(r)=0 by Phi_17.

A cross-layer direction in A proportional to (r,1,w) supplies orthogonality sum_s H_s(r)*zeta^(w*s)=0. For active r, this is the Fourier transform on F_17 of a nonzero complex vector supported on at most twelve positions. The [Tao–Chebotarev prime cyclic uncertainty inequality](https://arxiv.org/abs/math/0308286) gives at most eleven values of w where it vanishes. Thus at each active r the cross chart has a single palette of ≤11 finite slopes. This is a complex Fourier statement, distinct from the next rank assertion over F_17.

For any x and active r, select from A the twelve points with u=x+r*t and view them in their affine plane as the graph

    {(t,g_t(x+r*t)):t∈T}.

It has line cap five. All its slopes are in that same fixed chart palette. Both the abscissa set T and the palette stay fixed as x varies.

## 4. Small affine-rank atlas and determinant contradiction

The second complete finite gate proves: **for fixed 12-element abscissae T and fixed allowed palette at most eleven finite slopes, the cap-five graphs normalized by subtracting their height at a fixed t_0 number at most four and lie in an affine space of dimension at most two over F_17.**

Fix active r and interpret G_r(X)=(g_t(X+r*t))_{t∈T} as a vector of reduced F_17 polynomials of degree ≤16. All values of its row-normalized vector G_r(x)−g_t0(x+r*t0)*1 lie in an atlas cell. It follows:
1. Any pairwise coordinate difference takes at most four values of x.
2. The positive-degree coefficient rows of G_r have rank at most three. Indeed after subtracting the common t_0 coordinate, all values lie in an affine flat of dimension ≤2. Modulo its direction subspace the polynomial vector has constant evaluation on the whole field; since reduced representatives are unique, its nonconstant coefficients vanish in the quotient. Restoring the common scalar adds only the span of the all-one vector.

If two different t,u had deg g_t,deg g_u≤4, their shifted difference would have degree ≤4 for every active r but take at most four values. Any nonconstant polynomial with that degree over F_17 has image of size ≥ceil(17/4)=5 (each fibre has at most four roots). Thus the difference is constant in x for every active r. Comparing two distinct active r, and replacing x+r*t by y, gives g_u(y+d)−g_u(y) constant for all y, where d=±(r_1−r_2)(u−t) is nonzero. Iterating a nonzero translation generates the prime field, so g_u is affine, a contradiction. At least eleven g_t therefore have degree at least five.

At an active r, columns of distinct positive degrees are independent since translation preserves leading terms and the columns can be eliminated in descending degree. The positive coefficient matrix has rank at most three, so at most three distinct degrees occur among those eleven. Four g_t share a degree d≥5. Use their four columns and rows d,d−1,d−2,d−3. The determinant is a polynomial in r of degree ≤0+1+2+3=6. The coefficient of r^6 is the product of the four leading coefficients and binom(d,k), k=0,1,2,3, times a Vandermonde in the four distinct t's: in row d−k and column t, the leading term is lead(g_t)*binom(d,k)*t^k*r^k. Contributions of lower degree terms of g_t are lower degree in r. The leading coefficient is nonzero because 5≤d≤16, t's are distinct, and the modulus is 17. This contradicts vanishing at all ≥14 active r required by rank at most three. Thus 204 is excluded.

A reusable general version of the endgame: for an analogous layered pair with maximum row alphabet B and affine rank D, shifted coefficient rank is at most D+1. Degrees ≤floor((p−1)/B) can occur in at most one layer once two different active slopes exist. If the remaining layers outnumber (D+1)^2 and their degrees permit D+2 positive leading coefficient rows, a determinant of degree binom(D+2,2) obstructs more than that many active slopes. No untested condition on another field is claimed.

## 5. Finite certificate completeness and replay

Run standard-library Node file `verify-f17-204.cjs`. The optional --full switch emits all 353 atlas heights. See expected.json for exact output.

For the full 17-point graph gate, every five-element chosen missing set can be moved to contain 0,1,infinity. The code traverses the full normalized projective orbit by each ordered three preimages of 0,1,infinity. A projective map of slope labels is induced by a linear change of coordinates, preserving cap and direction count. Missing vertical implies a unique graph over all 17 x; normalize f(0)=0 by height translation. Three positive cap-five controls with exactly four missing directions verify direct line occupancies and show sharpness in this specific scope.

For the 12-point atlas, vertical is already missing. Pad the allowed finite palette to eleven if necessary. Its complement consists of six finite slopes. The transformation (t,y) -> (a*t+b,a*y) with nonzero a preserves all slope labels and gives 25 abscissa orbits covering all 6,188 sets. Independently (t,y) -> (t,a*y+b*t) changes slopes by any affine permutation without moving t, and gives 49 slope-complement orbits covering all 12,376 complements. This proves the Cartesian product of catalogues suffices. Vertical translations normalize the first height to zero. Changing t-representative, slope representative, or normalization base changes height vectors by invertible affine transformations of their normalized coordinate space (equivalently use the quotient by the common-constant direction). Count and affine dimension are preserved.

For each labelled problem, the DFS assigns one abscissa at a time, possibly choosing the smallest remaining domain. Each assigned point excludes exactly the candidates making a forbidden slope. It prunes only upon losing all possible heights at a required abscissa, or inserting a sixth point on a nonvertical line. Vertical cap is automatic. Every allowed normalized full vector is emitted once and its full line fibres are independently evaluated. Modular elimination of differences from the first vector yields the affine rank. The complete 25*49 atlas has 584,348 recursion nodes and 207 nonempty cells comprising 353 normalized shapes. Counts by (shape count, rank) are (1,0):84, (2,1):108, (3,2):7, (4,2):8. The original V8 run took 5.3 seconds. Independently coded forward-check graph search used a different 558,992 recursion nodes and matched all 207/353 census and rank counts.

Primary comparison: Ferguson–Iosevich settles p=11, not p=17; Malikiosis leaves the 204 full layers. Earlier Research1 p13 chart work used a fixed coordinate difference shared by compatible normalized shapes. Multiple p17 atlas cells have no such fixed pair, so the new varying-shear rank argument is essential. A bounded source search found no matching completed 204 primary theorem, which is not an absolute priority claim.
