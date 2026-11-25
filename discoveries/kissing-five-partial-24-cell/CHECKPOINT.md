# A partial 24-cell forces the 40-point bound in dimension five

Private research checkpoint, 25 November 2026. No global kissing-number improvement and no new 40-point construction.

**Theorem.** A kissing configuration in R^5 containing at least 22 vertices of a single equatorial regular 24-cell has at most 40 points. If a 40-point kissing configuration contains at least 23 vertices of a single equatorial 24-cell, it contains all 24 and is isometric to D5 or Leech's L5. The 22-point hypothesis covers every missing pair, including adjacent vertices. The equality classification is only claimed for 23.

## Source frontier and scope

The first unsettled dimension is five; see the [Boyvalenkov–Dodunekov–Musin survey](https://arxiv.org/abs/1507.03631) and [Musin's dimension-four theorem](https://arxiv.org/abs/math/0309430). The best established construction is 40, attained by several nonisometric codes. The refereed upper bound is 44 from [Mittelmann–Vallentin](https://arxiv.org/abs/0902.1105); [Tomar's November 2026 preprint](https://www.mefyi.com/papers/kissing-numbers-five-seven/) gives an exact rational certificate of 43, the strongest bound located. [Nyam-Osor's formalized results](https://aixiv.online/abs/2611.01607) restrict antipodal pairs of hypothetical 41–43 codes, not the geometry proved here. The present theorem only concerns codes with a prescribed common equatorial partial 24-cell.

[Cohn–Rajagopal, §2](https://arxiv.org/html/2412.00937v3) describe the two full 24-cell layer constructions and their triality; their Q5 and R5 demonstrate why a blanket D5/L5 claim about all 40-point codes would be false. [The earlier D4-containing obstruction](https://github.com/wustep/maths/blob/main/problems/kissing-5d/ATTACK.md) already bounds the full 24-cell extension by 40 using the orthoplex bound; nearby searches there and in [Kriebel's research dossier](https://github.com/AlecKriebel/Math/tree/main/kissing_number_5) use several distinct exact or grid-based hypotheses. The missing-one and missing-two continuous extension statements below were not found in a bounded targeted search. No exhaustive priority claim.

## 1. Coordinates and high/low split

Use the equivalent normalization in which every point has norm squared 2 and distinct points have inner product at most 1. In the equatorial R^4 let D be the 24 roots ±e_i±e_j for i<j. Retain D minus a set M of k=1 or 2 distinct roots. Every extra point is x=(u,t), with u∈R^4, ||u||²+t²=2. Its side is upper if t≥0, lower if t<0; in particular t=0 is counted once. Put τ=|t|. Call x high if τ≥1; call it low otherwise. A low must violate a missing-root inequality since ||u||>1 and the full D-polar inequalities are |u_i|+|u_j|≤1. In fact if a=max_i |u_i|≤1/2 then ||u||²≤1; if a∈[1/2,1] the other absolute coordinates are at most 1−a, so ||u||²≤a²+3(1−a)²≤1. Equality gives exactly the coordinate axes ±e_i and the sixteen half-cube vectors (±1/2)^4. Denote these 24 unit directions by F.

Each side has at most 8 high points. If u=0 for one, its vertical magnitude is √2 and it conflicts with every other high on that side. Otherwise the horizontals of any two highs have dot ≤1−ττ'≤0; their unit directions in R^4 are pairwise nonpositive. For m such directions with Gram matrix G, tr(G²)≥m²/4 by rank and trace, and ΣG_ij=||Σ directions||²≥0. For off-diagonal g∈[−1,0], g²+g≤0, so m²/4≤Σ(G_ij²+G_ij)≤2m. Equality m=8 requires every off-diagonal 0 or −1, forcing four orthogonal lines each containing both signs. Every high then has an orthogonal horizontal neighbor, so their full inner product of ττ'≤1 forces τ=1, and each horizontal has norm one. The eight horizontals are therefore a full orthoplex of unit directions in the polar of the retained roots.

If the missing roots do not form an antipodal pair, the ± horizontal pairs make the whole orthoplex compatible with D, since at least one representative of every antipodal pair of D was retained. If M={r,−r}, for the four unit basis directions w_i, Σ_i (r·w_i)²=||r||²=2, so at most one basis direction can fail |r·w_i|≤1. The other three lines lie in F. In F a coordinate axis is never orthogonal to a half-cube direction; half-cube vectors are orthogonal exactly at Hamming distance two, so their orthogonality requires equal parity of sign word. Hence the three valid lines belong to one of three fixed orthoplexes (coordinate, even half-cube, odd half-cube); their orthogonal complement is its fourth line. This restores the missing condition for the fourth. In all cases an eight-point high layer is exactly one of these three frames at height ±1.

## 2. An eight-point high layer forces its own low points to be exact missing roots

Rotate one of the three triality frames to the coordinate frame, preserving D. This is possible explicitly: a one-coordinate sign flip swaps the even/odd half-cube frames, and the normalized symmetric 4×4 Hadamard matrix with columns (++ ++), (++−−), (+−+−), (+−−+) swaps the coordinate and even frames; the latter permutes D since each sum/difference of distinct columns is a root. On the same side as this high layer any low has |u_i|≤1−τ from its contacts with (±e_i,±1), with the vertical sign chosen for that side.

Among the three perfect matchings of the four coordinate indices, some has all its root inequalities retained: the six unordered coordinate pairs partition into those three matchings, and at most two pairs are affected by deletion. For that matching, both horizontal pairs have sum of absolute values at most 1. If τ>1/2, then ||u||²<1, impossible. For 0≤τ≤1/2 put c=1−τ. Two nonnegative numbers each at most c and with sum at most one have squared sum at most c²+(1−c)² (maximize over the vertices of the truncated square). Therefore
    2−τ² = ||u||² ≤ 2((1−τ)²+τ²).
For 0<τ≤1/2 this contradicts 5τ²−4τ<0. Thus τ=0. Equality in the pair bounds now forces each pair to consist of absolute coordinates (1,0) in some order, so u is a D-root. Because it is an extra point it must be missing from the retained core. This proves the assertion before any full 24-cell insertion.

## 3. A single-root cone and all nonadjacent deletions

For a named missing root r, take coordinates r=e1+e2. Suppose a low violates r, and all eight roots e_i±e_j with i∈{1,2}, j∈{3,4} are present or obeyed. Then with p=u1+u2 and m=max(|u3|,|u4|), p>1, 0<u1,u2≤1 and m≤1−max(u1,u2). These horizontal conditions suffice. For p=2 the point is the apex (r,0). Otherwise set α=p−1 and write u=αr+(1−α)h. Here
    h1+h2=1; h1,h2≥0; |h3|,|h4|≤min(h1,h2).
The base described by these inequalities is the convex hull of e1,e2 and the four vectors (1/2,1/2,±1/2,±1/2). Indeed at width λ=2min(h1,h2), split weight λ into the square at (1/2,1/2), and the remaining weight into the axis e1 or e2 corresponding to the larger of h1,h2. Each base vertex has norm one. Viewing (u,2−p) as the same convex combination of (r,0) and (h,1) gives norm at most √2, hence τ≥2−p. In particular τ≥|u1−u2|+2m.

More generally whenever a low has p>1 and τ≥|u1−u2|+2m with τ>0, it belongs to the cone generated by the seven rays
    (r,0), (e1,1), (e2,1), and ((1/2,1/2,±1/2,±1/2),1).
To see this take A=(p−τ)/2>0 and H=(u−Ar)/τ. The last inequality exactly gives H in the stated base, so (u,τ)=A(r,0)+τ(H,1). Reflect vertically for the lower side. The seven generators each have norm √2 and every pair has dot at least 1. A nonapex low has τ>0 in this case, hence its representation includes the apex and a distinct lifted ray with positive coefficients. If a vector x of norm √2 is a nonnegative combination of such rays, its coefficient sum is at least 1 by the triangle inequality, strictly greater if two distinct rays have positive coefficients. Therefore the dot of two cone points on this sphere is at least the product of their coefficient sums, with strict >1 when one is multiray. The apex is the only ray at τ<1. No two distinct lows of this type can share a side.

The eight sufficient roots in this proof each have dot +1 with r. Thus if there is only one hole, or if two holes have dot −2, −1 or 0, every low can be assigned to one violated missing root and is in its seven-ray cone on its own side. This bounds lows by k on each side. Adjacent holes are the only remaining case.

## 4. The adjacent missing pair: one extra ray per cone

For raw dot +1, signed permutation gives r=e1+e2 and s=e1+e3. Assign a low to whichever missing root has the larger horizontal dot, breaking ties consistently. It suffices to consider the r chamber, p=u1+u2≥q'=u1+u3 and p>1; swap coordinates 2 and 3 for the other chamber. If q'≤1, every needed inequality in the preceding seven-ray proof holds, since the only deleted cross root for r was s. If q'>1, put (a,b,c,d)=u and κ=1−a. From the retained coordinate-4 root pairs and b+c≤1, and the two violations, we have
    a≥b≥c>κ≥0,  b+c≤1,  |d|≤κ.
Add the ray W=((1,1/2,1/2,0),1/√2) to the r cone. Its norm is √2 and its dot with every old ray is strictly greater than 1. We show the chamber is covered by this augmented cone.

Use the lifted rays (r,0), W, (e1,1), and H_±=((1/2,1/2,1/2,±1/2),1). For each h∈[2κ,2c], assign respective coefficients
    b−c, 2c−h, a−b−c+h/2, (h+2d)/2, (h−2d)/2.
All are nonnegative by the chamber inequalities and their horizontal combination is u. The resulting vertical height T(h) is increasing, with slope 3/2−1/√2, and T(2c)=a−b+2c. At h=2κ the coefficient sum is a+κ=1, so the full convex combination has norm at most √2 and T(2κ)≤τ. If τ≤T(2c), interpolate h to reach τ exactly. If τ>T(2c), then τ≥|a−b|+2max(|c|,|d|), so the original cone criterion also gives membership. This proves coverage without a grid.

In the augmented cone all rays still have pairwise dots at least 1, and the only possible single-ray lows are the apex and W, which themselves have dot 3/2. Two distinct lows in one side assigned to r must therefore conflict, whether multiray or single-ray. The symmetric s case gives the same result. Thus every side contains at most k lows also for the adjacent pair.

## 5. Finish the count and equality

If neither side has a full high layer, there are at most 14 highs and at most 2k lows, so the total is at most 24−k+14+2k=38+k≤40. If either side has 8 highs, each low on that side is a missing root by Section 2. Once such a root is inserted, no other low assigned to it on either side is allowed (it would violate its inner-product bound). For each remaining missing-root type all assigned lows lie on the opposite side, at most one each. Thus there are at most k lows overall and at most 16 highs, giving at most 40. The t=0 side convention ensures the inserted vertices are not counted twice.

When k=1, the first branch has at most 39; equality requires both high layers of eight and one low, exactly the missing root. A full D4 equator with both of those modes is D5 if the two frames agree and Leech's L5 if distinct. The sign flip and normalized Hadamard above generate all permutations of the three frames, so these are the only two isometry types. They differ, for instance, in the presence of antipodes outside D. Both are genuine 40-point codes, directly checked below. For k=2, equality either restores both missing roots with the two full high layers, or leaves the precise high/low split 7+7 and 2+2; the latter is not classified here. In particular the theorem gives no global upper bound below 43. Compactness only gives a qualitative stability neighborhood: there is some positive distance from a 41-code to retaining an aligned 22-root subset (no numeric radius is asserted).

## Exact controls

    python3 verify.py

The standard-library verifier checks the 276 missing pairs by raw dot, all polar frame identities and triality actions, the cone ray products over Q(√2), and the mixed ray W as an actual two-violation low compatible with all 22 retained roots. It checks exact sharp D5 and L5 witness coordinates and their distinct antipode counts. The continuum proof above uses these finite identities plus convexity and a rank bound, not a finite sampling of sphere points.
