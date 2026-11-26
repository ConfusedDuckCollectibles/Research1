# Four-petal sunflowers in triples: one exact construction gate and structural restrictions

Scoped research note, 26 November 2026. The general maximum remains between 39 and 49. No improvement to those unrestricted endpoints is claimed. All families are simple, all members have size three, and an empty common core counts. The forcing threshold is one more than the maximum.

## Primary baselines and scope

- Axante, Budala, Chitic, Dumitru and Nacu, [arXiv:2609.06175v1](https://arxiv.org/html/2609.06175v1), §§5.1–5.7: explicit 39-triple lower witness on 13 points and 49 upper (excludes 50), following earlier 38 and 69 in their Table 1. Their [public artifacts](https://github.com/bogdan27182/sunflower-paper) give computational source/ledgers for the upper. No newer primary f(3,4) endpoint was found. Their bound on the auxiliary matching-number-at-most-two family is 35 (Lemma 5.6).
- Kyle Wood, [arXiv:2611.13495v1](https://arxiv.org/html/2611.13495v1), proves a different extremum, ground-set span 39 (support counts distinct points, not sets). Its stated ensuing cardinality bound 130 does not improve Axante et al. This is not a 39-triple optimality proof.
- [OpenAI, An exponential bound for sunflowers, Nov. 24 primary version](https://github.com/openai/math/blob/3344da59991ed94f9ba711217b29c1d760badfa6/preprints/An-exponential-bound-for-sunflowers-November-24-2026/paper.pdf), Theorem 1.1, proves a fixed-petal exponential bound in maximum set size. It does not determine this finite maximum, and we do not reprove it. Babanskyy’s [November primary preprint](https://hexagonmath.org/pdf/2611.00552v1) concerns size four and three petals, not our reversed parameters.

## 1. Exact bounded trade result

Let F0 be precisely the family of Axante et al. Thm 5.1: three disjoint three-point bundles, the three bundle triples, nine triples made by a common point a and a pair inside a bundle, and 27 triples made by one point from each of two bundles and their pair-specific fresh point l_ij (the three l's and a distinct). No valid family of 40 contains at least 34 members of F0, even with arbitrary new ground points. This is an obstruction to one natural extension neighborhood, not a proof that f(3,4)=39.

Proof/computation. For k=0 to 5, delete exactly k members and add exactly k+1 different outside triples. The new triples together use at most 3(k+1) exterior points, so after relabelling enumerate all triples on 13+3(k+1) labels. The fact that some original point is absent after a deletion does not remove its original name from this domain. Encode each triple by a 32-bit integer and each selected old deletion index by a BigInt bitset; at the largest k the largest used label is 30, thus no signed wrap. Check F0 by all quadruples. List all 910 triples of old members sharing a core c, with union U. With one additional t this becomes a sunflower exactly when t intersect U equals c. For each t, recursively branch on the three available deletions for the first uncovered bad old triple, and when all are hit fill redundant deletions to exact k. This includes every size-k hitting set by induction and deduplicates. Invert to the deletion bins, each with the set of all individually valid insertions. Finally search increasing candidate indices. At each admission the only unchecked newly formed quadruples contain the new triple and at least one previously added triple, since validity with any three old triples is already certified. Scan all such quadruples. This is exhaustive on the prescribed domain; bit masks only encode finite names, not geometric or numerical guesses. If a family of more than 40 retained 34 or more, a 40-subfamily retaining 34 or more would contradict the same check.

Run `node verify-near39.cjs`. Observed row tuples (k, individually valid triples on that expanded domain, nonempty deletion bins, largest individual pool, bins with pool >= k+1, accepted partial nodes in joint search):
    
    0          0       0      0    0      0
    1          0       0      0    0      0
    2         18      54      1    0      0
    3         27    1998      3    0      0
    4         27   34038      4    0      0
    5        450  353349     20  297   1080

For k=5 the 450 are 126 entirely on old labels and 324 with one exterior label. The 324 are 18 old-pair patterns, each printed in all 18 fresh labels; no candidate with 2+ fresh survives individually. No witness for any k. Main V8 runtime about four seconds. Independent reconstruction and exact controlled first-occurrence fresh-name search matched the bin counts and excluded every joint candidate; its files are separately named INDEPENDENT-LOCAL39-*.

## 2. Analytic restrictions on a possible 49-member family

Use Axante et al. notation. For a triple T, N(T) is the number of family members meeting T, including T. For each x in T the exact-singleton link, consisting of pairs outside T joined to x, has maximum degree <=3 and matching number <=2, because T supplies another link edge. Thus it has at most 7 edges. Let s be the number of family members meeting T in exactly two points (not T). If maximum vertex degree is at most 9, summing the three exact-singleton sizes bounds N(T) by

    1 + s + min(21, 24 - 2s) <= 23.

In general Axante Eq. (5.1) gives N(T)<=25, including when Delta=10.

A family of size 49 must have a matching of size three (for example by their weaker Lemma 5.6); since 4 is forbidden, choose any such matching (T1,T2,T3). In fact every member can be an anchor: its disjoint family has at least 49-25=24 edges and cannot be intersecting, since intersecting families have at most 10 by their Lemma 5.4. For rooted-only block sizes a_i, double-anchor sizes c_i indexed by the missed anchor, and triple-anchor count z, the identity in their Eq. (5.8) is

    sum_i (N(T_i) + a_i - 2) = 2*49 - 6 + z = 92 + z.

Each a_i<=10 (the rooted blocks are intersecting), and a pair of size-10 blocks has c<=2 by their computer-assisted Lemma 5.9. The latter is applied only to the subfamily missing the third anchor, of matching number exactly two due to that disjoint third anchor. If all three N(T_i)<=23, two size-10 blocks would force the remaining anchor T_k to have N(T_k)>=49-10-10-2=27. Thus at most one a_i equals 10, and the left side is at most 3*23+28-6=91, impossible. Therefore every three-matching contains a triple meeting at least 24 members. In particular any family of 49 has maximum degree exactly 10, and every three-matching uses at least one degree-10 point.

More specifically, if N(T)>=24 write D=sum_{x in T}(10-d(x)). The exact degree count is N(T)=28-s-D. The separate singleton bound gives s>=2, so D<=2. In the extremal D=2, s=2 and all three singleton link graphs have seven edges; the degree patterns on T are (10,9,9) or (10,10,8). For the former its two exactly-two-point triples meet T in the two pairs through the degree-10 vertex; for the latter they share the pair of degree-10 vertices. Neither claim restricts the third points beyond what is proved upstream. This is a branching constraint, not a classification of possible 49-member families.

## 3. Matching-at-most-two cutoff 33, versus prior 35

Suppose 34 members and matching at most two (larger reduces to 34). They cannot all intersect. Choose two disjoint anchors A,B. There are no members disjoint from both; the A-only and B-only classes are intersecting and have sizes a,b<=10. Since N(A)=34-b<=25 and N(B)=34-a<=25, a,b>=9. A size-10 side would force at most six cross members by Axante Lemma 5.9, contradicting 34. Hence a=b=9 and N(A)=N(B)=25.

By Axante Lemma 5.5, each anchor has exactly three two-point meeting members, all with the same third point outside that anchor. Each exact-singleton link is R=K5−(P3 disjoint union K2), on five active vertices, and excludes its own anchor's common point. If A's common point y were outside both anchors, then the only A-rooted members would be A and its three two-point members: any singleton must intersect the two-point member supported on the other two A-points, hence contain y, which its singleton link forbids. Thus y in B and likewise B's common point in A. Consequently both rooted blocks consist of their anchor and eight one-point members. The three cross triples attached to A give a triangle in the exact-singleton link at y relative to B. Any B-rooted member other than B through y corresponds to an edge disjoint from this triangle in the five active vertices, so there is at most one. At least seven B-rooted members therefore avoid y.

Make graph G of the distinct outside pairs in the eight A-rooted triples other than A. It has at least three edges because pair degree is at most three. Every B-rooted member avoiding y has one of the two other B points and an outside pair Q that is a vertex cover of G: if some A-rooted triple {x} plus e misses Q, then that triple, the B-rooted member, and (A minus {x}) plus {y} are three disjoint members. Its outside pair can expand through at most two anchors, thus there are at least four distinct 2-element covers Q. If G contains disjoint edges ab and cd, all four covers must be the four cross pairs. The only 2-element vertex covers of these four pairs are ab and cd, so these would be the only edges of G, impossible. Otherwise G is a triangle or star (the elementary classification of pairwise-intersecting 2-sets); the triangle has only three 2-element covers. The star has at least three leaves, so all covers contain the center, which would be incident to all eight A-rooted and all seven B-rooted triples avoiding y, contradicting the degree bound 10. This proves the cutoff.

The only imported computer result in this proof is Lemma 5.9's bound for a 10-vs-at-least-9 rooted pair. In particular it does not depend on the finite route-cover table from Lemma 5.7. No proof that 33 is attainable; 33 is an upper bound.

## 4. Optional strengthening when matching <=2 and degree <=9

Under both conditions the matching-constrained size is at most 31. If it were 32, choose any disjoint anchors. Since N(T)<=23, both rooted sizes are 32 minus the other anchor's N, hence at least 9; a 10-rooted block is impossible using imported cross bound <=6 as above. Thus both rooted sizes are 9, and both anchors have N=23. Every triple occurs in a disjoint pair, since N<=23<32; consequently every triple has N=23.

The maximum in 1+s+min(21,24-2s) is attained only at s=1 or 2. For s=1 no pair of T appears in two other triples. If s=2 and some pair of T occurs twice more, its endpoints have exact-singleton sizes at most 6 each (degree at most 9, excluding T and the two other members), and the third endpoint has size at most 7, giving N(T)<=1+2+6+6+7=22. Hence every pair in every triple has total codegree at most 2. All exact-singleton links now have maximum degree at most two and matching number at most two. Such a graph has at most six edges by the component path/cycle classification (each component has at most 3 times its matching number). Thus N(T)<=1+2+3*6=21, contradiction. This scoped result does not close the size-32 or -33 degree-10 branches.

One exact next route if warranted: determine a new cross bound at (9,9) and (8,10) by the author's rooted catalogs, then incorporate a global constraint from each possible new maximum matching rather than merely retuning their scalar 50-profile scripts. A route to an actual record is either a 40 construction needing at least six deletions from F0, or a contradiction across all 49 profiles including the high-degree branch. Neither has been established.
