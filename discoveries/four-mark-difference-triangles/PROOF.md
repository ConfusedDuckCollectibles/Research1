# Construction and verification

## 1. Simultaneous equality

A four-mark ruler is written as (0,a,a+s,a+s+b) with positive integer gaps a,s,b. Its six differences are a,b,s,s+a,s+b,s+a+b. Write S for the sum of the largest marks of n rulers in a difference triangle set, and U for the maximum of these marks. All 6n positive differences from the n rulers must be distinct.

Consequently U is at least 6n. When U=6n, the differences are exactly [1,6n]. For each ruler, the sum of its six differences is 4(a+s+b)−a−b. All 2n outer gaps a,b are distinct positive integers. Thus, for every four-mark difference triangle set even without assuming U=6n,

    4S >= sum(1,...,6n) + sum(1,...,2n)
       = 20n^2 + 4n.

Equality holds exactly when the positive differences comprise [1,6n] and the two outer gaps per ruler collectively comprise [1,2n]. In particular, in the equality problem the other four differences of each ruler lie in [2n+1,6n]. This is the equality criterion from the lower-bound argument of [1, Proposition 7].

## 2. Exact-cover search and independent check

Up to reversal of each individual ruler, take a<b (equality is forbidden by distinctness). Enumerate

    1 <= a < b <= 2n,
    2n+1 <= s <= 6n−a−b.

Each row is the six-value set {a,b,s,s+a,s+b,s+a+b}. It has six distinct entries because a<b and the two intervals do not meet. An exact cover of [1,6n] by these rows gives a witness. Conversely every simultaneously optimal set has this representation. The search can therefore test these rows directly after imposing the two equalities. No structural condition besides necessary saturation and orientation was imposed.

The original search uses bit masks and branches on an uncovered value in the fewest remaining rows. The seven witnesses are supplied as plain integer coordinates. The separate checker recomputes the actual six pairwise differences by enumerating pairs of marks, not by trusting the search test. It also recomputes the sum, all the outer gaps, and every positive difference in one nontrivial product. The following successful search runs recover the five original case files (random generator seeded; row shuffling can be standard-library-dependent). A cutoff is not an exclusion.

    g++ -O3 -std=c++17 search.cpp -o search
    ./search 16 3500000 1
    ./search 17 2500000 8
    ./search 18 3500000 1
    ./search 20 200000 3
    ./search 21 700000 7

An additional exact-cover search produced the size-19 file. It used a separate smallest-column dancing-links implementation with xorshift32 row ordering; the search itself is not needed for verification. Both independent distance replays check the supplied file directly, including the full intervals 1–114 and 1–38, scope 114 and total span 1824. A further search for size 22 retained ten compatible rows from shifted size-21 seeds and completed the cover around the pinned ruler (0,43,88,132). The supplied witness is validated directly without relying on this construction heuristic.

## 3. Infinite perfect sequence with smaller total span

This is an application of the Wild product recorded in [1, Theorem 5, Proposition 6]. Let A and B be perfect families of n and b four-mark rulers with total spans S_A and S_B. Let v=12b+1. Include B and the scaled rulers vX for X in A. For each X in A and Y in B, additionally take the following 12 rulers, recentered to start at zero:

    v*X_i + Y_{pi(i)},     i=0,1,2,3,

where pi ranges over the affine permutations of the four-element field. These permutations map each ordered pair of distinct inputs to every ordered pair of distinct outputs exactly once. Coordinates in one product ruler remain increasing since the added perturbations are at most 6b in absolute difference, less than v.

For each positive difference d of A, coming from one specific ruler and one pair of indices, the offsets in the 12b corresponding differences run over all nonzero integers in [-6b,6b] exactly once. The unperturbed scaled ruler supplies the zero offset. Thus each d supplies the interval [vd−6b,vd+6b]. These intervals, joined with the initial [1,6b] from B, partition [1,6(vn+b)] without gaps or repeated values.

For each X,Y, the sum of the endpoint perturbations over pi is zero: each Y mark occurs three times at each endpoint. Hence the total span is exactly v² S_A+S_B and there are vn+b rulers.

Use the checked b=22 seed, S_B=2442, on the right at every step. With j=1 for the seed, the recurrence gives

    n_j = (265^j−1)/12,
    S_j = 2442 * (265^(2j)−1)/(265^2−1)
        = (5+1/133)*n_j^2 + (1−22/133)*n_j.

Every family is perfect. The j=2 case has 5,852 rulers, 35,112 distinct differences and total span 171,491,892, all verified by replay. The formula for all j follows by the product proof. For nonseed j this does not claim the exact total-span lower bound, and the product itself is not new.

## 4. Prior literature and scope

1. M. Shehadeh, F. R. Kschischang, A. Y. Sukmadji and W. Kingsford, [*Higher-Order Staircase Codes*, IEEE Trans. Inf. Theory 71 (2025), 2517–2538](https://doi.org/10.1109/TIT.2025.3544168); [arXiv v5, §VI-C](https://arxiv.org/html/2312.13415v5#S6.SS3). Proposition 7 gives the total-span lower bound, the paper supplies equality examples for 6–15 rulers, and Proposition 8 obtains the infinite-sequence coefficient 5+1/91 from its size-15 seed.
2. H. Liu, T. Feng, X. Wang and M. Zhang, [*Perfect difference families, perfect systems of difference sets and their applications*, JCTA 223 (2026), 106199](https://doi.org/10.1016/j.jcta.2026.106199), [arXiv:2510.20446v1](https://arxiv.org/html/2510.20446v1). Their Theorem 6.4 proves **scope** 6n is possible for every n except 2,3. Their concluding §7 asks about simultaneous attainment of 5n²+n. These are different problems: for example, the published scope-only size-16, -17, -18, -19 and -20 constructions have, respectively, printed outer gaps 49, 68, 40, 43 and 60 bigger than 2n (App. A, Lemma 4.2 or App. B, followed by Lemma 1.2).
3. M. Shehadeh, W. Kingsford and F. Kschischang, [*New Difference Triangle Sets by a Field-Programmable Gate Array-Based Search Technique*](https://doi.org/10.1002/jcd.22009), and T. S. Orden, [*Six optimal difference triangle sets for k=4*, Zenodo (2026)](https://doi.org/10.5281/zenodo.24096220) / [author's detailed account](https://cyborgmathematics.com/results/difference-triangle-sets-k4), cover nearby constructions with **five** marks per ruler; the six 16–21 constructions in the latter are scope-only, not total-span equality. Here there are **four** marks per ruler.
4. P. Wild, [*Combining perfect systems of difference sets*](https://doi.org/10.1112/blms/18.2.127), Bull. London Math. Soc. 18 (1986), 127–131. See [1] for the published DTS formulation of the product and its total span.

We do not claim to have exhaustively checked all historical constructions for accidental additional equality cases, or that the new asymptotic coefficient is globally best.
