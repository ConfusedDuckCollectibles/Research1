# Two local barriers for improving the Calderbank–Fishburn 236-cap

Let S be the classical Calderbank–Fishburn cap in the seven-dimensional affine space over the field with three elements. Two exact computational checks restrict how S could be changed to obtain a larger cap:

1. **Any improvement must delete at least 17 points of S.** More precisely, any cap containing a point outside S and deleting at most 16 points of S has size at most 227.
2. **Any improvement must change points whose affine span is the entire space.** No replacement confined to one proper affine subspace can increase the size.

These are restrictions relative to this particular cap (and its affine images), not a proof of the global optimum. They neither construct nor rule out a 237-point cap in general, and do not show that 17 deletions suffice.

## Replay

Compile a C++17 standard-library program and run it:

```sh
c++ -O3 -std=c++17 verify.cpp -o verify
./verify
```

There are no external libraries, input files, random choices, or network calls. The program reconstructs the known 236-point cap, checks all 27,730 pairs, and checks the full finite cases used below. It should reproduce [replay.txt](replay.txt). To print the line packings and coloring as well, run:

```sh
./verify --certificate
```

Each `plane` line contains the normal, layer label, old count, and allowed-pool size, followed by a sequence of triples of point IDs. A `color` line contains an outsider's point ID and its color (0–6). IDs encode the seven trits in base 3, with the first coordinate most significant. The verifier constructs these objects deterministically and checks their properties before reporting success.

## Definitions and baseline

A cap is a set containing no three distinct collinear points. In this space a distinct triple is collinear exactly when its coordinatewise sum is zero. Thus the third point on the line through distinct x and y is −x−y, with operations in the field.

We use S in [Holden's 2026 preprint](https://arxiv.org/abs/2611.19543), equation (1), which gives the classical [Calderbank–Fishburn construction](https://doi.org/10.1007/BF01388452) in the model of [Edel and Bierbrauer](https://www.yvesedel.de/Papers/smallCaps.pdf), Theorem 5. The code generates the classical support design from two orbits of the coordinate cycle (2 3 5 6 4), starting from {1,2,3} and {2,3,6}. It takes the signed supports and even-parity full supports on one outer layer, their complementary support design and the same full supports on the other outer layer, and signed axes in the middle. No precomputed baseline table is included.

## 1. Deletion bound

For an outsider x, each baseline pair whose third point is x must lose an endpoint before x can be inserted. These pairs form a matching: for fixed x and one endpoint u, the other endpoint is uniquely −x−u. The minimum outsider secant count is 10, attained at 24 outsiders. Thus adding even one outsider requires at least 10 deletions.

For two outsiders x and y, one baseline deletion is forced if their third point −x−y belongs to S. After accounting for it, discard all edges incident to the deleted point. The remaining edges of the two baseline secant matchings must be covered by baseline deletions. The residual graph has degree at most two, and cycles are even because the edges of the two matchings alternate. Every nontrivial component is a path or an even cycle; on a component with b vertices the minimum vertex-cover size is ⌊b/2⌋. This follows by taking alternate vertices for an upper bound and using a matching of that size for a lower bound. Summing over components and adding the forced deletion gives the **exact** cost of adding just x and y. The replay tests all 1,902,225 outsider pairs and finds a minimum of 16.

For up to 15 deletions, at most one outsider can be added, so a cap containing an outsider has size at most 236−10+1 = 227. For exactly 16 deletions, consider the graph of the 1,616 outsiders whose individual secant counts are at most 16. Join two when their exact pair cost is at most 16. Any set of outsiders that is simultaneously addable must be a clique in this graph: every pair must be compatible with the same deletion budget. This uses a necessary condition only. The replay recomputes all 888 edges and checks a proper 7-coloring, so a clique has at most seven vertices. A cap after 16 deletions therefore has size at most 236−16+7 = 227. With no outsiders there is no improvement. This proves the first claim.

## 2. Affine-span bound

It suffices to test affine hyperplanes, since every proper affine subspace lies in one. Fix a hyperplane H and freeze all points of S outside it. The lines with exactly one point in H and two frozen points are precisely those formed by one point from each of the other two parallel layers. Delete their third points from H to obtain the allowed pool P. A line cannot have exactly two points in H, so any cap in P gives an admissible replacement, and every admissible replacement is a cap in P.

For each hyperplane, put s = |P|−|S∩H|. A packing of s point-disjoint collinear triples in P proves that a cap in P has at most |P|−s = |S∩H| points, since at least one point must be omitted from each triple. When s = 0 the pool size already suffices. The replay normalizes each nonzero linear form by its first nonzero trit and tries all three fibers, covering all 3,279 hyperplanes. It checks the requisite packings in every case: s = 0, 1, 2, 13 occur 1,669, 1,220, 388, and 2 times, respectively. This proves the second claim.

## Prior art and status

Thackeray's [*The cap set problem: Up to dimension 7*](https://arxiv.org/abs/2206.09804), Theorem 4.2, gives the global upper bound 288; Proposition 4.1 uses cross-section midpoint constraints in a different argument. Holden's Table 1 records the 236-point lower bound. Neither bound changes here. The result is a local, reproducible obstruction to one family of searches, not a claim of a breakthrough.

This note and its verifier were prepared with AI assistance and are not peer reviewed. No claim of priority is made.
