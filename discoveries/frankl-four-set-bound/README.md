# A finite upper bound for Frankl-complete four-set configurations

**Theorem.** Any 24 distinct four-element subsets of a ten-element set
form a Frankl-complete configuration. Thus **FC(4,10) ≤ 24**. Combined
with Mingchang Liu’s existing 19-block non-FC witness, **20 ≤ FC(4,10) ≤ 24**.

A configuration is *Frankl-complete* (FC) if every finite union-closed
family containing its blocks has an element of the configuration’s own
support in at least half its members. Adding blocks cannot destroy
Frankl-completeness. This is a finite configuration bound; it does not
settle the exact value of FC(4,10) or Frankl’s full conjecture.

This is an AI-assisted, computer-assisted note, not peer reviewed. No
claim of priority. The new computation is fully replayable below, given
the explicitly cited upstream classifications.

## Prior results used

- [Pulaj and Wood, *Local Configurations in Union-Closed Families*,
  Theorem 2.3](https://arxiv.org/html/2301.01331v2#S2): FC(4,8)=12.
  Their Theorems 2.6–2.7 give an earlier lower bound and the averaging method.
- [Liu, *Frankl-complete sunflowers and extremal families of four-sets*,
  V1](https://doi.org/10.5281/zenodo.22702735) and its
  [exact companion release](https://github.com/michaeliu4/frankl-complete-four-sets/releases/tag/V1):
  FC(4,9)=16 and the complete eight-vertex non-FC lists of sizes 9 and 10
  (52 and 13 types respectively). We use the pinned
  [catalogue](https://github.com/michaeliu4/frankl-complete-four-sets/blob/V1/fc49/data/catalog8.txt)
  (Git blob `dc493f56112d292588aefa8e79512899ae9cfc90`). The upstream
  [driver](https://github.com/michaeliu4/frankl-complete-four-sets/blob/V1/fc49/verify.py)
  validates these inputs; the present replay does not rerun Liu’s entire
  certificate suite. It independently gets the 11-block layer by extending
  all 13 ten-block seeds by one edge, checking every single-edge deletion
  against the ten-block layer. The 28 successful labeled-template
  extensions are covered by exactly the three listed 11-block types.
- Liu’s publicly released
  [19-block dual certificate](https://github.com/michaeliu4/frankl-complete-four-sets/blob/V1/finite_bounds/nonfc_n10_19blocks.json)
  provides the lower bound. It concerns abundance on the original support,
  not a counterexample to Frankl’s conjecture. Liu’s
  [V2](https://doi.org/10.5281/zenodo.22734469) treats the separate asymptotic
  problem. The public V1 materials state no upper bound of 24; the
  immediate averaging bound from FC(4,9)=16 is 26, since 26·6 > 10·15.

## Proof

Use the hypergraph language: points are vertices and the distinct
four-sets are edges. A *card* is obtained by deleting a vertex and all
edges through it. The *codegree* of two vertices is the number of edges
containing both.

Suppose there is a non-FC 24-edge hypergraph G on ten points. All its
subconfigurations are non-FC. Each nine-vertex card has at most 15 edges,
so all degrees are at least 9. Their sum is 96, so every degree is at most
15, and the set L of degree-9 vertices has size `ell >= 4`.

**Generate candidate 15-cards.** Each vertex of L leaves a 15-edge card
with maximum degree at most 15. In such a card each degree is at least 4
or there would be a forbidden eight-vertex restriction of size 12 or more.
Its minimum degree d is at most floor(60/9)=6. Remove a minimum-degree
vertex, putting its 15−d surviving edges on a labeled representative of
the complete 9-, 10-, or 11-block eight-vertex lists. Add the removed vertex
with all possible selections of d distinct incident triples.
`cards15.cpp` keeps a partial branch only if the endpoint degree intervals
can meet [d,15] and all eight-vertex restrictions met so far remain
possible. The latter criterion is hereditary under adding edges, since
FC cannot become non-FC by adding edges. This produces 242 labeled
candidates, an over-approximation; it does not purport to decide whether
every retained card is non-FC.

**Case 1: two degree-9 vertices have codegree c >= 3.** Their common
deletion has 24−9−9+c=6+c edges, hence 9–11 edges and c in {3,4,5}.
Both single-deletion cards have 15 edges. In `prepare24.py`, for
every possible common deletion, transport every admissible card link
to its fixed template using *all* isomorphisms, not merely one labeling.
Every isomorphism preserves the degree of each vertex in the common
deletion, so the code only needs to enumerate bijections within
matching degree classes. Use all pairs, with repeats allowed and
order immaterial because the two deleted vertices have the same role.
Combine their two links with the common deletion. The remaining c
hyperedges contain both distinguished points and correspond exactly to
c distinct ordinary pairs on the other eight points. `verify24.cpp`
tests every c-subset of these 28 pairs for the remaining roots without
any partial-search pruning. Each completed configuration with minimum
degree at least 9 has a forbidden eight-point restriction. For rejected
roots, either the known blocks already have such a restriction, or the
required remaining degrees exceed an individual capacity c or total
capacity 2c. These rejections cannot remove a valid completion.

**Case 2: all low–low codegrees are at most 2.** For every edge e let
`t_e = |e intersect L|`. Then `sum t_e = 9 ell` and
`sum C(t_e,2) <= ell(ell-1)`. If ell=5, the pointwise inequality
`C(t,2) >= t-1` gives a contradictory lower bound of 21 > 20. If
6 <= ell <= 9, use `C(t,2) >= 2t-3` to get
`18 ell-72 > ell(ell-1)`; ell=10 is impossible since the degree sum
would be 90. Only ell=4 remains. In that case the two sides of the first
inequality coincide at 12, so every edge meets L in one or two points
and every low pair has codegree exactly 2. The other six vertices
have degree 10. Delete one low vertex. The three remaining low
vertices have degree exactly 7 in that 15-card. Their mutual codegrees
are still 2 since no edge could contain three low vertices. Every
surviving edge meets the marked three vertices in one or two points.
The complete candidate search checks every such three-mark in every
card of maximum degree at most 10: none satisfy the two constraints.

Both cases are impossible, proving FC(4,10) ≤ 24. As a consequence,
for n >= 10, averaging on ten-point restrictions gives
`FC(4,n) <= 1 + floor(23*C(n,4)/210)`; in particular the resulting
upper bound at n=11 is 37. This averaging corollary is not an
asymptotic improvement on Liu’s theorem.

## Exact offline replay

Run `sh replay.sh` using Python 3 and GCC or Clang with C++20 and unsigned
128-bit integer support. No packages, network or solver are needed.
Checks reject disabled assertions. An ordinary run takes seconds and
builds and writes outputs only in a temporary directory.

Expected counts:
- 242 labeled 15-cards after 1,079,242 recursive states.
- Exceptional case: 2,728 marks, 2,423 failed edge intersection and 305
  failed codegree; zero surviving marks.
- Case 1: 24 overlap templates, 12,018 unordered link pairs. Root checks
  reject 4,009 for degrees and 7,509 for already forbidden restrictions.
  Among the remaining 500 roots, test 23,542,155 ordinary-pair subsets;
  41,091 pass the degree condition, zero pass all eight-vertex restrictions.

The small data file is redistributed from Liu V1 under its original MIT
license. See `LIU-LICENSE.txt`. The upstream manuscript is not included.
