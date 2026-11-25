# Independent audit: no spectral cardinalities 104 or 117 in F13 cubed

**Conclusion.** The two cardinality exclusions check. This is a two-case spectral-to-tiling reduction, not a classification of all spectral sets or full Fuglede over F13 cubed. The input from prior work is Ferguson–Iosevich's geometric and moment template, not an already-generalized p=13 theorem. An independently implemented, exact three-branch planar replay is included.

## Spectral bounds and transfer (Ferguson–Iosevich Lemma 8.2)

In the space V=F13^3, set p=13, let (E,A) be a putative spectral pair of common cardinality mp, and let U(E) consist of projective directions v for which all p plane fibers x·v=t contain exactly m points of E. D(A) is the set of nonzero projective secant directions of A. Orthogonality gives a Fourier zero for every secant. Because the integer count polynomial of degree ≤p−1 vanishes at a primitive p-th root, it is a scalar multiple of the p-th cyclotomic polynomial. Thus D(A)⊂U(E); unitary transposition proves D(E)⊂U(A).

For an affine line of direction d containing a points of E, the p+1 planes containing that line meet off it disjointly and have occupancy total mp+pa. Their normals sweep the projective line corresponding to d^⊥. A complete projective line in U(E) would give pa=m, impossible for 2≤m≤p−1. Any plane containing p+1 points of E realizes every secant direction of its own direction plane (p parallel lines in each direction). This would give a full projective line in U(A), so the plane cap is p on both partners. The m p points of A occupy p parallel planes with normals d, hence two have difference in d^⊥ and give a uniform plane through the line of E. Consequently a≤m, and the pencil bound mp+pa≤m+p² gives a≤p−m after taking the integer floor. This reasoning uses only orthogonality of the standard nondegenerate pairing, not an orthogonal direct-sum decomposition; isotropic directions are allowed.

For **m=8**, take any two-space P of secant directions. At least one of the p affine translates of P contains eight points of the partner. The initial line cap is five. Gate (8,5,5) says that any chosen eight points there give six or more directions in P, so every projective line has ≥6 points of U for each partner. Every affine pencil therefore has at least six 8-point planes and can have at most eight others with ≤13 points each. Thus 104+13a≤6·8+8·13=152, so a≤3. Reapply the slice step with gate (8,7,3) to get ≥8 points of U on every projective line.

For **m=9**, the initial line cap is four. Any parallel family has a slice of at least nine points, and gate (9,7,4) yields ≥8 directions of U on every projective line. Both cases give 8≤|U∩L|≤13 for every projective line L. Optional arithmetic check: at this stage m=8 and m=9 both improve to line cap two, since 104+13a≤8·8+6·13=142 and 117+13a≤8·9+6·13=150 respectively. The initial line cap p−m was already enough for the terminal proposition.

## Independent exact gates

Below each triple records (point count, maximal allowed number of directions, maximal allowed number of collinear points):

| Gate | Unique normalized feasible point sets |
|---|---:|
| (8,5,5) | 0 |
| (8,7,3) | 0 |
| (9,7,4) | 0 |
| (8,6,5), neighboring positive control | 72 |
| (8,7,4), neighboring positive control | 9 |

Cardinality exceeds line cap in each case, so any candidate includes at least one noncollinear triple. Any ordered such triple can be sent to the three fixed points (0,0),(1,2),(3,1) by an invertible affine map, preserving line occupancies and globally permuting directions. To exhaust the remaining points the checker uses increasing order, never revisiting a sampled point or skipping a possible subset. The candidate line cap is monitored when adding point i: for each direction, count previously selected points on that line through i, rejecting exactly when the count reaches the cap (because adding i exceeds it). Old lines remain unchanged. A second monotone branch constraint is the popcount of the **actual** union of secant directions, independently of the author enumeration over fixed *allowed* direction masks. Slopes are computed by exact modular multiplication, with a vertical sentinel. The 9 positive sets in the last row each use exactly 6 directions, so a search over 7-element allowed masks reports 9·8=72 set/mask pairs, not 72 unique sets; this does not affect any zero gate.

Run `c++ -O3 -std=c++20 gates.cpp -o /tmp/f13-lift-audit && /tmp/f13-lift-audit` without disabling C assertions. Independent executable runtime is under 0.1 s on the audited system. The log includes positive point-coordinate witnesses and per-level node counts.

## Degree-10 rigidity (Ferguson–Iosevich Lemma 8.5)

Write N=|U| and n_L=|U∩L|. There are 183 projective lines, all of length 14; each point belongs to 14 lines, and one line contains any two points. Consequently Σ_L n_L=14N and Σ_L n_L²=N²+13N. If a nonzero homogeneous polynomial of degree ≤10 vanishes on U, multiply it by a power of a nonzero linear form to make the degree 10 if needed. A line on which the restriction is a nonzero binary form contains at most ten projective roots, including the possible point at infinity. Otherwise its rational linear defining form divides the polynomial (formal polynomial restriction), so there are at most ten such lines. A nonfactor line contributes (n_L−8)(n_L−10)≤0, and a factor line contributes at most (13−8)(13−10)=15. Summing gives
`N^2−239N+14640 = Σ_L (n_L−8)(n_L−10) ≤ 150.`
The quadratic's integer minimum is 360 at N=119 or 120, excluding it. The possible N real minimum 1439/4 is also sufficient. No assumption that U equals the secant direction set is used.

For 1≤j≤10, the homogeneous polynomial M_j(v)=Σ_(x∈E)(v·x)^j vanishes on U because m Σ_t t^j=0 in F13. Rigidity makes each M_j formally zero. (Invertible multinomial coefficients recover all mixed moments if using the literal statement of Proposition 13.1, but their individual names are unnecessary here.) M_0=mp=0 in the field as well.

## Terminal affine plane profiles (Ferguson–Iosevich Proposition 13.1)

For **each** nonzero normal v, including those outside U, the integer plane occupancies n(t), 0≤n(t)≤13, have moment equations through degree ten over F13. The evaluation space of polynomials of degree ≤10 on all 13 field values is dimension 11; by finite field power sums the evaluations of 1 and t are in its annihilator, are independent, and span it. Thus n(t)≡at+b mod 13 on the full field, not just on occupied planes. If a≠0, residues 0–12 occur once, their least possible integer sum is 78, and only the zero residue may be raised by 13. The actual total must be 78 or 91, not 104 or 117. If a=0 with nonzero residue then actual counts all equal m. If a=b=0, counts are 0 or 13 with exactly m occupied planes.

Choose a secant line meeting E in r≥2 points. Its 14 incident planes are nonempty, so their occupancies are m or 13. If t of them have occupancy m, the pencil identity gives
`(13−m)t=13(14−m−r)`.
Since 13 is prime and m=8 or 9, t is a multiple of 13 and lies in [0,14]. If t=0 then r=14−m, strictly above the p−m line cap. If t=13 then r=1. Both contradict the chosen secant. This closes all cases even where the Fourier support fails to be three-fold blocking.

## Priority and scope

- [Ferguson–Iosevich, arXiv:2610.30712v1](https://arxiv.org/html/2610.30712v1), Lemmas 8.2 and 8.5 and Proposition 13.1, provides the general ingredients that are instantiated here. Its F11 results are not an F13 enumeration. Its normalization of planar searches at Lemma 8.7/10.1 is direct precedent.
- [Fallon–Mayeli–Villano, arXiv:1902.02936v2](https://arxiv.org/abs/1902.02936), proves p=5,7; their general restriction and parallel-plane ideas are antecedent.
- [Malikiosis, Theorem 1.4 and §6 Eq. (6.1)](https://link.springer.com/article/10.1007/s43670-023-00078-2) must both be cited. The unique unconditional cutoff there is p²−p√p+√p≈125.733 at p=13, which excludes 130/143/156 but neither 104 nor 117. The §6 **conditional** alternative already excludes both if the Fourier support is three-fold blocking, since the cutoff 13(13−36/(√34+1))≈100.488. It leaves the alternative of containment in m parallel planes. No complete prior closing these two exact cardinalities was found in the bounded primary search. The added arithmetic excludes all directions uniformly rather than assuming three-fold blocking; it closes the prior alternative too.
- Do not conflate these cardinalities with the separate prior size-26 note, the full finite-field Fuglede conjecture, or universal portability of these planar gates to other primes. The novelty here, insofar as established, is the exact p=13 gates and consequent unconditional two-case exclusion.
