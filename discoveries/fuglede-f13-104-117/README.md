# Two more exclusions for finite Fuglede at p = 13

No spectral subset of **F₁₃³** has size **104** or **117**. The proof rules
out every set of those sizes, not merely sets having prescribed symmetries
or slice patterns. This note treats only those cardinalities; the [subsequent full-grid proof](../../breakthroughs/fuglede-f13/) settles the finite instance.

[Proof](PROOF.md) · [Exact Python replay](planar.py) ·
[Second C++ replay](independent.cpp)

A finite planar calculation forces at least eight uniform normals on each
projective line for these two putative sizes. The projective incidence bound
then forces all moments of degree at most ten to vanish, reducing their
occupancy profiles to a contradiction. There is no floating point, solver,
external dataset, or partial search.

Run 'sh replay.sh' with Python 3.10+ and a C++17 compiler. The expected
results are in [the recorded output](replay.txt).

## Role in the full result

Every tile in F₁₃³ is spectral, and every spectral set of cardinality 1,
13, 169 or 2197 tiles. The standard reductions leave sizes 13m with
2≤m≤12; see [Aten et al.](https://arxiv.org/abs/1509.01090) and
[Ferguson–Iosevich, §§8 and 13](https://arxiv.org/html/2610.30712v1).
[Malikiosis, Theorem 1.4](https://doi.org/10.1007/s43670-023-00078-2)
excludes 130, 143 and 156 because each exceeds 169−12√13. His §6, Eq. (6.1)
also gives a stronger **conditional** cutoff
13(13−36/(√34+1)) ≈ 100.49 if the Fourier support is a three-fold blocking set.
For p=13 it already excludes 104 and 117 in that branch, leaving a parallel
full-layer alternative. The argument here is unconditional and does not
assume away that alternative. The [separate size-26 note](../fuglede-f13-26/)
rules out 26, and the [subsequent size-39 proof](../fuglede-f13-39/)
rules out 39. The [updated ledger](../fuglede-f13-39/STATUS.md) records the classical determinant exclusion at 65 and 91, and the later proofs at 52 and 78.

The primary model is Ferguson–Iosevich’s proof at p=11, specifically
Lemma 8.2, Lemma 8.5 and Proposition 13.1. This is an AI-assisted research
note, not peer reviewed. We make no priority claim.
