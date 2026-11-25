# No 52-point spectral set on the 13 × 13 × 13 grid

**Theorem.** There is no spectral subset of cardinality 52 in the additive group $\mathbb F_{13}^{3}$.

A set is spectral if you can choose as many Fourier waves as there are points, with the waves mutually orthogonal when restricted to those points. The Fuglede conjecture asks whether that condition is equivalent to being able to tile by translation. A set with 52 points cannot tile this 2,197-point grid by size alone, so any spectral example of size 52 would disprove this instance of the conjecture. The proof rules out every such example.

This module closes **52**. The [subsequent 78-point exclusion](../../breakthroughs/fuglede-f13/#the-78-point-case) completes the entire finite grid case; see the [complete ledger](../fuglede-f13-39/STATUS.md).

## Proof map

1. [Initial reductions](01-reductions/PROOF.md) bound line and plane intersections, prove mixed-moment identities through degree five, and organize the large uniform pencils through exact plane-count profiles.
2. [Four complete layers](02-layers/PROOF.md) are impossible: the Fourier constraints impose small differences on too many shifts of a pair of slices, and an exact finite search excludes such a pair.
3. [Nonlayer structure](03-nonlayer/PROOF.md) shows that a nine-uniform pencil on either side would occur on both, with a single explicit projected shadow.
4. [Nine-uniform pencils](04-both-k9/PROOF.md) are impossible: all 432 anchored height choices produce inconsistent exact moment systems. The remaining branch still has moments through degree seven.
5. [Final branch](05-endgame/PROOF.md): exact plane-count and projective incidence bounds force 13 collinear slices with at most three directions. Their centered second-moment matrices cannot add to zero, except in an all-parallel case that has too many difference directions.

The proof uses integer and finite-field arithmetic. The [replay guide](REPLAY.md) lists the source, expected output and independent checks. No floating-point threshold or external dataset is needed. Independently replayed and audited within this project; not peer reviewed.

## Relation to prior work

[Fallon, Mayeli and Villano](https://arxiv.org/abs/1902.02936) established the finite-field geometric setup and verified the smaller primes 5 and 7. [Ferguson and Iosevich](https://arxiv.org/html/2610.30712v1) proved the prime-11 case. Their projection, moment and incidence methods, notably Proposition 10.6, Lemma 10.7 and Lemma 13.2, are close precedents. [Qihang Wang](https://arxiv.org/html/2609.20495v1) already gives the planar four-special-direction obstruction at mass 52; step 1 applies the same profile argument to the constrained weighted projection. The conclusion here is three-dimensional, not that planar theorem or the already-known arithmetic obstructions at sizes 65 and 91. No unrestricted priority claim is made.
