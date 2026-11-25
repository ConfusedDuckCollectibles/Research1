# Fuglede's conjecture on the 13 × 13 × 13 grid

**Theorem.** In the additive group $\mathbb F_{13}^{3}$, a nonempty set tiles by translation if and only if it is spectral.

A spectral set is one on which a set of mutually orthogonal Fourier waves forms a basis. The conjecture says that this wave property is equivalent to a geometric one: copies of the set can cover the grid by translation, without overlap. This note proves the equivalence on the entire 2,197-point grid. It does not settle Fuglede's conjecture in other finite spaces or in Euclidean space.

**What changed:** Earlier results left one possible counterexample size, 78. The new proof rules it out. It forces any hypothetical 78-point counterexample into six parallel slices on each side of the Fourier pairing. The last step shows that two ways of cutting those slices would make one of them a straight line, violating a proved line bound.

- [Complete size ledger and main proof](PROOF.md)
- [The 78-point case, step by step](#the-78-point-case)
- [Run the exact checks](REPLAY.md)
- [Independent proof-chain and final check audit](04-final/CHAIN-AUDIT.md)

## The 78-point case

1. [Neither member layered](01-neither-layer/PROOF.md): integer projection profiles and a weighted projective count rule this out.
2. [Exactly one member layered](02-one-layer/PROOF.md): a finite graph check rules out every permitted four-layer pattern with too few difference directions.
3. [Both members layered](03-linked-layers/PROOF.md): a line argument forces line cap five and shows that there is a perpendicular pair of layer normals. This gives linked graph coordinates. Across each active cut, one common palette of exactly five slopes governs all six points.
4. [Final contradiction](04-final/PROOF.md): for every permissible set of layer labels, the palette fixes the height difference of one fixed pair of labels. Two active cuts force an entire layer graph to be affine, which is impossible.

The proofs and scripts were developed with AI assistance and checked independently within this project using separate exact implementations. They have not been peer reviewed. No claim of publication priority or external confirmation is made.
