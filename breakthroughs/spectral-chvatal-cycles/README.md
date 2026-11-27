# Weighted Spectral Chvátal: the large-size cycle family

**Scope: a proved special case of Conjecture H, not a proof of Conjecture H.** This Research1 note is not peer reviewed, and does not claim priority over other work.

Think of a box of sticker cards. A valid hand is one in which every pair of cards shares a sticker. The best obvious hand takes every card with one chosen sticker. Conjecture H asks whether a particular mathematical upper bound always matches that strategy for a broad class of weighted boxes. We prove it for a family built from loops, once the box has at least 32 original sticker types and enough points to fit the loops.

## Exact result

Take any simple 2-regular graph on r ≥ 5 anchors, with possibly disconnected cycles. On n ≥ max(32, 2r+1) old points, assign a distinct tip to each anchor edge. Delete the triple consisting of that edge and tip; give surviving triples weight 1. Add the r pairs joining each anchor to the same fresh marker, of weight 2. The weighted Lovász theta of the disjointness graph is exactly the largest star weight, `binom(n−1,2)`. The certificate also proves equality for all nonnegative weight ratios in the associated two-layer downset. In particular, for at least 16 anchors the result applies at every admissible size.

The proof separates the exact weighted-star equations from positivity. A projection and a norm bound construct a certificate for the **whole** matrix, including the interaction between the different parts. It is not a numerical extrapolation. The analytic proof reduces the small-size work for this family to 468 cases up to cycle partition and old-point count for (5 ≤ r ≤ 15 and 2r+1 ≤ n ≤ 31), before crediting separately settled cases. The small cases are not certified by the large-size proof, and unrelated weighted downsets remain outside its scope.

## Proof and replay

- [Full proof of the large-size theorem](proof.md)
- [Exact correction-equation proof](range.md), which by itself does not establish positivity
- [Conjecture H and context](https://arxiv.org/html/2609.28404v3#S5), by David Ellis, Yuval Filmus, and Ehud Friedgut

With Node.js installed, run these **from this directory**. No packages or network are needed.

```sh
node replay.cjs | diff - expected.json
node range-replay.cjs
```

The first command reproduces the exact rational cutoff inequality and smaller algebraic controls. The second checks the correction equations on connected and disconnected examples, and verifies an exact negative minor for one sparse correction to show why equation-solving alone is insufficient. The proofs, not sampled controls, establish the universal claims.

## Credit and evidence boundaries

James Bazan directed the investigation and iterated on the problem with Gemini and OpenAI's Dot. Gemini Carbon supplied the all-size three-anchor certificate (with calculations independently checked by Dot), explored the square, and developed a bridge formula for general cycles. [Triangle report](https://share.gemini.google/uSwHOzDoj5t7) · [Carbon’s cycle report](https://share.gemini.google/DVKV2JFxKpUc). Gemini Neon [included the corrected square verifier and reported its results](https://share.gemini.google/n5eDqTJABc29), and separately developed a range argument. The square verifier was inspected but not replayed end to end for this release; the square claim and other individual Gemini checks are not dependencies of the large-size theorem.

The Research1 range proof and analytic large-size positivity proof were developed with Dot and checked separately against their exact arithmetic replays. The original research note and its replay code are included here. The linked Gemini reports remain separate source material. No claim of independent human peer review is made.
