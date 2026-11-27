# Weighted Spectral Chvátal: the cycle family

**Complete for the family below. Full Conjecture H is still open.** This Research1 work is not peer reviewed and makes no priority claim.

Imagine cards with stickers. You want the best-scoring pile in which every pair of cards shares a sticker. The obvious strategy is to pick one sticker and take all the cards that have it. Conjecture H asks whether a mathematical upper bound always matches the score of the best such strategy for weighted boxes of cards with a particular closure property. We prove equality for a family whose deletion pattern consists of one or more loops.

## Exact theorem

Take any simple 2-regular graph on r ≥ 3 anchors: a cycle, or a disjoint union of cycles. Choose any n ≥ 2r+1 old points, assigning a distinct old tip to each edge. Remove each triple made from the endpoints of an edge and its tip, give all surviving old triples weight 1, and add a weight-2 pair from each anchor to the same fresh marker. In the graph where two sets are adjacent when they are disjoint, the weighted Lovász theta is exactly the largest star weight, `binom(n−1,2)`. The same equality holds at every nonnegative weight ratio in the associated two-layer downset defined in [the proof](proof.md). This includes the smallest allowed n and disconnected anchor graphs.

The theorem does not cover arbitrary deletion patterns, arbitrary weights, or all downsets.

## How the proof fits together

- **All r ≥ 5:** [analytic proof for large n](proof.md), [sharper cutoff](cutoff-26.md), and [exact closure of the remaining 141 cases](completion.md). Of those 141, a rational bound certifies 133 and rational residual repair certifies eight. Floating-point routines only propose certificates; exact arithmetic accepts or rejects them.
- **Four anchors:** [independent exact square proof](square.md), using coupled rational repair and bounds, independent of Gemini’s earlier square certificate.
- **Three anchors:** [Gemini Carbon’s all-size triangle construction](https://share.gemini.google/uSwHOzDoj5t7), checked independently in [triangle-replay.py](triangle-replay.py). The checker tests the exact star equations and the base cases n=7–11, then verifies 46 polynomial minors for all n≥12 by exact forward differences. Its integer seed is from Carbon’s construction; the verifier was reconstructed independently.

The [range theorem](range.md) explains why the star equations can be solved, but is not a positivity proof by itself. Background: [Ellis, Filmus, and Friedgut, §5](https://arxiv.org/html/2609.28404v3#S5).

## Run the checks

Run the Node.js commands from this directory. The full actual-defect check can take several minutes and accepts optional start/end n for sharding. Expected outputs and summary counts are included.

```sh
node replay.cjs | diff - expected.json
node range-replay.cjs
node cutoff-replay.cjs
node actual-defect.cjs
node small-repair.cjs
node square-small.cjs
node square-actual.cjs
node square-typed.cjs
```

For the triangle, run `python triangle-replay.py` with SymPy (checked with SymPy 1.14). See each proof note for the distinction between analytic results, exact certificate checks, and diagnostic output.

## Credit

James Bazan directed the investigation and iterated with Google Gemini and OpenAI’s Dot. [Gemini Carbon](https://share.gemini.google/DVKV2JFxKpUc) developed the triangle and square proofs and found a closed-form bridge step for general cycles. [Gemini Neon](https://share.gemini.google/n5eDqTJABc29) shared the corrected square source and a separate range argument. Dot developed and checked the analytic large-size proof, the sharper cutoff, an independent range proof, and the exact small-case closure, including the independent square proof included here. No external human peer review is claimed.
