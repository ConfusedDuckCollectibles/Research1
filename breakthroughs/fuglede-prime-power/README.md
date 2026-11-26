# Fuglede on a circle: the prime-power case

Can a shape that supports a complete set of mutually orthogonal waves also cover space exactly by translation? This note settles a specified infinite class of the one-dimensional Fuglede problem and describes the wave frequencies in that class. It does **not** settle one-dimensional Fuglede without restrictions.

## What the proof establishes

- **Finite cyclic groups, any ambient order:** For a set with a prime-power number of points, being spectral is equivalent to tiling by translations. In fact, a spectral pair projects injectively to the cyclic factor for that prime. All sets sharing the spectrum admit the same explicit tiling complement.
- **Integer sets:** For each finite tile of prime-power size, the proof characterizes and counts all spectra containing frequency zero. The possible least denominators of its frequencies agree exactly with the possible fundamental periods of its tilings.
- **The real line:** A bounded spectral set of measure one tiles if it has a spectrum with a prime-power period. This is a spectral-to-tiling result only. Arbitrary periods are not covered.
- **The main reduction:** An unused odd prime-power factor can be deleted from a spectral pair when that prime divides neither the remaining group order nor the size of the pair. The matrix form applies to any root-of-unity Hadamard matrix with the stated twist.

For example, the four-point integer set `{0,3,10,13}` has exactly five spectra containing zero, taking frequencies modulo 1. Their least denominators are 4 or 20. Its tiling complements have exactly 4 or 20 as their possible fundamental periods.

## A separate sharp counterexample

The [companion proof](MAXIMAL-ORTHOGONAL.md) answers a different question: can a set of mutually orthogonal waves be impossible to extend even though it is not a full basis? For every number of points at least four, some cyclic tile has such a maximal incomplete family. For tile sizes one, two, and three, this cannot happen. This corrects a proposed generalization in the literature, not Fuglede itself.

## Follow-up: twice-prime sizes

The [twice-prime proof](twice-prime/) extends the one-dimensional results to all finite cyclic orders at sizes twice an odd prime, and to the matching spectral periods on the real line. The two open-ended families are distinct; neither implies the unrestricted conjecture.

## Check the result

Read the [proof](PROOF.md), then run the optional finite replays with Node.js:

```sh
node replay.js
node independent-exact.js
node maximal-exact.js
```

The proofs are algebraic and do not depend on these runs. The first replay checks ten examples and seven matrix-twist families. The second checks all 64,191 sets containing zero in nine small cyclic cases, recovering 9,226 normalized spectral pairs and confirming the exact phase criterion on these inputs. The companion check exhausts the potential extra frequencies in its examples for sizes 4 through 60.

## Prior work and status

The proof extends [Tan–Zhang’s odd-prime deletion](https://arxiv.org/html/2610.02218v1#S2) from an odd prime to its arbitrary powers. The two-prime step follows the exact positive-mask identity in [Malikiosis, Theorem 6.1](https://doi.org/10.1017/fms.2022.14). It uses the finite-to-real framework of [Dutkay–Lai](https://arxiv.org/abs/1301.0814), [Fu–Song’s rationality theorem](https://arxiv.org/html/2607.19858v1), and the tiling invariants of [Hu, Wang, Wang and Zhu](https://arxiv.org/html/2611.22849v1). See the proof for the division of work.

AI-assisted research note, with independent computational cross-checks. Not externally peer reviewed. No priority claim.
