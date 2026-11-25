# Finite Fuglede at 13: ruling out 39 points

The finite Fuglede problem asks when a set that can be tiled by translation also has an orthogonal basis of character waves, and vice versa. This note rules out one possible counterexample size: **no set of 39 points in F₁₃³ has such a basis**.

The proof is computer-assisted, using exact integer and finite-field arithmetic. The main idea is that a hypothetical example has a particular 39-point planar projection. Directions in that projection force certain low-degree moments to vanish, which forces a symmetry incompatible with its character matrix. There is no large incomplete search in the proof; the decisive lift check uses one seven-point slice.

[Full proof](proof.md) · [Independent audit](independent/AUDIT.md)

## Replay

- Primary: `sh verify39.sh` (C++17, Python 3.10+, Node.js); compare with [expected.txt](expected.txt).
- Separate checks, from `independent/`: `node finish.js`, `node planar.js`, `node skeleton.js`, and `node grids.js k8-fixture.txt`. See [the recorded output](independent/EXPECTED.json). The grid check validates the resulting supports; the primary generation and a separate exact profile recomputation check census coverage.

With this result and the separate [26-point](../fuglede-f13-26/) and [104- and 117-point](../fuglede-f13-104-117/) proofs, four candidate sizes are still unexcluded in this field: **52, 65, 78 and 91**. The full conjecture over F₁₃³ remains open. This AI-assisted note has not been peer reviewed. It makes no claim of priority.
