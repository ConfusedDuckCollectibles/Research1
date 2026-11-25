# A 26-point exclusion for finite Fuglede at p=13

There is no 26-point subset of F_13^3 whose characters contain an orthogonal basis. This rules out one of the candidate non-tiling spectral cardinalities for the three-dimensional field of order 13. This note treats only size 26. The [subsequent complete proof](../../breakthroughs/fuglede-f13/) settles the full finite grid.

The [proof](PROOF.md) combines exact projection checks with a geometric argument using **both** members of a hypothetical spectral pair. A promising collection of projections actually survives the counting equations. It cannot occur on the partner where a nonuniform normal forces an injective projection, and the remaining mask types are incompatible with a six-point planar cap.

From this directory, run:

    sh replay.sh

This uses a C++17 compiler and Python 3.10 or newer with no extra libraries. The full runs use exact modular/integer arithmetic with no solver, randomness, or search cutoff. This note and code were AI-assisted and have not been peer reviewed. No priority claim is made.

Related work: [the size-39 exclusion](../fuglede-f13-39/) and [the size-104 and 117 exclusions](../fuglede-f13-104-117/).
