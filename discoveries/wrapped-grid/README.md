# Least periods for power-of-two wrapped grids

On a wrapped rectangular grid, the opposite edges are joined. Let τ(w,n) be the largest number of grid points that can be selected without putting three on a wrapped line. At fixed width w = 2^a, the sequence τ(w,1), τ(w,2), … has **least period w²/2**, for every a ≥ 1. For example, width 8 has least period 32.

The [proof](PROOF.md) shows that no smaller height can attain 2w points. Michael Skotnica’s 2016 thesis already constructs 2w points at height w²/2 and shows that w²/2 is a period; the addition here is the proof of minimality. The earlier construction is reproduced in bit form and credited to the thesis.

## Reproduce

From this directory, with Python 3.10 or newer:

```sh
python3 verify.py
```

No packages, network, or input data are needed. The verifier checks all 3,152,140 triples in the credited constructions for a = 1 through 7. It also checks that a deliberately corrupted construction fails, and independently counts the short-direction families and their incidences on ten small tori. The theorem itself is proved analytically; the finite checks do not stand in for the proof.

## Scope and status

This is a specialized result about lines on discrete tori, not a solution to the ordinary square-grid no-three-in-line problem. The note and verifier were prepared with AI assistance and are not peer reviewed. A targeted search did not find a matching proof of minimality; no priority claim is made.

The [MIT license](LICENSE) applies to the original verifier code only.
