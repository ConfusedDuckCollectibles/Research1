# No 204-point spectral set on the 17 × 17 × 17 grid

A 204-block shape cannot tile this wraparound grid: 204 does not divide 4,913. We prove it cannot pass Fuglede’s wave test either. This rules out one complete candidate size, not the entire 17-grid theorem.

**Proof map**

1. **Deduction:** published finite-geometry bounds force any hypothetical shape and its spectrum into twelve parallel planes with 17 points in each.
2. **Exact checks + deduction:** a planar direction check and a small table of quadratic profiles let us choose perpendicular layer coordinates. Fourier orthogonality then makes all cross-layer slopes in an active slice share a single small palette.
3. **Exact check:** across all 1,225 normalized domain-and-palette pairs, a cell contains at most four shapes, and their height vectors have affine rank at most two.
4. **Deduction:** the rank restriction forces a polynomial determinant to be zero at too many field elements to be possible.

[Full proof](PROOF.md) · [Independent reduction audit](INDEPENDENT-AUDIT.md) · [Independent final atlas and rank audit](INDEPENDENT-FINAL-AUDIT.md)

Run with Node.js, no dependencies or downloads:

```sh
node verify-f17-204.cjs
node independent-planar17.cjs
node independent-planar12.cjs
node independent-atlas.cjs
```

The 12-point direction test is an additional check, not a premise of the final contradiction. JSON expected outputs are alongside the code. The mathematical argument and the exact searches were reproduced independently within this AI-assisted project. No external peer review or priority claim.

The sizes still not excluded by the checks in this repository are **34, 68, 102, 136, 153, and 170**. See the [earlier ledger and general cubic constraint](../fuglede-prime-cubic-moments/).
