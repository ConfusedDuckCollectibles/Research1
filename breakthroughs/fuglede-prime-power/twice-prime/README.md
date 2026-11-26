# Fuglede for twice-prime sizes

This is a follow-up to the [prime-power result](../). It treats a different infinite family: sizes with **twice an odd prime** points, such as 10, 14 or 22. No bound is placed on the size of the ambient finite cyclic group.

- In any finite cyclic group, a set of size twice an odd prime is spectral exactly when it tiles by translation. Moreover, all sets sharing the same spectrum can use one explicitly constructed tiling complement.
- For a bounded measurable real-line set of measure one, a spectrum of period twice an odd prime guarantees a tiling. The same complement works for all the finite fibers, so they assemble into one real-line tiling.
- For a size that is the product of two distinct odd primes, the corresponding finite result holds in odd cyclic groups. In even ambient order at odd pq, tiling still implies spectrality; the reverse direction and common complement are not settled in general.
- The proof also counts all spectra containing zero in the core cyclic group whose order uses just the same two primes.

The underlying argument is a proof, not extrapolation from small cases. Its key step is that a vanishing sum of exactly pq phases involving just those two primes consists entirely of p-cycles or entirely of q-cycles. This forces a digit grid and then the Fourier zeros needed for the common complement. The [full proof and exact scope](PROOF.md) spell out how it works.

Run the finite checks with Node.js: `node replay.js`. The script tests eight model pairs and enumerates all 358,281 six-point subsets containing zero in two small cyclic groups. These controls support the general argument but do not prove it.

## Prior work and status

[Lam–Leung](https://arxiv.org/abs/math/9511209) supply the nonnegative cycle decomposition. The product complement for the appropriate pure and mixed zeros follows the [Coven–Meyerowitz framework](https://doi.org/10.1017/fms.2022.14). [Bell, §10](https://arxiv.org/html/2611.00307v1#S10) already proved the matched zeros, common complement and real-line lift at size six; the advance here is the unbounded twice-prime family and the core spectrum count. The [prior prime-power note](../) supplies the unused odd-coordinate deletion. The [full proof](PROOF.md) distinguishes other overlapping cases.

This AI-assisted note passed an independent proof audit and exact checks. It has not been externally peer reviewed. No priority claim. The unrestricted one-dimensional conjecture remains open.
