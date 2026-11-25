# When a plasma echo looks like a gravitational lens

Fast radio bursts can arrive as two copies of the same signal. A compact gravitational lens is one explanation. This note proves that an idealized cloud of plasma can imitate a particular point-mass echo pair across an entire octave: the time between the two strongest arrivals, their relative strength, and their relative wave phase all agree in geometric optics.

The construction counts **all three** rays rather than throwing one away. In the reference example, the third is below **0.6%** of the weaker principal ray in voltage (below **0.0036%** in intensity). It gets fainter as the core gets steeper. Neither the curvature nor the finite-frequency approximation error is bounded uniformly in that limit.

This is a model result, **not** a claim about any observed burst, how common such a screen is, or any dark-matter bound. It compares relative features, not a known absolute waveform. Exact calibrated wave intensity cannot agree on an open frequency interval, and the two models can also differ in where the images appear. With the geometry, delay, band, support and maximum column fixed, there is even a strictly positive minimum calibrated full-wave intensity mismatch. The proof does not say whether that gap is large enough for a telescope to see.

- [Proof and limitations](PROOF.md): smooth compact positive-column construction, complete stationary-ray count, and exact-wave distinction.
- `node replay.js`: exact rational checks and a separately labelled floating-point control. No data or packages required.
- [Fixed-band wave mismatch](FIXED-CLASS-GAP.md): qualitative uniform gap, with no instrument sensitivity bound.
- [Independent audit](AUDIT.md): second derivation and reproduction checks. Run `node independent.js` for the separate exact certificate.

Closest prior: [Beniamini and Kumar (2026)](https://arxiv.org/html/2609.31842v1) analyze a Gaussian plasma lens with a faint central ray and an approximately geometric delay, and [CHIME’s interferometry search](https://doi.org/10.1103/PhysRevD.106.043016) already discusses plasma multipath. The narrower result here is simultaneous exact matching of all the relative geometric-optics observables with a globally counted screen. This note has not been peer reviewed.
