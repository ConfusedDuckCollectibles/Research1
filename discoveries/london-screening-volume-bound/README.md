# How much superconducting material can a magnetic reading prove is there?

**In a spherical London model, a thin outer shell can give exactly the same exterior magnetic field as a fully superconducting ball.** We prove a sharp upper limit on screening for a fixed amount of superconducting material, even if its location inside the ball is unknown.

## The result

Suppose the specimen fits in a known spherical region, its local superconducting stiffness has an independently known upper limit, and its weak-field response follows the static, vortex-free London model. Measure its induced magnetic dipole for three perpendicular directions of a uniform field and average the three responses. At any fixed material budget, the strongest possible average comes from concentrating the material in an outer shell at the stiffness limit. The proof allows any measurable layout, including disconnected material and continuously varying stiffness. It supplies a sharp threshold: an average above the shell’s response rules out all smaller material fractions in this model.

For one exact example, a shell occupying **5.8808%** of the ball and a fully active ball with slightly lower stiffness produce the **same magnetic field at every exterior point**, under any uniform applied direction in the model. The shell reaches about 96.925% of the ideal perfectly screening sphere’s dipole. A calibrated field gradient distinguishes those two particular radial profiles; it does not recover arbitrary interiors.

## Why it matters

When a measurement is used to infer how much of a specimen is superconducting, a strong exterior signal can still be consistent with far less than a fully active sample. The bound says exactly how strong that ambiguity can be under the stated geometry and stiffness limit. A single-axis reading, hysteretic or trapped flux, a conducting path, or a sample inside an unmodeled pressure cell does not supply the conditions needed to apply the bound. No real specimen was fitted, and the note does not decide a claim about any specific superconducting material.

## Proof and replay

- [Sharp volume bound and exact field counterexample](PROOF.md)
- [Stability of the maximizing layout](STABILITY.md)
- [How a higher magnetic harmonic separates the radial pair](QUADRUPOLE.md)
- [Context, related results, and scope](SCOPE.md)

Run from this folder with Node.js. No packages or experimental datasets are required.

```sh
node replay.cjs
node replay-quadrupole.cjs
node replay-independent.cjs
```

Expected outputs: [dipole](replay-expected.json), [quadrupole](replay-quadrupole-expected.json), and [independent interval check](replay-independent-expected.json). The independent check uses constant-potential bounds on 1,024 radial slices and a separate solid-sphere power series, rather than the closed-form shell response. All three scripts use exact arithmetic for their certified bounds.

Research direction and public archive: James Bazan / ConfusedDuckCollectibles. AI-assisted derivation, code, and independent computational and proof checks within the project. Not externally peer reviewed.
