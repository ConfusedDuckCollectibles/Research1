# Can two thermal gratings rule out positive-path heat transport?

A heat grating is a pattern of alternating warm and cool stripes. In a crystal its amplitude can oscillate either because heat-carrying phonons move together or because they travel at finite speeds without enough collisions to lose their initial motion. An oscillation by itself does not settle which effect produced it. This note gives a test that does not require an oscillation or a sign change.

**Result.** At the same elapsed time, compare the energy amplitudes of two gratings on the same axis. Subtract from each the collisionless response predicted by the same independently known velocity distribution. If heat transport can be represented by bounded-speed carrier trajectories preserving that distribution, the two differences cannot vary independently: their ratio lies in a sharp interval. No scattering times are fitted. More generally, the whole set of possible responses at any number of same-time gratings is characterized exactly.

**What makes the constraint physical?** A carrier can change velocity along its trajectory. Its average velocity, which sets how far it moves, is less spread out in a precise convex sense than the distribution of instantaneous velocities. Gratings probe this common distribution of averaged velocities at different spatial scales. The converse shows that no additional same-time test follows from bounded positive paths and the fixed velocity distribution alone; stronger claims require stronger dynamical assumptions.

**Worked example.** A standard collective phonon model conserves momentum in normal collisions and dissipates it through a weaker resistive collision. In its continuous-angle gray Debye version, both sampled grating amplitudes are positive and both exceed their ballistic predictions, yet the two-point constraint is violated. An exact rational calculation certifies a normalized separation greater than 0.0013 and tolerates an absolute error of 0.0004 in each ballistic-relative amplitude when the dimensionless wave numbers and coefficients are exact. This is a model witness, not an estimate of instrument precision or a measurement on graphite.

**Do not infer the cause from a violation alone.** The result tests the specified positive-path null, which includes ordinary energy-conserving relaxation-time kinetics with matched preparation and readout. It does not rule out coherent phonon effects, calibration error or changes in how heat is deposited or measured. The two-grating rule is sharp as a *closed cone of residual directions* for a uniform projected velocity law; it is not the full set of attainable amplitudes. The full characterization is a fixed-time result for general stationary paths, not a fit by one specific collision model at all observation times. No existing material measurements are reinterpreted here.

## Read and replay

- [Proof and continuum witness](PROOF.md)
- [Exact converse, more than two gratings, and calibration requirements](CONVERSE.md)
- [Prior work and provenance](SOURCES.md)

Offline, with Python 3 and its standard library:

```sh
python3 verify.py | diff -u replay.txt -
python3 verify_converse.py | diff -u converse_replay.txt -
sha256sum -c manifest.sha256
sha256sum -c extension_manifest.sha256
```

The first check bounds the exact continuous-angle model without discretizing velocities. The second checks an atomic stationary periodic realization and an exact example where the short two-grating cone is insufficient. Neither downloads experimental data.
