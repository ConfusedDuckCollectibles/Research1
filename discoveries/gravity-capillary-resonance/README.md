# Why exact resonances alone cannot sustain a wave cascade in this model

Surface waves are governed by gravity at large scales and surface
tension at small ones. A common idealization keeps only interactions
whose frequencies match exactly. We prove that in a deep, periodic
rectangular basin with rational dimensionless surface tension and
rational squared aspect ratio, those exact matches cannot carry waves
from a finite starting set to arbitrarily short wavelengths.

The frequency matches split into arithmetic groups. Each group has
only finitely many modes, because its allowed radii map to integer
points on an elliptic curve. Even when an interaction involves several
groups, a separate conserved quantity keeps initially inactive groups
inactive in the exact-resonant polynomial model.

This does not rule out near-resonant transfers, real-water turbulence,
finite-depth effects, or a cascade from data with infinite support.
[Shao's earlier work](https://arxiv.org/abs/2301.00115) used the elliptic
curve method for pure capillary waves on a sphere; the result here
concerns mixed gravity and capillarity on a periodic plane.
See the [proof](proof.md) and [independent audit](audit.md). Not peer reviewed.

### Reproduce

```sh
python3 verify.py
python3 independent.py
```

Both scripts use only the Python standard library. They verify exact
resonances, detunings, and mapping arithmetic; finiteness itself comes
from the mathematical proof, not the search.

### Finite-depth companion

[At fixed algebraic normalized depth and tension](finite-depth/ADDENDUM.md), exact frequency matches balance separately at every radius. The note gives an explicit mixed triad at a tuned transcendental depth and shows how pure-gravity quartets can get exponentially close to resonance. It credits Kartashova's finite-depth approach and related normal-form work; it is not a PDE or priority claim.

[Separate finite-depth result](../dense-water-wave-resonances/): at fixed positive tension, directly coupled noncollinear cubic resonance depths are dense even though an algebraic depth at algebraic tension has no such exact resonance. The follow-up also quantifies near resonances without drawing full dynamics conclusions.
