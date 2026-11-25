# Three points for three-direction polygons

Take a collection of convex shapes drawn with sides in at most three fixed directions. Suppose every four shapes include three that overlap at a common point. **Three points are enough to touch every shape**, even if some pairs are disjoint. Three is best possible.

A classical result handles copies of one triangle that have only been moved or resized. The proof here allows both triangle orientations, parallelograms, and non-similar hexagons of arbitrary widths. The key step is to keep track of which groups of shapes overlap while reducing a difficult case to a line segment.

This settles the fixed-three-direction case of the planar Hadwiger–Debrunner `(4,3)` problem. It does not settle the unrestricted problem or improve its general nine-point upper bound. The note credits the classical triangle proof and overlap arguments it adapts. It is AI-assisted and has not been peer reviewed.

- [Theorem and proof](PROOF.md)
- [Exact replay](exact_replay.py), run with `python3 exact_replay.py` (standard library only). It checks the nine-hexagon example that needs three points, the segment reduction in a nonpairwise example, and a near-miss that fails the overlap premise. The finite checks illustrate the theorem; the proof gives the general upper bound.
