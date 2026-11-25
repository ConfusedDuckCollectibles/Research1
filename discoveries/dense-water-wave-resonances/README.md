# Near-perfect wave matches at finite depth

On a square periodic water surface, three waves can exchange energy through a cubic interaction when both their frequencies and their wavevectors match in the required way. We fix a positive surface tension and construct exact three-wave matches arbitrarily close to **any** positive depth. The interaction coefficient does not vanish. The proof also counts these matches within the constructed family.

The same family reveals a sharp distinction. Some depths have no exact matches between different wave radii at any finite order, yet have cubic mismatches so tiny that one standard mathematical simplification breaks down: there is no continuous cubic Hamiltonian on any finite-order Sobolev space that removes all of those terms through the usual coefficient equation. Such depths form a residual set, meaning a countable intersection of open dense sets after excluding the exact resonances, while having Hausdorff dimension zero. For almost every depth in the measure sense, the family instead obeys an inverse-power bound.

This is a result about a particular periodic model and a particular transformation. It does not establish an energy cascade, failure of the water-wave equations, or a laboratory effect. Earlier authors already found dense depth resonances for standing-wave harmonics; the result here requires exact wavevector matching with noncollinear modes and nonzero coupling.

- [Proof and scope](proof.md)
- [Primary sources](SOURCES.md)
- [Exact interval replay](verify.py): `python3 verify.py`
- [Independent proof audit and exact fixed-frequency counts](audit/AUDIT.md): `python3 audit/exact_replay.py`

### Fixed-depth rounding

[The follow-up](fixed-depth/ADDENDUM.md) holds the depth fixed and quantifies how close the same family comes to resonance. It proves a rational–irrational split in the scaled mismatches, including an exact positive floor on one side and smaller-than-inverse-power-one mismatches infinitely often on the other. It concerns this family only.
