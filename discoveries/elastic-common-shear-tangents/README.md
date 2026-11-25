# What changes when you stiffen only one phase?

Two solids with the same shear stiffness but different compressibility hide their microstructure from the initial effective stiffness. Change just one phase's shear stiffness, keeping its shape, fraction and compressibility fixed, and some of that hidden structure becomes visible. This note gives the exact set of possible initial responses for three measurements: the bulk stiffness slope, the slope of coupling between compression and shear, and the curvature of orientation-averaged shear stiffness. Any two of the three determine the third up to the specified coupling norm and material constants, and the proof supplies geometries realizing the entire locus, including anisotropic examples in 2D and 3D.

When the phase with the lower bulk stiffness is stiffened in shear and the initial coupling slope is zero, a larger bulk gain necessarily means a larger second-order loss relative to the *linear* shear gain. This does not mean the shear stiffness itself falls. In the synthetic half-and-half example with baseline common shear 1, the shear curvature magnitude at zero coupling ranges exactly from 3/11 to 21/44, a factor of 7/4. A finite test excludes one pair of targets and supplies separate two-phase controls with finitely many interfaces; this test illustrates the sensitivity result and is not an improvement over the literature's static endpoint results.

- `PROOF.md`: assumptions, exact all-dimensional tangent formula and locus, construction, and finite example. Starts from the scope of a perfectly bonded linear elastic periodic model. No dynamic, nonlinear, fatigue, interface or sample manufacturing claims.
- `FROZEN.md`: target and finite benchmark recorded before constructing the proof.
- `SOURCES.md`: closest known results, including the **classical** confocal shape map and stronger known static edge bounds.
- `verify.py`: replay with standard Python 3. Checks an independent direct closed-form layer energy on four contrasts, exact rational finite margins, and two explicitly confocal positive-geometry checks. This checks arithmetic, sign and normalization, not the universal theorem by enumeration.
- `REPLAY.txt`, `MANIFEST.sha256`: expected results and checksums.

Run `python3 verify.py` and `sha256sum -c MANIFEST.sha256`. No data, external packages, or access required.

AI-assisted research note, not peer reviewed; no priority claim. [Independent proof audit and replay](independent/AUDIT.md).
