# When syndromes can reveal a coherent error

A usual odd repetition-code syndrome history can be exactly the same for coherent rotations and their Pauli approximation. Waiting for more ordinarily measured rounds does not help in the ideal model. A two-slot comparison changes what can be seen: physically keep or omit a small number of intermediate checks, randomize a calibrated sign flip in the second slot, and use only the check records. An arbitrary stochastic Pauli model with the same noise law under both settings gives no signal, even when its error rates and correlations are freely refit.

The exact limit is set by the smallest connected block left by the measured intermediate checks. If it contains r qubits, the first visible echo signal scales as the angle to power 2r. The ideal needed shot count has the matching inverse exponent. With fixed, invertible but full-support reported-bit noise, the required shot exponent doubles. To test an isolated coherent source on one interior site of the nearest-neighbor chain, withholding either neighboring check alone is insufficient; withholding exactly both is sufficient. Nothing here requires a final logical-qubit readout.

This is an ideal, one-noise-axis calibration model. Real echo controls can themselves cause setting-dependent errors, and a positive contrast does not exclude that alternative unless matched controls bound it. No hardware ran; the matched ideal commands and their readout assumptions are part of the theorem.

- `PROOF.md`: complete joint probability theorem, exact tight ideal and noisy shot scalings, local readout crossover, and calibration limits.
- `SOURCES.md`: the closest active-probe, repeated-code, coherent-path and hardware literature.
- `FROZEN.md`: bounded target and assumptions.
- `verify.py`: standard-library, exact Gaussian-rational two-slot amplitude propagation, independently checked against the character theorem. Tests heterogeneous unequal slots and complex logical input, representative echo twirl, zero and positive single-source tests, exact cancellation at r=1, 2 and 3, and one ideal finite-budget synthetic power control. Decimal output is only formatting; assertions are rational.
- `REPLAY.txt` and `MANIFEST.sha256`: expected replay and integrity.

Replay: `python3 verify.py`, then `sha256sum -c MANIFEST.sha256`.
