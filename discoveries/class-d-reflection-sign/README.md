# Three measurements to determine a superconducting reflection sign

A finite superconducting device can have a reflection sign that ordinary conductance does not reveal. We construct two physically allowed scatterers that give **exactly the same probability for every single-particle input and output at zero energy in both leads**, but have opposite reflection signs. No hidden channel is needed for the ambiguity.

The way around it is to inject two indistinguishable electrons in two orthogonal input modes, rotate the two outgoing orbitals with an ordinary calibrated mixer, and measure what comes back together. Three specific mixer settings suffice. In the same local counting experiment, **any two fixed settings can fail** even if every possible return count is known. The score uses all injections, including when one or both particles transmit to the other lead.

This identifies the sign of a finite reflection block, not whether a material has protected Majorana modes or whether it is topological in the bulk. It requires a coherent elastic model, matching of electron and hole orbitals, indistinguishable inputs and calibrated energy and charge-count readout. No device was tested.

[Proof and measurement assumptions](PROOF.md) · [Primary literature](SOURCES.md) · [Frozen scope](FROZEN.md)

Run the exact, offline standard-library check: `python3 -B verify.py | diff -u replay.txt -`. It verifies the two exact counterexamples, the count formula and the calibration example; it uses no third-party data or optimizer.
