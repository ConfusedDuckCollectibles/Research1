# A 21-line pattern rules out 58 equally angled lines in dimension 18

How many lines can meet at the same angle in 18 dimensions? We know 57 are possible and 60 are impossible. This exact computer-assisted proof rules out 58 or 59 **when the set contains one of the specified 21-line Pasch-wing patterns**. Outside those 21 lines, the others may move freely. They need not lie in a fixed lattice or follow the same construction.

The nearest earlier lattice argument used those 21 lines plus four fixed lines from a Latin-square rectangle. The result here removes those four extra anchors. It does not improve the unrestricted bound of 57–59.

- [Proof and prior work](PROOF.md)
- `python3 replay.py` runs the primary exact replay (requires a C++17 compiler).
- [Independent audit](audit/AUDIT.md): `python3 audit/independent.py` needs only Python's standard library.
