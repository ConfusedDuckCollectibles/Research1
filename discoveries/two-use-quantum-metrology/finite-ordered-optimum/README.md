# Optimal quantum sensing with repeated noisy calls

**How much can you learn about a phase if you can use a noisy quantum sensor repeatedly?** For ideal controls between calls, this note gives an exact way to compute the best sensitivity for any finite ordered sequence in one standard noise model. One number measures the information carried forward, and each call requires an optimization over two real controls. This is an alternative to general quantum-comb optimization, not a claim about computational complexity.

- **A policy that works at the end can backfire in the middle.** An exact three-call example beats every policy that uses the maximal-overlap encoding at the middle call. Keeping some phase information even when an excitation is lost can help the calls that remain.
- **Memory need not grow to reach the known asymptotic limit.** A constructive policy reaches the established leading sensitivity coefficient with three reusable qubits in total and no saved outcome history. It uses ideal reset and control between calls, not unitary-only control. The bound is local to a chosen operating phase, not a lab result.
- The result does not solve the one-ancilla question, cover all types of noise, or claim a better asymptotic limit than the one already known.

[Proof and prior work](PROOF.md) · [Independent finite proof audit](audit/finite/AUDIT.md) · [Independent fixed-memory proof audit](audit/bounded-register/AUDIT.md)

Run `python3 verify.py` with SymPy for exact checks. The audits have separate verifiers, also using SymPy. For an optional independent numerical comparison, `python3 compare_comb.py --fetch` reads a pinned published dataset by SHA-256 and compares three entries using SciPy. No third-party source code is fetched or executed.
