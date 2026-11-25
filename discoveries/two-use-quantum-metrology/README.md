# A two-use advantage for quantum sensing with feedback

A noisy quantum sensor can lose an excitation before it finishes sensing a phase. With two calls available, should the sensors run in parallel, or should the second wait for information from the first?

In this two-qubit amplitude-damping model, **feedback strictly beats every parallel strategy** whenever both calls retain some signal and at least one has nonzero loss. Each strategy gets one call to each of the two specified channels, and the parallel comparison allows arbitrary entangled preparation and noiseless ancillas. A simple strategy checks whether the first excitation was lost, then adjusts the input to the second only on that branch. Every run uses both calls. Nothing is discarded from the result, and the environment is not measured.

At survival probabilities 64% and 81%, the replayed construction gives **6.70% more quantum Fisher information** than the exact parallel optimum. The proof covers all survival pairs in the stated model, including unequal losses. Prior work already observed the advantage numerically for equal losses and for a particular unequal-loss family, and conditional loss recycling is not new in general. The result here is an analytic proof for the full two-call domain.

This is a local precision result in an idealized model. It does not give the optimal feedback strategy, a lab result, a general sensor claim, or an asymptotic advantage.

- [Proof, exact examples, and sources](PROOF.md)
- [Exact replay](verify.py): run `python3 verify.py` with SymPy.
- [Independent proof audit](audit/AUDIT.md) and [second exact verifier](audit/verify.py): `python3 audit/verify.py` uses only the standard library. It also checks a fixed input for which reversing the order improves the specific feedback protocol; it does not prescribe a globally best order.
