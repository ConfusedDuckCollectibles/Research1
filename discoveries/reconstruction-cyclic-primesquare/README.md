# Reconstructing a ring graph after one edge is removed

Imagine receiving every version of a network with one vertex missing, with the labels removed. Can you recover the original network? This is the classic graph reconstruction problem.

This note proves recovery for a specific infinite family: start with a connected circulant on p² vertices, for any prime p, and remove one edge. The resulting graph is uniquely determined by its vertex-deletion cards, even when competing graphs have no circulant structure. The proof also covers a separate prime-power family with a stated restriction on allowed step lengths.

The proof uses both local card structure and Fourier information. In one exact 289-vertex example, two different edge removals give the same adjacency spectrum but nonisomorphic graphs. The cards distinguish them. The full graph reconstruction conjecture remains open.

- [Proof, hypotheses, and primary sources](PROOF.md)
- [Independent audit](independent/AUDIT.md)
- Replays: `node check.js` and `node independent/check.js`, with expected JSON outputs.

This note has not been peer reviewed.
