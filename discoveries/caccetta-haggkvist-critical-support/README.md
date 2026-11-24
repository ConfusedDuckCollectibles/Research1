# Forcing short directed cycles at the one-third threshold

Imagine a network where every connection has a direction. The
Caccetta–Häggkvist triangle conjecture predicts that if every node points
to at least one third of the network, there must be a three-node
directed cycle. The unrestricted conjecture remains open.

Take a network with `n = 3r` nodes, at least `r` outgoing connections
per node, and at least `r−1` incoming connections per node. Call a node
critical if it has exactly `r` outgoing and exactly `r−1` incoming.
We prove that avoiding a directed triangle then requires **at least
seven critical nodes, with a directed cycle among them**. This requires
no assumption that the graph is regular.

The proof also finds positive weights on a subset of the critical
nodes. Every node has exactly as much weighted incoming support as
weighted nonadjacent support. The proof excludes every such support
of size six or less. The acyclic case works without a finite census.

This extends a known implication of the [September 2026 proof of
Seymour’s second-neighborhood theorem](https://github.com/openai/math/blob/main/lean/docs/173.md).
That theorem already handles the nonmultiple-of-three part of the
statement and supplies the degree-tight support method. The new
restriction is that supports must be source-free, cannot be tiny,
and must have the stated adjacency to outsiders. This is not an
improved unrestricted degree threshold or a proof of the conjecture.
See the [proof](PROOF.md) and [independent audit](independent/audit.md).
This note is not peer reviewed.

### Reproduce

```sh
python3 verify_cycles.py
python3 controls.py
bash independent/replay.sh
```

The first two commands use the Python standard library. The last also
uses a C++17 compiler and independently enumerates the six-vertex
orientations, then checks the relevant weighted kernels exactly.
