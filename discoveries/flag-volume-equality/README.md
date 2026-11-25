# When does the flag–volume bound become exact?

In a polytope, a *flag* is a chain that starts at a vertex and adds one dimension at a time, such as a vertex on an edge on a facet. Martin Winter found an upper bound on a volume product using the number of these chains. His paper asks exactly which shapes make the bound tight.

We proved a geometric test for equality. For every face, choose any linear measurement along that face. Equality holds exactly when each such measurement can be extended to the whole polytope without widening its range. For a centrally symmetric shape, this forces all vertices to lie on each chosen facet or on its opposite facet. It narrows the search to two-level polytopes and makes the test finite there.

For the entire family of Hansen polytopes, constructed from perfect graphs, the test shows that exactly the Hanner polytopes meet the bound. The full classification for arbitrary polytopes remains open. This is about equality in an **upper** bound; it is not a new proof of the symmetric Mahler lower bound.

- [Proof and prior work](PROOF.md)
- [What happens in low dimensions](DIMENSION-THREE.md), including the result for polytopes with just two more vertices than their dimension.
- [The simple case](SIMPLE.md): in any dimension, a simple polytope meets the bound exactly when it is an affine product of simplices.
- [Exact catalogue replay](CATALOGUE.md): every centrally symmetric equality case through dimension seven is Hanner, using the published two-level catalogue through dimension six. The dimension-eight check is reported by Winter, not independently replayed here.
