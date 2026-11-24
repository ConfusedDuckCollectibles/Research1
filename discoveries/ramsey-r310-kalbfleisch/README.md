# No 40-vertex (3,10) Ramsey graph contains the Kalbfleisch graph

A short structural reduction and three independently checked SAT refutations. AI-assisted research note; not peer reviewed. No priority claim.

A (3,10)-Ramsey graph is a simple graph with no triangle and no independent set of size 10. Let C be the Kalbfleisch graph on integers modulo 35, with edges of differences ±1, ±7, ±11, ±16. We prove that a (3,10)-Ramsey graph on 40 vertices cannot contain C, even as a subgraph. **This does not determine R(3,10).** Its published bounds remain 40 ≤ R(3,10) ≤ 41.

## Reduction

The graph C is 8-regular, triangle-free, and maximal triangle-free. These facts, including that every nonedge has a common neighbor, can be checked directly from the displayed differences. Hence a copy of C in any triangle-free graph is induced.

In a hypothetical 40-vertex (3,10)-Ramsey graph, the maximum degree is at most 9: the neighbors of every vertex are independent. The minimum degree is at least 4, because at least 36 nonneighbors of a vertex would contain an independent 9-set by R(3,9)=36, which together with the vertex would form an independent 10-set. Each vertex in C therefore has at most one neighbor in the five-vertex exterior.

Fill in exterior edges while keeping the exterior triangle-free. This cannot create a triangle meeting C, since such a triangle would have a core vertex adjacent to two exterior vertices. Adding edges never creates independent 10-sets. The resulting exterior is maximal triangle-free. On five vertices the only possibilities up to relabeling are K(1,4), K(2,3), and C5: a maximal triangle-free bipartite graph on five vertices is connected and complete bipartite, and if it is not bipartite a shortest odd cycle has length five, with no allowed chord. Now fill in missing core-to-exterior edges as long as no triangle is created. This preserves the degree bounds and the bound on the independence number.

## Exact finite encoding

For one exterior F, the 175 primary Boolean variables x(u,a) describe whether core vertex u sees exterior vertex a. The formula imposes:

1. At most one exterior neighbor per core vertex, and independence in C of the core neighbors of each exterior vertex.
2. For each independent set T in F and independent set S in C with |S|+|T|=10, at least one x(u,a) with u in S, a in T. The core has independence number 8, so only |T|≥2 contributes clauses.
3. For each exterior vertex a, its core degree is between 4−deg_F(a) and 9−deg_F(a), inclusive.
4. Each missing core-to-exterior edge has a common neighbor, by the saturation step. The clause for u,a allows x(u,a), or x(w,a) for some core neighbor w of u, or x(u,b) for some exterior neighbor b of a.
5. A fixed exterior vertex 0 has degree at most 3 within the labeled F (take a leaf in the star). Hence it has a core neighbor, and translation of C allows x(0,0) to be set to true.

The degree counters encode exact prefix thresholds. The independent source checks the graph clauses by re-enumeration and checks every counter gate against its Boolean truth table. The core counts relevant to the three formulas are 3,360 independent 8-sets, 13,760 independent 7-sets, and 22,995 independent 6-sets.

| Exterior | Variables incl. counters | Clauses | DRAT |
| --- | ---: | ---: | --- |
| K(1,4) | 1,486 | 104,457 | VERIFIED |
| K(2,3) | 1,379 | 33,039 | VERIFIED |
| C5 | 1,435 | 22,861 | VERIFIED |

The three refutations exclude every allowable saturated graph. The saturation and relabeling above send any proposed extension to one of these formulas, proving the claim.

## What this rules out

The published 39-vertex example contains the Kalbfleisch graph, so adding one vertex without changing its existing edges cannot produce a 40-vertex example. The theorem also applies to any other construction containing the same core. By the published uniqueness of the 35-vertex (3,9)-Ramsey graph, every set of 35 vertices in a hypothetical 40-vertex (3,10)-Ramsey graph must contain an independent 9-set. No five vertices can therefore meet every independent 9-set.

## Independent replay without large downloads

Use Python 3.10+, [Kissat v4.0.4](https://github.com/arminbiere/kissat/releases/tag/rel-4.0.4), and [drat-trim](https://github.com/marijnheule/drat-trim). Run in a writable copy of this directory:

    sh replay_c35.sh /absolute/path/to/kissat /absolute/path/to/drat-trim

This generates the three CNFs from the step set, independently audits all clauses, regenerates three ASCII DRAT refutations, and requires “s VERIFIED” from the checker for each. A timeout or missing proof is not a result. No third-party graph-coordinate file is needed. REFERENCE-DIGESTS.txt records the original CNF and proof hashes; another build or version can generate a different valid proof. The small text package intentionally does not include the large generated CNFs, refutations, or third-party source code.

## Sources and scope

- Geoffrey Exoo, [*On Two Classical Ramsey Numbers of the Form R(3,n)*](https://doi.org/10.1137/0402043), SIAM J. Discrete Math. 2(4), 488–490 (1989), constructs a 39-vertex witness. The [official graph artifact](https://cs.indstate.edu/ge/RAMSEY/r3.10.39) has an induced copy of C, verified separately by an explicit isomorphism. This result closes every one-vertex extension of that particular witness, and all other 40-vertex extensions containing the same core; it is not a search of all 39-vertex witnesses.
- Jan Goedgebeur and Stanisław P. Radziszowski, [*New Computational Upper Bounds for Ramsey Numbers R(3,k)*](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v20i1p30/pdf/), Electron. J. Combin. 20(1), #P30 (2013), Theorem 3 identifies C as the unique 35-vertex (3,9)-Ramsey graph. McKay's [graph archive](https://users.cecs.anu.edu.au/~bdm/data/ramsey.html) supplies its graph6 form.
- Vigleik Angeltveit, [*R(3,10)≤41*](https://doi.org/10.37236/12936), Electron. J. Combin. 32(4), #P4.30 (2025), proves the incumbent upper bound. The order-40 case is not settled here.
- The October 2026 [Unsolved Lab research report](https://unsolvedlab.com/problems/ramsey-r310), method M20, listed the induced-Kalbfleisch branch as open and described an earlier bounded attempt. Its known degree-4 and translation-matching exclusions are different scopes. No priority over unindexed or unpublished work is claimed.
