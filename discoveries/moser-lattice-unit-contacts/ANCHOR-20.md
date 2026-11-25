# Sharp 72 at order 25 when an extremal 20-point subset is present

**Proposition (using the public October 2026 values u(22)=60 and u(24)=68).** A set of 25 distinct planar points containing a 20-point subset that determines 54 unit pairs can determine at most 72 unit pairs. The bound is attained. This is not an exact value of the unrestricted u(25); a hypothetical 73-pair set may avoid all extremal 20-point subsets.

## 1. The only possible 73-pair extension would be an opposite-ended gem

By the published u(24)=68 and averaging vertex deletions, u(25)≤73. Suppose a 73-pair set has a 20/54 subset H. Every point of the full graph has degree at least 5 (otherwise its deletion leaves 69 or more edges). Also, no triangle has all three vertices of degree 5, since deleting them removes exactly 12 edges and leaves 61 on 22 points.

Write S for the five points outside H. The number of edges meeting S is 19, so Σ_{s∈S}d(s)=19+e(S). Since all five degrees are at least 5, e(S)≥6. If e(S)=6, all five degrees are 5, so G[S] is triangle-free by the previous paragraph. A triangle-free graph on five vertices with six edges cannot have an odd cycle (the only possible one is a 5-cycle, using all five vertices, and any chord would form a triangle). It is therefore bipartite, and six edges force K_{2,3}; this is not a planar unit-distance graph because two distinct circles intersect at most twice.

By the known u(5)=7, we have e(S)=7. The only 5/7 unit graph is the gem (the join of a 4-vertex path with a universal center); this is the n=5 row of the AMP exact classification and is also elementary from the unit-circle description. The degree total in S is 26, so one point has full degree 6 and the other four have full degree 5. The degree-6 vertex belongs to all triangles by the deletion constraint and hence is the gem center (the three triangles of the gem have only this vertex in common). The two endpoints of the 4-vertex rim path consequently have exactly three unit neighbors in H. They are opposite on the circle of radius one about the gem center: successive rim steps are rotations by +60 or −60 degrees, and reversing a sign would retrace a vertex. Thus the two endpoints are at mutual distance 2.

## 2. Every embedding of H lies in a four-item normalized list

Let eta=(1+i sqrt(3))/2 and w=(5+i sqrt(11))/6. Put
    A=(1,2,eta,1+eta,2 eta),       B=(1,2,eta,1+eta).
On the indexed product A×B, use the Cartesian unit graphs of these two triangular patches plus the extra cross-edge between the indexed points (2,1) and (1,2), zero-based. Exact distances in the embedding (a,b)↦a+w b give precisely these 54 edges, with 20 distinct points. The order-20 row of Alexeev–Mixon–Parshall Table 1 states that the 20/54 unit graph is unique up to isomorphism, so this labels every possible H.

Here normalization is of **all geometric embeddings**, not just the published one. The Cartesian 4-cycles are unit rhombi and satisfy the opposite-vertex sum identity. For every edge in the A factor, its oriented displacement is therefore the same in consecutive B fibers, and then in all B fibers by connectivity. Integrating along an A path shows that any embedding of the product is f(a,b)=P_a+Q_b. The A patch is rigid up to Euclidean isometry: fix its triangle on labels (0,1,3) and then place triangles (0,2,3) and (2,3,4). Each last vertex is forced to the other unit-circle intersection, since the previous opposite triangle vertex is already occupied and all points are distinct. The first two triangles give the same rigidity for B. After a common isometry and adjusting translation constants, P_a=a; then Q_b has the form rho b or rho bar(b), |rho|=1.

Write q=2−eta; |q|²=3. The cross-edge equation is |−q+rho q'|=1, with q'=q in the direct case and q'=bar(q) in the reflected case. Hence rho q'/q has real part 5/6 and is on the unit circle, so equals w or bar(w); q/bar(q)=bar(eta). The four possibilities are thus
    f(a,b)=a+w b,       a+bar(w)b,
    f(a,b)=a+w bar(eta)bar(b),       a+bar(w)bar(eta)bar(b).
All four pass exact graph and injectivity checks. It is unnecessary to claim these four are mutually incongruent.

## 3. Exact tiny circumcenter gate and sharp lower control

For each listed embedding of H, take every unordered triple of its 20 points, find its exact circumcenter if noncollinear and keep those of radius one. Any exterior point adjacent at unit distance to three H-points must occur in this complete enumeration (three distinct points on a circle are automatically noncollinear). The Cramer check uses the (X,Y) coordinates and Q(sqrt(33)) rational-pair arithmetic from MOSER-CONTACTS.md, with denominators included. In each of the four embeddings all 1,140 triples were tried. Exactly 306 had radius one; 302 triples have their centers in H. There are exactly four distinct exterior centers, each with three H-contacts. None of the six pairs of possible exterior centers is at distance 2. This contradicts the two endpoints from §1.

The known lower bound 72 is attained even in this restricted class: use all 25 points {a+w b : a,b∈A}. The exact code asserts 72 unit pairs; the 20 obtained by restricting b to B give the exact extremal subset above. This is a simple replay of the published lower count, not a claim of new record construction.

**Replay:** run `node verify-anchor20.cjs`. It uses integer rational arithmetic and needs no external data. In addition to the exterior-center test, it verifies 20 distinct core points, the exact 54 edges, the labelled Cartesian graph plus cross-edge, every accepted circumcenter against the three original distances, and the 25/72 positive control.

**Precise references.** Alexeev–Mixon–Parshall, [arXiv:2412.11914v2](https://arxiv.org/html/2412.11914v2), Table 1 rows 5 and 20; §4 uses the underlying rhombus identities for a different purpose. Artificial Wasteland, [u(22)=60](https://doi.org/10.5281/zenodo.22468977), [u(24)=68](https://doi.org/10.5281/zenodo.23272696), source computations not externally refereed. Engel et al., [arXiv:2406.15317v3](https://arxiv.org/html/2406.15317v3), §2 Moser basis and Table 2 known 25/72. This argument is separate from the near-Moser classification and does not assume an arbitrary 25-point set lies in a Moser copy.
