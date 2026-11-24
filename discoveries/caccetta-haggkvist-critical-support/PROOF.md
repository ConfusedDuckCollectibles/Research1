# A weighted obstruction for one-unit semidegree loss in the directed-triangle problem

Research note, 24 November 2026. Uses the September 2026
Seymour second-neighborhood proof as a stated dependency. Not a resolution of the unrestricted
Caccetta–Häggkvist triangle conjecture; no priority claim.

## Theorem

Let D be a nonempty finite oriented graph: no loops, no multiple
arcs, and no two opposite arcs. Let n=|V(D)| and r=ceil(n/3). Suppose
every vertex has outdegree ≥r and indegree ≥r−1. If n is not divisible
by 3, D contains a directed triangle. If n=3r, put
```
L = {v : d⁺(v)=r and d⁻(v)=r−1}.
```
Then D contains a directed triangle if L induces an acyclic digraph,
or if |L|≤6. The latter includes the more concrete condition that all
but at most six vertices have indegree at least r. No outregularity
or strong connectivity assumption is needed. In a general loopless
simple digraph our conclusion is a cycle of length at most 3:
if both opposite arcs occur, there is already a 2-cycle. In that
setting do not assert an oriented triangle itself.

The structural statement behind the theorem: under a hypothetical
counterexample there is a nonzero nonnegative weight p with support
S⊆L, |S|≥7, such that D[S] has minimum indegree at least one; each
v∈S reaches **every** nonneighbor (other than itself) at directed
distance two; and at every vertex j, including outside S, the
p-mass of N⁻(j) equals that of the nonneighbors of j. Every vertex
outside S has both a support in-neighbor and a support nonneighbor.
The acyclic conclusion is entirely analytic; the seven-vertex
support threshold uses a tiny exact finite check only at order 6.

## 1. Dual support of the ordinary second-neighborhood matrix

Write U(v) for all vertices other than v that are adjacent to v in
neither direction. In a directed-triangle-free oriented graph,
the directed distance-two set N²(v) is a subset of U(v).
Define the integer matrix B by
```
B[v,j] = [j in N⁺(v)] − [j in N²(v)].
```
The vertex-weighted second-neighborhood theorem says that
for any strictly positive x there is a row (Bx)_v≤0.
It follows from the unweighted theorem by independent blow-up:
replace vertex j by x_j independent clones for positive integer x_j,
complete each directed bipartite arc class; directed distance-one
and distance-two sets have the intended sizes. There is no
length-two return to a clone class because the original graph
has no opposite arcs. Strict inequalities for arbitrary real
positive x would persist under rational approximation and
common scaling, so the integer version suffices. This weighted
equivalence is prior art (Seacrest; also [OpenAI §5][snc]).

For the uniform vector 1,
```
(B 1)_v = d⁺(v)−|N²(v)| ≥ 2d⁺(v)+d⁻(v)−(n−1) ≥ 3r−n.
```
If n<3r, strict positivity contradicts the theorem directly.
Now n=3r. Let J be the rows with (B1)_v=0; necessarily J⊆L
and N²(v)=U(v) on J. At every remaining row the slack is positive.
If there were a real vector z with (Bz)_v>0 for all v∈J,
then x=1+εz would remain strictly positive and give Bx>0 in
every row for sufficiently small positive ε. This is impossible.
By the real alternative of Gordan (or separation of the row
convex hull from zero), some nonzero nonnegative combination of
the rows in J is zero. Write its coefficients as p, extended by
zero elsewhere, and let S={v:p_v>0}. This set is nonempty.
Because each row indexed by S uses N²(v)=U(v), the columnwise
identities pᵀB=0 say precisely
```
p(N⁻(j)) = p(U(j))    for every j.                             (1)
```
No claim that this weight is positive on all of D is needed.
All its coordinates on S are positive.

## 2. No common dominators and no source in the support

If an outsider receives no arc from S, (1) at that outsider
says it has no missing adjacency to S either; it must send
to all S. Let A consist of such outsiders. No vertex of S
sends to A. No other outsider u sends to A, either, as u
receives an arc from some s in S and vertices of A send to s:
an arc u→a with a in A would give s→u→a→s. Hence any
vertex of A has indegree at most |A|−1, while |A|≤r−1
because A sends to every support vertex. This contradicts
the global indegree floor r−1 if A is nonempty.
Thus each outsider has a support in-neighbor, and its
positive in-mass equals its missing mass, so it has a
support nonneighbor as well.

If s is a source in D[S], (1) at s says that it sends to
every other support vertex. If z→s from outside S, no
vertex of S can send to z: an arc from s would make an
opposite pair and an arc from another u would make
s→u→z→s. This was excluded above. Thus d⁻(s)=0, so r−1=0.
The case r=1, n=3 already cannot be a counterexample:
one outgoing arc per vertex and no opposite arcs force
the directed triangle. Consequently D[S] has no source
and hence has a directed cycle. In particular D[L]
is not acyclic.

## 3. Support orders 4, 5 and 6

Since D[S] is source-free and triangle-free, |S|≥4.
A shortest directed cycle is induced: a chord in either
direction shortens it. If |S|=4, the support is that induced
directed 4-cycle. Equation (1) at its vertices shows that
the weight of each predecessor equals the weight of the
opposite vertex, hence all four weights are equal.

No source-free support of order 5 exists. A shortest 5-cycle
would be induced, and summing the balances inside it would
give total p-mass equal to twice itself. Otherwise let
c₀→c₁→c₂→c₃→c₀ be an induced cycle and let x be the fifth
vertex. Summing the four balances at the cycle vertices
cancels their internal terms, so x has as many outgoing arcs
as nonedges to the cycle. Its indegree inside S is positive.
If x has no outgoing arc, it receives from all four and
violates its own balance. Otherwise it sends to exactly
one c_j, misses one, and receives from two. To avoid a
triangle the missing vertex must be c_(j+1), since that
successor cannot send to x. The other two send to x. The
balances at c_(j+2) and c_(j+3) show p(c_(j+1))=p(c_(j+2)).
The balance at x would then equate
p(c_(j+2))+p(c_(j+3)) with p(c_(j+1)), impossible.

At order 6, the supplied standard-library certificate
`verify_cycles.py` enumerates every possibility by planting
a shortest induced cycle of length 4, 5 or 6 in the first
coordinates and enumerating all remaining pairs. All
source-free, triangle-free cases are covered. The counts
before identifying orbits are respectively 2757, 91, 1.
For each adjacency matrix it verifies exact rational kernel
existence from the principal balance equations on S. The
two isomorphism classes with a strictly positive kernel are
below; their weights are unique up to scale. Encode each
unordered pair (i,j), i<j, in lexicographic order as
0=j→i, 1=nonedge, 2=i→j.
```
code             p                  nested external in-neighborhoods
000120101201010  (3,1,1,1,3,3)      N⁻_ext(2) ⊆ N⁻_ext(1)
000010122012022  (4,1,1,1,1,4)      N⁻_ext(1) ⊆ N⁻_ext(0)
```
For each possible outsider j its adjacency type partitions
S into I (sending to j), O (receiving from j), and M (missing).
Equation (1) requires p(I)=p(M). Triangle-freeness requires
no arc in S from O to I. Enumerating just 3^6 types gives
the claimed two nestings, with respectively 15 and 27
feasible local types. But the first nesting has internal
indegrees 1<2 and the second has 1<4. The total indegrees
of all support vertices are the same r−1. In either case
the vertex on the left of the nesting needs *more*
external in-neighbors, a contradiction.

## 4. Excluding a directed four-cycle support

This step uses the independent proof auditor's shorter attachment
argument. It needs only exact degree counts, not a numerical
general triangle bound or the distance-two saturation.

Suppose S={s₀,s₁,s₂,s₃} is the uniformly weighted directed cycle,
indexed cyclically. Its members have outdegree r and indegree r−1
in D, and both degrees one in S. Every outsider has equally many
I and M by (1), and I is nonempty by §2. There are just two types:
- Type 1: I={s_i}, M={s_(i−1)}, O={s_(i+1),s_(i+2)}.
  Indeed the predecessor of the unique incoming vertex
  cannot be in O without a directed triangle. Let X_i be
  this class, of size x_i.
- Type 2: |I|=|M|=2 and O=∅. Let T be their union.

Indegree r−1 at s_j gives
```
1 + x_(j−1) + x_(j−2) = r−1.
```
Hence x₀=x₂=a, x₁=x₃=b, and a+b=r−2. Counting outsiders
now gives |T|=3r−4−2(a+b)=r>0.

Every vertex v of T has no arc to S and has at most r−1
possible outneighbors in T. Thus it must send to some X_i.
An arc from v to X_i is possible only if the two incoming
support vertices of v are exactly {s_(i−1),s_i}; otherwise
v receives from a support vertex in the outgoing pair
of X_i and the three vertices form a directed triangle.
In particular every vertex of T has a consecutive incoming
support pair and can send to only the corresponding X_i
outside T. Each pair includes one even-indexed and one
odd-indexed cycle vertex.

The even-indexed support vertices each send to one other
support vertex and to exactly a vertices in the corresponding
X_i. Their remaining outgoing arcs (r−1−a each) all go to T.
Counting these arcs in T gives r=2(r−1−a). Similarly
r=2(r−1−b), so a=b=r/2−1 (if r is odd, this is already
impossible). Consequently every vertex in D[T] has outdegree
at least r−a=r/2+1, since at most a outgoing arcs can
leave T. But an oriented graph on r vertices has average
outdegree at most (r−1)/2. Contradiction.

This rules out every weight support of order at most six.
Together with §§1–3 it proves the theorem.

[snc]: https://github.com/openai/math/blob/main/preprints/A-proof-of-Seymours-second-neighborhood-conjecture-September-23-2026/paper.pdf
