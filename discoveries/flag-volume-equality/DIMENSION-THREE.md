# Winter's non-symmetric equality in dimensions at most three

This companion to PROOF.md uses its Criterion (A): a full-dimensional
polytope is sharp for Winter's flag–volume bound iff the difference
body of each face is the corresponding linear section of the full
difference body. The equivalent width extension will be useful below.

A sharp polytope is edge-antipodal in the classical sense of Talata
and [Swanepoel, 2007, §2](https://arxiv.org/abs/math/0601638).
In fact any edge $e=[u,v]$ has $u-v\in\partial(P-P)$ by (A)
on $e$, and any supporting functional there gives parallel
supports through $u,v$. This uses the stronger face condition to
imply the old edge condition, not a new characterization of
edge-antipodal polytopes.

In dimension two, the edge criterion itself gives the needed shape
directly. Put an edge of length $\ell$ at height zero, with the
polygon above it. If its endpoint angles are $\alpha,\beta$,
then at all sufficiently small positive heights $t$ its
parallel chord has length $\ell-t(\cot\alpha+\cot\beta)$.
It cannot exceed $\ell$, hence $\sin(\alpha+\beta)\ge0$
and $\alpha+\beta\le\pi$. Summing this inequality over edges
and using the angle sum $(n-2)\pi$ gives $n\le4$.
If $n=4$, every adjacent angle sum is $\pi$, so both
pairs of opposite edges are parallel. Thus only triangles and
parallelograms can be sharp; both are faces of Hanner polytopes.

In dimension three, [Bezdek, Bisztriczky and Böröczky, 2005](https://library.slmath.org/books/Book52/files/08bezd.pdf),
Theorem and remark on p. 130, give at most eight vertices, only
an affine cube for eight, and state that all edge-antipodal
3-polytopes are antipodal. For the latter they refer to a
supplementary analysis without including it; see also
[Bisztriczky's later same assertion](https://wpsites.ucalgary.ca/tibor-bisztriczky/).
The four-vertex case is a tetrahedron. The other cases use the
stronger range extension criterion:

* **Five vertices.** More generally let a sharp $d$-polytope
  have $d+2$ vertices. Its unique affine relation (up to scale)
  can be written $\sum_{i\in A}a_i v_i=\sum_{j\in B}b_jv_j$
  with positive coefficients on two disjoint supports and both
  sums of coefficients 1. Any other vertex has coefficient zero.
  Each support has at least two members by extremality.
  For any $i\in A,j\in B$, the remaining $d$ vertices
  are affinely independent, since otherwise there would be
  a second relation. The affine height of $v_i,v_j$ from
  their hyperplane has the same nonzero sign by the relation,
  so the remaining vertices form a simplex facet.
  Set an affine functional on this facet equal to 0 for its
  $A$-vertices and 1 for its $B$-vertices, with arbitrary
  0/1 values on its zero-coefficient vertices. Any extension
  confined to $[0,1]$ must have
  $a_i f(v_i)-b_j f(v_j)=1-b_j$, so $a_i+b_j\ge1$.
  Taking minimum coefficients on both sides forces
  $|A|=|B|=2$ and every positive coefficient to be 1/2.
  Hence the four vertices of the circuit are a parallelogram,
  and all zero-coefficient vertices are successive pyramid
  apices: at each step the apex is affinely independent of
  the others by uniqueness of the relation. For $d=3$ this
  is a square pyramid, up to affine equivalence. Conversely,
  for $k\ge1$, the $k$-fold pyramid over a parallelogram
  is the join of that parallelogram and a $(k-1)$-simplex.
  This is a facet of the free sum of a 3-cube and a $k$-cross-polytope
  (take the join of respective facets). For $k=0$ the base
  itself is Hanner. Thus all the $d+2$-vertex examples
  described are sharp by heredity.
* **Six vertices.** [Schürmann and Swanepoel, 2006](https://msp.org/pjm/2006/228-2/pjm-v228-n2-p08-p.pdf),
  Theorem 5.7, puts them on two parallel planes with three
  per plane. No triple in an antipodal set has an interior
  point of a segment between the other two, so these are
  triangular facets $A,B$. Range extension in both
  directions in (A) forces $A-A=B-B$: all vertices lie on
  those levels, and a normal tilt changes only the constant
  on the other level. A triangle is determined by its
  difference hexagon up to translation and reflection
  through a point. Indeed its three unoriented edge vectors
  are the three opposite vertex pairs of that hexagon.
  Choosing signs to close a triangle fixes all the relative
  signs, since the three nonparallel edge vectors have
  exactly one linear relation. Rearrangements only give a
  translation or central reflection. Thus the polytope is
  a triangular prism or an affine octahedron.
* **Seven vertices.** Schürmann–Swanepoel Theorem 5.8 puts
  the vertices at the corners of a cube except that
  $(0,0,0),(0,0,1)$ are replaced by $(0,0,\lambda)$,
  $0\le\lambda\le1$. The triangle in $y=0$ faces the
  entire square in $y=1$. The functional $x+z$ has
  triangle width $2-\lambda$ and square width 2, failing
  whenever $\lambda>0$; at $\lambda=0$, use $x-z$,
  triangle width 1 and square width 2. Again these are the
  only two levels so a normal tilt cannot reduce the
  opposite width. No such polytope is sharp.

All five surviving three-dimensional types occur as faces of
Hanner polytopes and therefore are sharp by heredity. The
tetrahedron is a facet of the 4-cross-polytope; the square
pyramid is a facet of the free sum of the 3-cube and a segment;
the triangular prism is a facet of the 3-cross-polytope times
a segment; the remaining two are Hanner themselves. This covers
Conjecture 3(ii) in dimensions at most three, using the cited
classical assertion on edge-antipodal 3-polytopes. The separate
$d+2$-vertex argument does not depend on that assertion.
