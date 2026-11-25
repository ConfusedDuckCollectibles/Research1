# An equality criterion for Winter's flag–volume inequality

AI-assisted Research1 note, 25 November 2026. Not peer reviewed. The argument starts from Martin
Winter, [arXiv:2610.09990v1](https://arxiv.org/html/2610.09990v1),
Theorem 1 and Lemma 4. A flag is a maximal chain of nonempty faces,
one per dimension from vertex through the full polytope. Volumes for
a face and its intrinsic polar are taken within its own dimension
after translation; the product is affine invariant. Write
$$
 M(E)=\operatorname{vol}(E)\operatorname{vol}((E-E)^\circ),\qquad
 B(E)=\frac{\#\operatorname{flags}(E)}{((\dim E)!)^2}.
$$
Say that $E$ is *sharp* if $M(E)=B(E)$. For $P=-P\subset\mathbb R^d$,
this is exactly equality in Winter's (2), because
$(P-P)^\circ=\tfrac12P^\circ$. This is equality in an upper bound relative to the flag count. Having
the minimum flag count $2^d d!$, together with the symmetric Mahler
lower bound, implies sharpness. The converse is the separate question
studied here.

## 1. Exact criterion with no symmetry assumption

For a nonempty face $F$ set $L_F=\operatorname{span}(F-F)$. We prove
$$
 P\text{ sharp}\quad\Longleftrightarrow\quad
 (P-P)\cap L_F=F-F\quad\text{for every nonempty face }F. \tag{A}
$$
Equivalently, every affine functional on every face $F$ extends to
an affine functional on $P$ with the same minimum and maximum.
The criterion concerns every face, not just facets of $P$.

**Flag criterion.** For a flag $C=(F_0,\ldots,F_d)$, define
$$
 m_i(t)=\min_{v\in F_i}t(v),\qquad
 g_C(t)=(m_0(t)-m_1(t),\ldots,m_{d-1}(t)-m_d(t)),
$$
where $t$ is linear (an additive constant can be immaterially added).
On the open cone of strictly positive gaps, $g_C$ is a bijection
onto $(0,\infty)^d$. Indeed, set the value on $F_0$ to zero; once the
affine functional on $\operatorname{aff}F_{i-1}$ is fixed, all
extensions to $\operatorname{aff}F_i$ differ by a multiple of a
height $h$ vanishing on $F_{i-1}$ and positive at each vertex of
$F_i\setminus F_{i-1}$. The minimum over these outside vertices is
the minimum of finitely many affine functions of the height
coefficient with positive slopes. It is continuous, strictly
increasing from minus infinity to plus infinity. There is a unique
coefficient giving the desired strictly reduced $m_i$. This
works whether or not the minimizing vertex is unique. Every affine
functional differs from exactly one linear one by a constant.

Write $W_P(t)=\max_Pt-\min_Pt$ and
$\Delta_d=\{x_i\ge0,\,\sum_i x_i\le1\}$.
Winter's proof partitions $Q=(P-P)^\circ=\{W_P\le1\}$ into
finitely many closed vertex-order chambers. On each chamber $g_C$
is his matrix $D_v$; determinants and pulling triangulations give
$$
 M(P)=\frac1{d!}\sum_C \operatorname{vol}_d\big(g_C(Q)\big).
$$
Here images of different chambers have disjoint full-dimensional
interiors (as in Winter, Lemma 4), or directly by the uniqueness
just shown for positive gaps. The flag image is compact, contained
in $\Delta_d$, and a proper closed subset of $\Delta_d$ misses
positive volume, since $\Delta_d$ is the closure of its interior.
Therefore equality holds iff every flag image is $\Delta_d$.
For a strictly positive gap vector the unique source $t$ has
$\sum g_C(t)=t(F_0)-\min_Pt$. Normalize by this positive sum.
Its normalized gap vector is on the top facet of $\Delta_d$;
it belongs to the image exactly if $W_P(t)=\sum g_C(t)$.
Conversely, if this identity holds for all positive gap systems,
the inverse of every interior simplex point belongs to $Q$;
closedness fills the boundary. Hence sharpness is equivalent to:
$$
 \text{if all flag gaps are strictly positive, then its initial
 vertex }F_0\text{ maximizes }t\text{ on }P. \tag{B}
$$

**Heredity.** If a flag in a facet $F$ has positive gaps for an
affine $f$ on $F$, extend $f$ to $\bar f$ on $P$ and subtract
a sufficiently large positive multiple of a height $h$ vanishing
on $F$ and strictly positive at the other vertices. Then the last
gap is positive, and (B) says the initial vertex maximizes $f$
on $F$. Therefore sharpness of $P$ implies sharpness of each
facet, and iterating proves it for every face.

**A facet extension.** Suppose $P$ is sharp, let $F$ be a facet
and $f$ be generic on its vertices (distinct values). Pulling $F$
by increasing $f$ gives at least one nondegenerate top-dimensional
simplex: the vertex minima of its indexing flag are distinct, so the
flag gaps are positive. This follows by recursively coning the
minimum vertex to the facets not containing it; no assumption of a
simplicial facet. By heredity and (B) on $F$, the initial vertex
maximizes $f$ there. Put $a=\min_F f,b=\max_F f$, choose an affine
extension $\bar f$, and use the same positive height $h$ as above.
At the threshold
$$
 \lambda_*=\max_{v\in\operatorname{vert}P\setminus F}
                  \frac{a-\bar f(v)}{h(v)}
$$
all outside values of $\bar f+\lambda_* h$ are at least $a$.
For every $\lambda<\lambda_*$, some outside value falls below $a$,
giving a full positive-gap flag in $P$ on appending $P$ to the
fixed flag in $F$. From (B), all values on $P$ are at most $b$.
Letting $\lambda\uparrow\lambda_*$ proves an extension with range
exactly $[a,b]$. Finally approximate any tied $f$ by generic ones,
choosing $\bar f$ continuously in $f$ by fixing its value at one
point off $\operatorname{aff}F$. The displayed finite maximum is
continuous, so the conclusion holds without genericity. In
dimension one the initial 0-face flag has vacuously positive gaps
and the same threshold argument applies.

**Width versus sections.** Set $L=L_F$ and $D=P-P$, and let
$\pi$ be restriction of covectors to $L$. The identity
$\pi(D^\circ)=(D\cap L)^\circ$ follows by taking polars of both sides:
for $x\in L$, $\sup_{t\in D^\circ}t(x)\le1$ iff $x\in D$.
These convex sets are compact with zero in their relative interiors.
Thus the same-width extension property is equivalent to
$\pi(D^\circ)=(F-F)^\circ$, equivalently $D\cap L=F-F$.
Matching widths implies matching both endpoints after adding a
constant, which does not change the restriction as an affine
functional on $F$ once one value is anchored. Applying the facet
result successively through a face chain gives (A) in the forward
direction, since the direction spaces are nested.

**Converse.** Suppose (A) holds. Within any facet $F$, a face
$E\subseteq F$ has
$(F-F)\cap L_E=(P-P)\cap L_F\cap L_E=E-E$, so by induction $F$
is sharp. Let a flag on $P$ have positive gaps under $t$, and
let $F$ be its penultimate face. Its initial vertex maximizes
$t|_F$ by (B) on $F$. By (A), there is an extension $u$ of
$t|_F$ whose range is exactly that on $F$. The two affine
extensions differ by $t-u=\lambda h$. Some outside value is below
the common minimum on $F$, and no outside value under $u$ is;
thus $\lambda<0$. It follows that $t\le u$ on $P$, and the
initial vertex maximizes $t$. Use (B). A point is taken sharp with
its usual intrinsic volume 1. This completes the proof of (A).

## 2. Symmetric equality forces two levels in every facet direction

Let $P=-P$ be sharp and choose coordinates $(x,z)$ so a given
facet $F$ is at $z=1$ and $-F$ is at $z=-1$. Let $F_+$
be the horizontal projection of $F$. For any horizontal linear
functional $f$, put $a=\min_{F_+}f,b=\max_{F_+}f$. By the extension
criterion an affine functional agreeing with $f(x)$ at $z=1$
has range $[a,b]$ on $P$. It has form $f(x)+c(z-1)$. Its
range on $-F$ is $[-b-2c,-a-2c]$, an interval of the same
width, so inclusion in $[a,b]$ forces $c=-(a+b)/2$. Consequently
every $(x,z)\in P$ satisfies
$$
 f(x)\le \tfrac{1+z}{2}b-\tfrac{1-z}{2}a.
$$
This holds for every $f$; the right side is the support function
of $\tfrac{1+z}{2}F_+ + \tfrac{1-z}{2}(-F_+)$, exactly the horizontal
projection of the slice $z$ of $\operatorname{conv}(F\cup -F)$.
Hence $P=\operatorname{conv}(F\cup -F)$. As $F$ was arbitrary,
$P$ is 2-level. This does not invoke Mahler's inequality, the
conjectured equality classification, or any numerical calculation.

For central 2-level $P$, the facet conditions of (A) at the top
rank are automatic: $P\cap L_F$ is the central slice of
$\operatorname{conv}(F\cup -F)$, namely $\tfrac12(F-F)$. Thus
a central 2-level polytope is sharp iff all its facets are sharp.
For a general 2-level $E$ and facet $F$, let $G$ denote its
opposite face (possibly lower dimensional). By the extension
criterion, the corresponding facet condition in (A) is simply
$G-G\subseteq F-F$. In fact an extension shifts the values on
the other level by a common constant; it can fit inside the
range on $F$ exactly iff the opposite width is no greater.
This is an exact finite test on every face–facet incidence,
using vertex differences and convex-hull membership.

## 3. Hansen polytopes: Conjecture 3(i) holds on this entire class

Use Freij, Henze, Schmitt and Ziegler (2013),
[arXiv:1201.5790](https://arxiv.org/abs/1201.5790), §2, §3,
and note their Hanner criterion is for *Hansen* polytopes (threshold
graphs), not for the different unconditional construction (cographs).
Write $H(G)=\operatorname{conv}\{\pm(1,\mathbf1_S):S
\text{ stable in }G\}$ for a perfect graph $G$.
Its facet at first coordinate 1 is $\operatorname{stab}(G)$,
and $\operatorname{stab}(K)$ is a face thereof for any induced
subgraph $K$, by setting the coordinates outside $K$ to zero.
By the cited Hansen duality $H(G)^\circ$ is linearly equivalent
to $H(\overline G)$. Sharpness is invariant under polarity on
origin-symmetric polytopes by Winter's (2) and flag duality.

For the path $1-2-3-4$, in its stable set polytope take the facet
$F:x_2+x_3=1$ and opposite face $G_0:x_2+x_3=0$.
Along $x_1+x_2+x_3+x_4$, widths on $F,G_0$ are respectively
1 and 2. For the cycle $1-2-3-4-1$, take $F:x_1+x_2=1$
and $G_0:x_1+x_2=0$, with functional $-x_1+x_3-x_4$.
Values on $F$ lie in $[-1,0]$, both attained; on $G_0$ they
span $[-1,1]$. The listed hyperplanes really are facets with
four affinely independent stable 0/1 vectors on each (check in
the companion replay). Both obstructions fail (A) strictly.

Thus a sharp Hansen polytope's graph has no induced $P_4$
or $C_4$, and no induced $2K_2$ because its complement
would contain $C_4$. These are exactly the forbidden induced
graphs for threshold graphs, hence $H(G)$ is Hanner by the
cited classification. Conversely Hanner polytopes have known
volume product $4^d/d!$ and $2^dd!$ flags, giving equality.
(The unweighted threshold graph characterization is stated in
the Freij et al. paper, Definition 3.2 and the following text;
the Hansen/Hanner classification is their Theorem 3.3.)

## Current scope / priority

* The 4D symmetric lower Mahler bound was preempted by the
  [September 22, 2026 all-dimensional preprint and accompanying formalization](https://github.com/openai/math/tree/main/preprints/The-symmetric-Mahler-conjecture-and-its-equality-cases-September-22-2026).
  None of the arguments above uses that announced result.
* Winter's [November 2025 primary slide deck, last slide](https://martinwintermath.github.io/pdf/slides/2025_kolkom.pdf)
  reports that they checked centrally symmetric 2-level types through
  dimension 8 and found only Hanner equality cases for the flag–volume
  bound. Coupled with §2 above, this gives the all-central-type
  conclusion in that range **conditional on their reported
  computations**. It is not an independent reproduction for 7 or 8.
* Exact population-level obstruction (A), its 2-level consequence,
  and the Hansen specialization are new candidates relative to the
  sources checked, not a claim of priority or of solving either of
  Winter's general Conjecture 3 clauses. The full centrally symmetric
  classification now reduces to centrally symmetric 2-level
  polytopes. The non-symmetric clause still lacks even a general
  2-level reduction.
