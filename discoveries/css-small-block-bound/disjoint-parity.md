# CSS large-block bound with disjoint residual parity constraints

In the intrinsic pair-block contraction of the [original CSS note](https://github.com/ConfusedDuckCollectibles/Research1/tree/41abe5308b40ffa4fc7c8d426eb1ee1768528ed7/discoveries/css-small-block-bound), \(A,B\) are the surviving X and Z matrices, and every block in \(L\) contains at least \(D\) physical qubits. A compressed X representative has minimum physical lift weight equal to its number of blocks, while a compressed Z representative occupies each block in full. Let \(d=d_X\), \(K=\operatorname{row}(B|_L)\), and \(h=k-\ell\), the dimension of the X logical image of \(K^\perp\) with zero extension off \(L\). In particular, \(K^\perp\) itself may be larger than \(h\). Stabilizer membership below always means membership in the **full** \(\operatorname{row} A\) after zero extension, not the projection of \(A\) to \(L\).

## Recognizing the residual parity constraints

Let \(U\) be the span of the weight-one and weight-two vectors of \(K\). On \(L\), connect the two coordinates of any pair in \(U\) and pin every component containing a unit in \(U\). Then \(U\) contains every word on each pinned component and exactly the even words on each unpinned component, independently on different components. This follows by adding pairs along paths and, for the pinned components, adding the one unit. Conversely its generating units and pairs have this form, so no other words occur.

Let \(C_1,\ldots,C_M\) be the unpinned components and \(c_i=|C_i|\) their numbers of blocks. The map
\[
\pi:K\longrightarrow\mathbb{F}_2^M,\qquad
\pi(z)_i=z\cdot 1_{C_i}
\]
has kernel \(U\): zero parity on each unpinned component leaves precisely the subspace just described. Set \(R=\pi(K)\). For any residual equation one can subtract an element of \(U\) from a representative in \(K\) and leave one coordinate in each affected unpinned component and none anywhere else. Thus no nonzero word of \(R\) has weight one or two; such a representative would itself be a unit or pair in \(K\), hence in \(U\). A vector in \(K^\perp\) is zero on pinned components and has one common bit per unpinned component. Those bits must lie in \(R^\perp\).

If each physical Z generator has weight \(W_Z<4D\), its projection to \(L\) involves at most three blocks, and its residual equation at most three components. If in addition \(\dim R=1\), a presented generator must give the unique nonzero residual equation, so it has weight exactly three. This argument does not imply disjointness at higher rank. In general the disjoint-triple criterion is intrinsic: enumerate weight-three words of \(R\) and test that they are a pairwise disjoint basis. In fact, if a disjoint triple basis exists then the complete list of weight-three words is precisely that basis, because adding two or more disjoint triples has weight at least six. Zero residual rank recovers the original pair-generated case.

## Disjoint-triple bound and exact restricted minimum

Suppose \(R\) has a basis of \(r\) triples with pairwise disjoint supports. Then
\[
 n_L \ \ge\ D\bigl[hd-\min(r,\lfloor h/2\rfloor)\lfloor d/2\rfloor\bigr].
\]
Add \(n_0+n_S\) to obtain the physical bound on \(n\). No claim of disjoint logical images is needed.

In component coordinates, \(R^\perp\) is the direct sum of one-dimensional unaffected components and the even subspace on each triple. In any order, map the regions into the full X logical quotient and write \(j\) for the increment in the dimension of the accumulated image due to the current region. These increments add to \(h\). An increment-one unaffected component must be a nontrivial logical and uses at least \(d\) blocks. A triple with increment one contains some nontrivial pair, so uses at least \(d\) blocks (in fact \(d+1\) counting its third component). With increment two, the triple's two-dimensional kernel injects into the quotient modulo the earlier images. Thus **all three** pairs are nontrivial even in the full logical quotient. Their weights \(c_1+c_2,c_1+c_3,c_2+c_3\) are at least \(d\), so the region uses at least \(\lceil 3d/2\rceil=2d-\lfloor d/2\rfloor\) blocks. At most \(r\) and at most \(\lfloor h/2\rfloor\) triples have increment two. Summing on disjoint block supports, each of cost at least \(D\), proves the claim.

The exact minimum X weight supported in \(L\) is the least size \(c_i\) of an unaffected indicator or size \(c_i+c_j\) of a same-triple pair which is not in the full \(\operatorname{row}A\); it is infinite if there are no such candidates. Indeed a shortened word decomposes into these pieces on disjoint supports, so a nontrivial word contains a nontrivial piece of no larger weight. This need not be the global \(d_X\) when \(S\) is nonempty.

## Equality, degeneracies and necessity of disjointness

For \(d\ge2\), integer \(D\ge1\) and \(h\ge 2r\) with \(h>0\), equality is attained by a repetition construction. Make all original blocks have width \(D\), with adjacent X pairs inside each. For each triple take component lengths \((\lfloor d/2\rfloor,\lceil d/2\rceil,\lceil d/2\rceil)\) in blocks; put consecutive Z pairs of blocks inside each component and a Z triple covering one block from each. Include \(h-2r\) independent chains of \(d\) blocks, with consecutive Z pairs. The resulting compressed \(A\) is zero; the shortened X words are whole chains and pairs of components in triples, so \(k=h\) and \(d_X=d\). A Z block unit is not stabilized: within a triple each row word has the same block-count parity on all three component chains, whereas a unit changes just one; within an unaffected chain every row word is even. Therefore \(d_Z=D\). The physical Z check maxima are \(3D\) when \(r>0\) and \(2D\) otherwise. For example an odd-distance triple has lengths \(1,2,2\), giving \([[5D,2,(3,D)]]\).

Neither the displayed equality nor its construction is claimed for \(d=1\) or \(h<2r\). The main bound still holds; keeping the positive size of unused components gives the following stronger optional bound. For \(d\ge2\), the costs per triple of increments \(0,1,2\) are at least \(3,d+1,2d-\lfloor d/2\rfloor\), respectively. Above a baseline \(hd\), these corrections are \(3,1,-\lfloor d/2\rfloor\); they decrease with larger increments. If \(h<2r\) a minimizing allocation uses all \(h\) increments on triples. Replacing two 1-increments by a 0 and a 2 never increases the bound, since \(3-\lfloor d/2\rfloor\le2\). For \(h\ge2r\) put increment two on every triple. In formulas, put \(p=\min(r,\lfloor h/2\rfloor)\), \(s=\min(r,\lceil h/2\rceil)-p\); then
\[
 n_L/D\ \ge\ hd-p\lfloor d/2\rfloor+s+3(r-p-s).
\]
For \(d=1\), each triple still occupies at least three blocks and contributes at most two dimensions, giving \(n_L/D\ge\max(3r,h+r)\). The case \(h=0\) has no shortened X logical and the restricted minimum is infinity.

Disjointness is necessary for the asserted inference, not just for this proof. On six size-\(D\) blocks, use X pairs within each and compressed Z checks on \(\{1,2,3\},\{3,4,5\},\{1,5,6\}\). The seven nonzero words in \(K\) have weights three or four, so \(U=0\), and \(R\) has rank three with \(W_Z=3D\). The six columns of the displayed check matrix are nonzero and pairwise distinct, with a dependent triple, so the X distance is three. A Z unit is not stabilized, so the parameters are \([[6D,3,(3,D)]]\). Substituting rank three in place of disjoint triple count would incorrectly require \(n_L\ge8D\).

## Optional general disjoint parity regions

The exact restricted-minimum test applies to pair pieces inside **any** disjoint residual parity supports, since an even set splits into disjoint pairs. A simple count also applies to regions of size \(t\ge3\). Let \(j\) be a region's increment as above. Choose any component \(v\) having the smallest length \(a\ge1\), label its own component by zero, and label component \(i\ne v\) by the image of the pair \(\{i,v\}\) after quotienting by preceding regions' images. These labels span dimension \(j\); at least \(j\) are nonzero. Each nonzero-labeled component has length at least \(\max(a,d-a)\), while the remaining components each have length at least \(a\). The total is at least
\[
 (t-j)a+j\max(a,d-a).
\]
For \(d\ge2,\ j\ge1\), minimize over positive integer \(a\) to get \(jd+\min(t-2j,(t-2j)\lfloor d/2\rfloor)\). On \(1\le a\le\lfloor d/2\rfloor\) the objective is linear with slope \(t-2j\). On \(a\ge\lceil d/2\rceil\) it is \(ta\); for odd \(d\), its first value there is larger by \(t-j\) than the value at \(a=\lfloor d/2\rfloor\). For \(j=0\) or \(d=1\) the cost is at least \(t\). Sum these costs and \(d\) per contributing unaffected component, then multiply by \(D\). A rank-independent bound can be read off by minimizing that sum over integer allocations of the \(h\) increments within their dimension limits. No global equality assertion is made in the partially contributing case. When all parity regions contribute their full \(t-1\) dimensions and enough independent chains are included, the same construction with component lengths \(\lfloor d/2\rfloor,\lceil d/2\rceil,\ldots,\lceil d/2\rceil\) attains the full-increment local terms for \(d\ge2\).

The local two-dimensional cost is the dimension-two case of the [classical Griesmer bound](https://www.bitsavers.org/pdf/ibm/IBM_Journal_of_Research_and_Development/045/ibmrd0405M.pdf). The complementary logical dimension is the CSS [cleaning identity](https://arxiv.org/pdf/1011.3529), and the equality controls use standard [CSS repetition copying](https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.5.040302). The point here is the intrinsic residual criterion and its finite charge on the already-contracted large region, not a new general Griesmer inequality or a code performance record.
