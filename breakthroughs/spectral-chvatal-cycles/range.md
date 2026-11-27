# Constructive range theorem for the arbitrary-cycle Chvátal equations

## Abstract range statement

Let the old coordinates be partitioned into tight coordinates I and tips J, with |J|≥5. Let U be a family of old triples containing every triple with one coordinate in I and two in J. Other triples can be absent. Write X for the U-by-I incidence matrix. Define the linear map A(Δ)=ΔX for symmetric matrices Δ supported on disjoint pairs of triples. Then its image over the rationals or reals is exactly the U-by-I matrices B satisfying these two conditions.

- B(U,i)=0 whenever i belongs to U.
- The square matrix XᵀB on I is symmetric.

Both are necessary. A disjoint triple cannot contain any coordinate of U, and XᵀΔX is symmetric.

For sufficiency, use the reference triples {i}∪P with i in I and P a two-subset of J. In each non-reference row U, for every i outside U in I use one edge to any reference triple {i}∪P disjoint from U. Such a P exists since at most three tips are in U. Put the desired B(U,i) on that edge. The edge touches no other tight column in its non-reference row. Thus the resulting error is supported on reference rows. Subtracting this partial symmetric correction preserves the two compatibility conditions.

Now fix two distinct coordinates i,j in I. The remaining entries in column j of the i-reference rows and in column i of the j-reference rows have equal totals by symmetry of XᵀB. Edges between these two groups form the bipartite disjointness graph of two copies of the tip pairs. It is connected. Any two left pairs sharing a tip have a common right neighbor because their union has three tips and at least two other tips are available. The graph of left pairs linked by overlap is connected, and every right pair has a left neighbor. On a connected bipartite graph, arbitrary signed vertex-sum demands are realized by edge weights exactly when the left and right totals match. Choose a spanning tree and eliminate leaves, assigning the leaf's residual demand to its unique remaining edge. Equal totals leave zero at the root. Different unordered tight coordinate pairs affect disjoint sets of reference equations, so these tree solutions combine. This proves sufficiency and supplies an explicit right inverse with no pseudoinverse or denominators beyond those in B.

## Apply it to the arbitrary 2-regular anchor problem

Start with n≥2r+1 old points, r anchors, distinct tips for the r edges of a simple 2-regular graph on the anchors, and all remaining points untouched. For each edge delete the old triple consisting of its anchor endpoints and its tip. All surviving old triples have weight 1. Add one bridge for each anchor, all sharing the common marker and with weight 2. The tight set I consists of all anchors and untouched points. Every tight star has weight K=binom(n−1,2). Each tip has old star weight K−1.

Use the normalized Gram convention with diagonal 1/w, zero on distinct intersecting faces, and the equation mapping each weighted tight-star indicator to the all-ones vector. If this normalized matrix is PSD then K times it minus the all-ones matrix is PSD by Cauchy–Schwarz and certifies the residual theta bound. On the bridge block the normalized Gram is I/2 because all bridges intersect. Let C be its triple-by-bridge block and Y the bridge-by-tight-coordinate weighted incidence, so Y has entry 2 for its own anchor and zero otherwise. The bridge row condition is
    CᵀX+(1/2)Y=J_bridge,tight.
There is a simple explicit C. For every bridge i and old triple U not containing i, start with 1/q where q=binom(n−2,2). Use zero for triples containing i. For every other anchor j, add (2−e_ij)/(q binom(r,2)) to each of the reference triples {j} plus two tips, where e_ij is one when ij is an edge and zero otherwise. The uniform part misses 2−e_ij removed triples in the j-star and none in any untouched star. The correction is supported on reference triples with only j tight, so it restores exactly those deficits. This proves the entire bridge-row equation directly.

For any symmetric initial triple block H0 with diagonal 1 and zeros on distinct intersecting triples, the remaining right-hand side B is J_triple,tight−CY−H0X. If i is in U, the H0 term is exactly 1 and the bridge term is zero, making B(U,i)=0. Also
    XᵀCY = (J_bridge,tight−Y/2)ᵀY.
Writing y=Yᵀ1 and using Xᵀ1=K1−y, the difference XᵀJ_triple,tight−XᵀCY is KJ−y1ᵀ−1yᵀ+YᵀY/2. This is symmetric, and XᵀH0X is symmetric. Thus the range theorem applies to this C or to any correctly normalized C satisfying the bridge equations and intersection zeros. It produces a rational Δ for every simple 2-regular anchor graph with r≥5 and every n≥2r+1, including disjoint unions of cycles. It also applies to more general deletion patterns whenever the same incidence and compatibility hypotheses hold.

PSD is a separate requirement. The particular sparse tree correction at C5, n=11 has a negative triple-block principal 3-by-3 minor on masks 97, 386 and 400. Its scaled determinant is −19,125,227,520 over 7560³, exactly. This construction proves exact affine feasibility and leaves the coupled Schur-complement inequality open. The sparse tree solution is an explicit particular solution. Any search for a PSD solution can add symmetric disjoint-supported matrices annihilating X. Averaging a particular solution under permutations of untouched points preserves all constraints when starting from the explicit C and the uniform Kneser H0.

## Exact replay and scope

Run `node range-replay.cjs`. The replay takes H0=I+A_disj/binom(n−4,2) on surviving triples and checks all bridge and triple star equations with BigInt. It constructs valid corrections on C5 at n=11, C6 and C3+C3 at n=13, and C3+C4 at n=15. These are affine checks. The r=3,4 cases lie outside the five-tip range proof, and the square certificate reported by Gemini is not used here and was not replayed end-to-end for this release. [Current Carbon share](https://share.gemini.google/DVKV2JFxKpUc) · [Ellis, Filmus and Friedgut, §5](https://arxiv.org/html/2609.28404v3#S5)