# Three contacts with the Moser lattice

## Statement

Write eta=(1+i sqrt(3))/2 and w=(5+i sqrt(11))/6, and let M=Z[1,eta,w,eta w] denote the **additive** integer span of these four elements. This is the Moser lattice, not the larger Moser ring. The following are exact statements about the infinite additive group and arbitrary real-plane centers.

1. Every circle of radius one centered outside M meets M in at most three points. Four or more M-points on the circle force its center into M.
2. Three points occur. In every such case they are the vertices of an equilateral triangle of side sqrt(3). More precisely, off-lattice centers with three M contacts form exactly the four nonzero cosets of M represented, in the ordered basis above, by (0,0,1,1)/3, (0,0,2,2)/3, (1,1,0,0)/3, (2,2,0,0)/3. The results transport to any Euclidean isometric copy of M.

No claim is made about contacts with the full Moser ring, or about whether every dense finite unit-distance graph embeds in M.

## Exact reduction

For any coefficient tuple (a,b,c,d), let X=12 Re(z) and Y=36 Im(z)/sqrt(3) for z=a+b eta+c w+d eta w. Thus
    X=(12a+6b+10c+5d)−d sqrt(33),
    Y=(18b+15d)+(2c+d) sqrt(33).
By direct expansion, the norm squared is (3 X²+Y²)/432=(P+Q sqrt(33))/6, where
    P=6(a²+ab+b²+c²+cd+d²)+10(ac+bd)+5(ad+bc),
    Q=bc−ad.
This is the rational form used adjacent to Theorem 2.5 of Engel et al.

Assume at least three M contacts. Translate one contact to the origin, a translation preserving M. The origin and two other contacts are noncollinear because a line meets a circle in at most two points. Therefore the center is their unique circumcenter; solving the two linear dot-product equations puts its X,Y in Q(sqrt(33)), hence the center itself in K=Q(i sqrt(3),i sqrt(11)). Complex conjugation commutes with the field automorphism sigma fixing i sqrt(3) and negating i sqrt(11). Applying sigma to any equality z bar(z)=1 for a center-to-contact difference proves those differences have unit norm in the conjugate embedding as well. Consequently each remaining contact difference from zero has scaled squared lengths at most 4 in both complex embeddings, by the ordinary triangle inequality in each embedding.

For U=(a,b), V=(c,d), and E(s,t)=s²+st+t², the mean p=P/6 of those conjugate squared lengths satisfies
    p=(11/12) E(U+V)+(1/12) E(U−V)
     >= (1/12) [E(U+V)+E(U−V)]
     >= (a²+b²+c²+d²)/12.
The last step uses E(s,t)>=(s²+t²)/2 and the parallelogram identity for E. Since p≤4, each integer coordinate lies in [−6,6] (sqrt(48)<7). For this box, exact membership in the candidate difference list is P≤24 and (24−P)²≥33Q², equivalent to both conjugate length bounds. This is a finite list of 204 nonzero tuples, not a radius or span assumption in the physical plane.

For two members u,v of this list, put N_i=3X_i²+Y_i² and det=X_u Y_v−X_v Y_u. Skip det=0. If the center coordinates are (X/D,Y/D), the perpendicular bisector equations give
    D=6 det,
    X=N_u Y_v−N_v Y_u,
    Y=3(X_u N_v−X_v N_u).
The circle through zero has radius one exactly when 3X²+Y²=432 D². All operations are in Q(sqrt(33)); the script uses integer pairs and exact rational normalization, never evaluating a radical numerically. It merges precisely equal centers and records the union of their neighboring difference tuples across accepted pairs. A true center with >=3 contacts is counted by at least one pair; every other contact besides zero appears in a pair as the points are not collinear. Conversely an accepted pair supplies three distinct literal unit contacts. The union, plus zero, is therefore the entire intersection.

To test M membership or cosets, for normalized center coordinates X0+X1 sqrt(33), Y0+Y1 sqrt(33) (with rational coefficients), invert the displayed basis equations:
    d=−X1; c=(Y1+X1)/2; b=(Y0−15d)/18;
    a=(X0−6b−10c−5d)/12.
The enumeration keeps 2,460 pairs at 30 centers (18 lattice centers with 18 total contacts; 12 exterior centers with exactly 3). Each exterior case has all three pairwise squared distances 3 and falls in one of the four stated cosets; each coset occurs. Translation by any M element preserves the three-contact property, so this is the full classification, not merely a list of samples. A concrete sharp center with basis coordinates (−2,1,4/3,−2/3) has the three contacts with tuples (0,0,0,0), (−3,0,2,0), (−3,3,2,−2). It is not in M since its rational basis coordinates are not integral.

**Replay.** Run `node verify-moser.cjs`. The script is original, dependency-free and uses only the small finite calculation just described. The expected JSON is
    {"differences":204,"pairs":2460,"centers":30,"hist":[{"18":18},{"3":12}],"equilateral":12,"cosets":["0011","0022","1100","2200"]}

## Closest sources and result boundaries

- Engel, Hammond-Lee, Su, Varga and Zsámboki, *Diverse beam search to find densest-known planar unit distance graphs*, arXiv:2406.15317v3, https://arxiv.org/html/2406.15317v3, §2.1 Theorems 2.3 and 2.5, establish the basis and enumerate the 18 unit vectors from the origin. Their implementation and data are linked at https://codeberg.org/zsamboki/dbs-udg. The affine-center classification above is not that origin count. Their §4 finite lattice extension searches do not rule out arbitrary off-lattice points.
- The recent *A divisibility obstruction for joining two Moser lattices* companion (https://github.com/calinnedelcu/joining-rotations-artifacts ; https://doi.org/10.5281/zenodo.22855504) analyzes cosets and joins for chromatic-number applications. It is relevant arithmetic context; a matching arbitrary-translated-center contact theorem was not found in the accessible artifact summary. No priority claim is warranted without a full comparison.
- Alexeev–Mixon–Parshall, arXiv:2412.11914v2, https://arxiv.org/html/2412.11914v2, establish exact u(n) through 21 and a list of graph6 extremizers. The Appendix draws nearly all of them inside a larger arithmetic setting; it is not a theorem that dense planar configurations at n=25 reside in M.
- Modern status and the conditional exact n=25 consequences are in REDUCTION-25.md. These are finite geometric statements, independent of the disproof of the historical asymptotic conjecture.
