# Fuglede's conjecture holds in F₁₃³

Let $G=\mathbb F_{13}^{3}$ and $\zeta=\exp(2\pi i/13)$. A nonempty subset E is **spectral** if there is a set A of size |E| such that the restrictions to E of $x\mapsto\zeta^{a\cdot x}$ for a∈A are orthogonal. It **tiles by translation** if some T⊂G gives a unique representation of each element of G as e+t with e∈E, t∈T. The theorem is the equivalence of these two conditions in G.

## 1. The standard reductions, including both implications

For nonzero v, a Fourier zero of an indicator means that the 13 integer counts on the levels v·x=j are all equal: their generating polynomial of degree at most 12 vanishes at ζ, and the minimal polynomial is $1+z+\cdots+z^{12}$. This also shows that Fourier zeros are invariant under multiplication of v by any nonzero field element.

For a spectral pair with 1<|E|=|A|<2197, any nonzero difference in A gives equal counts on the levels of E, so 13 divides |E|. If |E|>169, each direction is determined by some pair in E (apply pigeonhole to the 169 lines in that direction). Then orthogonality of the rows, equivalently the spectrality of A with spectrum E, would force every nonzero Fourier coefficient of 1_A to vanish. Fourier inversion would make 1_A constant, impossible for a proper nonempty set. Thus apart from sizes 1 and 2197, the only candidates are multiples of 13 through 169.

The spectral sizes 1 and 2197 tile immediately. At size 13, a nonzero spectral difference gives exactly one point on each level, so E tiles by the kernel plane. At size 169, the same inversion argument guarantees a direction missing from the differences of E: otherwise A would have no nonzero Fourier coefficient. Thus E is a full system of representatives of cosets of that line, so tiles by it.

Conversely, if E tiles with T, then |E|·|T|=2197, so |E| is one of 1, 13, 169, 2197. Size 1 and the whole group are spectral. For size 13, T is proper, and some nonzero v has a nonzero Fourier coefficient of 1_T. The Fourier transform of the tiling identity forces 1_E to vanish there. Thus E has one point on each v-level, and the characters with frequencies on the line spanned by v are a spectrum. For size 169, T has 13 points. By the cone-invariance of its Fourier zeros, the same Fourier tiling identity holds with any nonzero scalar multiple cT replacing T. This implies the projective nonzero difference-direction sets of E and T are disjoint: if e−e′=c(t−t′) with distinct pairs, then e+ct′=e′+ct, contradicting uniqueness. Take any nonzero direction of T. E meets every coset of that direction line once, so the character frequencies in the annihilator plane of that line form a spectrum. These standard facts also appear in [Aten et al.](https://arxiv.org/abs/1509.01090) and [Ferguson–Iosevich, §8](https://arxiv.org/html/2610.30712v1).

## 2. All possible non-tiling spectral sizes

| Size | Exact exclusion |
| --- | --- |
| 26 | [Projection and planar-cap proof](../../discoveries/fuglede-f13-26/PROOF.md) |
| 39 | [Height and determinant proof](../../discoveries/fuglede-f13-39/proof.md) |
| 52 | [Layered and nonlayered moment proof](../../discoveries/fuglede-f13-52/) |
| 65, 91 | [Classical cyclotomic determinant obstruction](../../discoveries/fuglede-f13-39/STATUS.md#the-existing-obstruction-at-65-and-91) |
| 78 | [Entire exclusion in the four steps below](#3-the-78-point-case) |
| 104, 117 | [Planar caps and moment proof](../../discoveries/fuglede-f13-104-117/PROOF.md) |
| 130, 143, 156 | [Malikiosis, Theorem 1.4](https://link.springer.com/article/10.1007/s43670-023-00078-2): all lie strictly between $169-12\sqrt{13}$ and 169. |

For 65 and 91, the character matrix H has determinant Δ with $\Delta\overline\Delta=N^N$. Use q=5 for 65 and q=7 for 91. A power of q is −1 modulo 13, so complex conjugation fixes every prime over q in the unramified cyclotomic field and the left-hand q-prime valuation is even. The right has the odd valuation N, as q occurs once in N. This is an application of [Brock's self-conjugacy obstruction](https://doi.org/10.1016/0097-3165(88)90054-4), not a new result.

These eleven candidate cardinalities are all the non-tiling multiples of 13 strictly between 13 and 169. The exclusions therefore complete the other implication, spectral implies tile.

## 3. The 78-point case

The four linked parts give the full reduction and checked finite input, with independent implementations:

1. [Neither member layered](01-neither-layer/PROOF.md). Plane caps and the Radon profile identity leave only ten- and eleven-uniform pencil possibilities. The six-fibre label gate and incidence deficit eliminate them. **At least one** member is in six parallel planes, with 13 selected points per plane.
2. [Exactly one member layered](02-one-layer/PROOF.md). Three common missing directions in the layered member give a complete-map atlas. Every choice of six levels contains a three-term arithmetic progression. Only two anchored triples can have too few cross-layer directions; neither extends to a fourth layer. A weighted contradiction then forces **both** members to be layered.
3. [Both members layered](03-linked-layers/PROOF.md). The common-slice direction argument rules out a six-point line and ensures that **some** pair of layer normals is perpendicular. This permits simultaneous linked graph representations, with no assumption that arbitrary layer normals are perpendicular. Orthogonality puts every active cross-cut in a five-slope family, now for the entire union of cross-layer directions, not just each pair.
4. [Final labelled graph lemma](04-final/PROOF.md). Across all possible five-slope palettes, a single fixed pair of the six labelled positions has height difference uniquely determined by its palette. Two active cuts force a nonzero-step finite derivative to be constant, making a 13-point layer graph affine. This contradicts the five-point line cap.

## Relation to prior work

[Fallon–Mayeli–Villano](https://arxiv.org/abs/1902.02936) prove the smaller odd-prime cases 5 and 7 and develop the direction geometry. [Ferguson–Iosevich](https://arxiv.org/html/2610.30712v1), especially §§8 and 10–13, prove the 11 case and provide the projection, factor-line, moment and weighted incidence methods extended here. [Qihang Wang](https://arxiv.org/html/2609.20495v1) gives the planar 52-point obstruction used in the earlier 52 reduction. The active-column step invokes the classical prime Fourier minor theorem ([Frenkel](https://arxiv.org/abs/math/0312398); see also [Tao](https://arxiv.org/abs/math/0308286)). This proof settles the stated finite grid; it asserts neither the general prime-field conjecture nor priority over unknown or unpublished work.
