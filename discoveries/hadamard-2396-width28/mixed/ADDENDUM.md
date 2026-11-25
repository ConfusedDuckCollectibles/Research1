# Mixed scalar obstruction and fixed-sum quotient tables at Hadamard order 2396

This is a companion to [the pure-scalar width-28 restriction (#22)](https://github.com/ConfusedDuckCollectibles/Research1/blob/457c7d3facc564fc8f2f7b12cacc4d75d185d921/discoveries/hadamard-2396-width28/NOTE.md). It gives two exact statements for the published abelian index-eight Goethals–Seidel template:

1. A specified mixed affine symmetry cannot have scalar order divisible by \(9\) on the order-37 coordinate, even when it moves the other coordinate and permutes or negates the four blocks.
2. For elementary quotient \(Q=C_2^3\), restoring the actual trivial-character norm 599 to the known normalized quotient classification gives four row-sum profiles, 44,040,192 labeled signed quotient arrays, and 32 patterns under independent row translations and signs and row permutations. None is asserted to lift to a full core.

These do not settle order 2396, exclude a group type without a symmetry assumption, or classify all automorphisms of any resulting Hadamard matrix.

## Native template and compressed constraints

The starting condition is [Egan–Kotsireas–London–Ó Catháin, Theorem 7](https://arxiv.org/html/2610.27175v2) and [Hexagon construction note, §3.2 Eq. (4)](https://hexagonmath.org/pdf/2611.00637v1). For four sign functions \(a_f:G\to\{\pm1\}\) on an abelian group \(G\) of size 592, the target with index-eight subgroup \(N\) is
\[
 R(s):=\sum_f\sum_x a_f(x)a_f(x+s)=2396\mathbf1_{\{0\}}(s)-32\mathbf1_N(s)+4 .
\]
In particular it is **+4 off \(N\)**, \(-28\) on \(N\setminus\{0\}\), and 2368 at zero. Write \(G=A\times\mathbb F_{37}\), \(N=\langle e\rangle\times\mathbb F_{37}\), where \(|A|=16\) and \(e\) has order two. Let \(Q=A/\langle e\rangle\), \(z_f(u)=\sum_y a_f(u,y)\), and
\[
 \gamma_f([u])=\frac{z_f(u)+z_f(u+e)}2\in\mathbb Z,\qquad
 s_f=\sum_{c\in Q}\gamma_f(c),\qquad t_f=\sum_{c\in Q}\gamma_f(c)^2 .
\]
Fourier transform the native target on characters trivial on \(N\). The power sum on \(G\) is 2396 at the trivial character (because \(-32\cdot74+4\cdot592=0\)) and 28 at each of the seven others. The compressed half-sums therefore satisfy, with unnormalized transform,
\[
 \sum_fs_f^2=599,\qquad
 \sum_f|\widehat\gamma_f(\theta)|^2=7\ \ (\theta\ne1),\qquad
 \sum_ft_f=81,\qquad
 \sum_f(8t_f-s_f^2)=49 .
 \tag{1}
\]
The latter two follow from Parseval. Each row variance \(8t_f-s_f^2=\sum_{\{c,c'\}}(\gamma_f(c)-\gamma_f(c'))^2\) is nonnegative; the sum is over the unordered pairs of distinct cosets. Every compressed entry satisfies \(|\gamma_f(c)|\le7\). Otherwise that entry consumes at least 64 of the total square norm 81, each other entry in its row has magnitude at most 4, and their seven comparisons contribute at least \(7\cdot4^2=112>49\).

## Elementary quotient: normalization and restoration

Identify \(Q\) with \(\mathbb F_2^3\). Its Walsh coefficients are integers and each \(\widehat\gamma_f(\theta)\equiv s_f\bmod2\), because the character values are \(\pm1\). The only absolute-value pattern of four squares summing to 7 is \(2,1,1,1\). Thus one fixed even row has absolute value 2 at all seven nontrivial characters, and the three fixed odd rows have absolute value 1. The row variances \(8t-s^2\) are 28 and 7 respectively. Modulo 8 an odd row has \(s\equiv\pm1\bmod8\) and the even row \(s\equiv\pm2\bmod8\).

Each odd row is uniquely restored from a normalized row of sum 1 by an integer constant and the sign prescribed by \(s\bmod8\). The norm of the normalized row is \((1+7)/8=1\); it is a \(+1\) delta at any one coset. The even row is restored in the same way from a normalized row of sum 2 and norm \((4+28)/8=4\). Either it is a \(+2\) delta, or it has four distinct nonzero entries, three \(+1\) and one \(-1\). In the latter case all nontrivial autocorrelations must vanish: the squared Walsh coefficients are constant 4 including at the trivial character. At a nonzero difference from the negative support point to a positive point, the negative pair contributes \(-2\); the remaining two positives must make a pair with that same difference. Equivalently the four points sum to zero in \(\mathbb F_2^3\), i.e. they form an affine plane. Conversely each such plane pattern gives zero autocorrelation off zero, since the two pairs at each nonzero difference in its underlying plane have opposite sign products. Delta rows have the required spectra as well. There are eight odd normalized rows and \(8+14\cdot4=64\) even normalized rows (each of the seven linear planes has two cosets).

Of the 14 sorted nonnegative four-square representations of 599, the only ones with three entries \(\pm1\bmod8\) and one \(\pm2\bmod8\) are

| Absolute row sums | Even row sum magnitude | Row label assignments |
|---|---:|---:|
| \((1,7,15,18)\) | 18 | 24 |
| \((2,9,15,17)\) | 2 | 24 |
| \((6,7,15,17)\) | 6 | 24 |
| \((7,10,15,15)\) | 10 | 12 |

The complete enumeration uses a two-pair join on squares \(0^2,\ldots,24^2\). Conversely, for each profile, restore any independent normalized rows to any assigned signed row sums. The total squared trivial coefficient is 599 and the divided power sum for each other character is \(4+1+1+1=7\), exactly (1). No other condition at the quotient level is omitted. The count for fixed labeled \(Q\) is
\[
 (24+24+24+12)\cdot 16\cdot 8^3\cdot64=44{,}040{,}192 .
\]
All sums are nonzero, so signs can make them positive and row permutations can sort them. Each odd exceptional location can be translated to zero independently. For the even row there is one translated delta type or a plane through zero with its minority sign at zero. There are seven choices of underlying plane, unchanged under translation; at fixed sum the sign of the defect relative to the normalized shape is already determined. Hence precisely \(4(1+7)=32\) patterns under the given operations. The double odd sum 15 in the last profile does not further identify even-row shapes. Inversion acts trivially on \(Q\), so does not change this count.

### Optional ambient-relabeling refinement

These 32 patterns do not mod out by common quotient automorphisms. An elementary quotient implies \(2A\subseteq\langle e\rangle\), so the only possibilities are \(A=C_2^4\) with any \(e\ne0\), or \(A=C_4\times C_2^2\) with \(e\) twice a generator of \(C_4\). In the first case, a linear complement to \(\langle e\rangle\) lifts every automorphism of \(Q\), so the seven planes form one orbit and there are \(4\cdot2=8\) classes with common ambient automorphisms included. In the second case let \(H\subset Q\) be the image of \(A[2]\), a two-dimensional plane. Any lift preserves \(H\), since \(A[2]\) is characteristic. Conversely every linear map preserving \(H\) lifts by applying the prescribed invertible linear map on the two elementary summands and sending a \(C_4\) generator \(g\) to \(g+u\) for suitable \(u\) in those summands. Its quotient action can send \(g+\langle e\rangle\) to itself plus any vector of \(H\), and it fixes \(e=2g\). The stabilizer is transitive on the six other planes: it moves their lines of intersection with \(H\) transitively and, for a fixed line, a suitable \(g\mapsto g+u\) swaps the two possibilities. Thus there are \(4\cdot3=12\) classes here, one delta type and two plane types per profile. These are ambient-liftable equivalences of quotient arrays, not of possible lifts above them.

## Mixed affine scalar obstruction

Precisely, suppose the sign functions satisfy one symmetry of the following form:
\[
 a_{\pi(f)}\bigl(\phi_f(u),\,\lambda y+d_{f,u}\bigr)=\epsilon_f\,a_f(u,y),
 \qquad \lambda\in\mathbb F_{37}^{\times},\quad \epsilon_f\in\{\pm1\}.
 \tag{2}
\]
Here \(\pi\) permutes four blocks, each \(\phi_f\) is a permutation of \(A\) mapping \(e\)-pairs to \(e\)-pairs, and the induced permutations on \(Q\) are affine with respect to its abelian group structure. The translations \(d_{f,u}\) can be arbitrary for each slice. The same \(\lambda\) acts on every slice; signs are constant on each full block. We show that \(9\nmid\operatorname{ord}(\lambda)\). In the elementary quotient case the proof actually needs only that the slice permutations preserve the pairs, not the extra affine assumption on their quotient actions.

Indeed \(s_{\pi(f)}=\epsilon_f s_f\) by bijectivity on each block. If \(\pi\) had a 3-cycle, three equal squares in \(\sum s_f^2=599\) would leave a single square congruent to \(2\bmod3\), impossible. Thus \(\pi\) has order dividing 4; the signed block permutation has order dividing 8. If \(\operatorname{ord}(\lambda)\in\{9,18,36\}\), raise (2) to the eighth power. The resulting signed block permutation is the identity and its scalar has order 9. This power is closed in the stated form (the linear coefficient of each composed \(\mathbb F_{37}\) map is \(\lambda^8\)). The row slice maps now each have order containing at most one factor 3: their permutation of eight pairs has this property, since no cycle can have length divisible by 9, and the kernel permuting inside pairs is a 2-group. Raise once more by a common multiple of the prime-to-3 parts of those four slice orders, chosen coprime to 3. This leaves the scalar of order 9 and makes each slice permutation of order 1 or 3. It also keeps the block action unsigned and fixed and preserves affinity on \(Q\).

Each such pair action fixes at least two quotient cosets, since eight points can be split only into fixed points and 3-cycles. The two slices of a fixed pair must themselves be fixed (there can be no transposition in a permutation whose order divides three). On either fixed slice a field-affine map with scalar of order 9 has a unique center and four 9-cycles. Thus \(z_f(u)=\pm1+18j\) at each fixed slice; its paired quotient value has \(\gamma_f(c)\in9\mathbb Z+\{-1,0,1\}\). The uniform bound \(|\gamma|\le7\) forces these fixed values to lie in \(\{-1,0,1\}\), called light.

- **\(Q=C_8\) or \(C_4\times C_2\).** Its affine group is a 2-group, so the induced action of order dividing 3 fixes every quotient coset. Indeed \(|\operatorname{Aut}(C_8)|=4\). For \(C_4\times C_2\), generator images have the form \((a,b),(2c,d)\), where the map on \(Q/2Q\) is invertible exactly when \(a\) is odd and \(d=1\), giving eight automorphisms. These choices are sufficient: the images span the quotient and include the nonzero element of \(2Q\) by doubling the first generator image. The affine groups have eight times the stated sizes. All quotient entries would be light, so \(\sum t_f\le32<81\), impossible.
- **\(Q=C_2^3\).** Each odd row is constant except at a single coset, which must be fixed by its invariance. That exception and at least one other fixed coset are light, forcing the entire row light. The row is of integer baseline plus or minus a single delta, so \(|s|\in\{1,7\}\). Thus the three odd row sums cannot include magnitude 15, present in all four profiles. Equivalently, without the four-square enumeration, the even sum has residue \(\pm2\bmod8\), magnitude at most 24 and square at least \(599-3\cdot49=452\). It must have magnitude 22, leaving 115 for three terms in \(\{1,49\}\), again impossible.

This kills only (2) with scalar order divisible by 9. It does not treat arbitrary position-dependent signs, unrelated slice scalars, or a general signed matrix automorphism. Nor does it extend by itself the earlier pure order-4 and order-6 obstructions to mixed actions.

## Prior work and exact controls

[Hexagon, Remark 7, Proposition 8, and Appendix C Table C1](https://hexagonmath.org/pdf/2611.00637v1) already classifies admissible order-eight quotient tables **modulo separate block constants**: for \(Q=C_2^3\) it reports 2,097,152 labeled tables, eight classes under core moves, and two after including common automorphisms of \(Q\). The normalized shapes and counts above reproduce those, since \(4\cdot16\cdot8^3\cdot64=2{,}097{,}152\) before the physical row sums are selected. The contribution here is restoring constants at the actual trivial spectrum 599, distinguishing which quotient relabelings lift to the ambient pair, and the mixed order-9 scalar obstruction. The frame transport in Proposition 8 is not itself an assertion of ambient liftability. The [Davis–Astra 2396 report](https://raw.githubusercontent.com/1337isnot1337/solved-conjectures/main/results/hadamard-2396/README.md) treats a different global signed automorphism restriction and width-4/width-12 or Legendre-pair routes. [Ramos's length-333 Legendre-pair report](https://github.com/Arthur742Ramos/hadamard-668-multiplier-obstructions/blob/main/lp333/report.md) includes related mixed scalar obstructions with an order-37 factor for order 668. No direct same-template counterpart to (2) found in this bounded comparison; no priority claim.

    python3 -B verify.py | diff -u replay.txt -
    sha256sum -c SHA256SUMS

The original exact replay independently inverts every one of the 128 Walsh sign choices per normalization, verifies the row shapes against all restored single-row spectra, counts translation classes directly, checks all three affine groups and the 225 elementary affine actions of order dividing three, and explicitly lifts the stabilizer of \(H\). It also verifies the 1-plus-4-by-9 affine orbit shape after removal of the 2-part on all 888 scalar/translation controls, and the earlier #22 Singer compressed witness with native off-pair correlation \(37\cdot4=148\) and physical profile \((1,7,15,18)\). That witness is not a full-core lift.
