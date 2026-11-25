# A multiplier obstruction for width-28 cores at Hadamard order 2396

**Result.** In the published index-eight abelian width-28 Goethals–Seidel template at order 2396, a common multiplier that fixes the 2-primary coordinate and acts by a scalar on the order-37 coordinate can have order only $1,2,$ or $3$. No 2-primary group type or subgroup direction is assumed. This does not exclude width-28 cores without this symmetry, other construction routes, or a Hadamard matrix of order 2396.

## Exact template and compressions

Use the published **positive four off the subgroup** correlation, [Egan–Kotsireas–London–Ó Catháin, Theorem 7](https://arxiv.org/html/2610.27175v2) and [Hexagon construction note, §3.2, Eq. (4), Theorem 5](https://hexagonmath.org/pdf/2611.00637v1):
$$
 R(s):=\sum_{f=1}^4\sum_{x\in G}a_f(x)a_f(x+s)
 =2396\,\mathbf1_{\{0\}}(s)-32\,\mathbf1_N(s)+4,
 \qquad a_f:G\longrightarrow\{\pm1\}.
 \tag{1}
$$
Thus $R(0)=2368$, $R(s)=-28$ on $N-\{0\}$, and $R(s)=4$ off $N$. Here $|G|=592=16\cdot37$, $|N|=74$, and primary decomposition gives $G=A\times\mathbb F_{37}$, $N=\langle e\rangle\times\mathbb F_{37}$, where $|A|=16$ and $e\in A$ has order two. We assume only (1), not a particular border frame. Proposition 8 of the Hexagon note supplies frames for all three quotient types of order 8 after its indicated core-preserving moves.

The power spectral sum is $\sum_f|\widehat a_f(\chi)|^2$, where $\widehat a_f(\chi)=\sum_xa_f(x)\chi(x)$ is unnormalized. Fourier transform of (1) gives 28 at every nontrivial character trivial on $N$, and 2396 at every other character, **including the trivial character**. To see the cancellation at the trivial one, $-32\cdot74+4\cdot592=0$.

First compress over $\mathbb F_{37}$ to $z_f(u)=\sum_{y\in\mathbb F_{37}}a_f(u,y)$. The summed autocorrelation of $z$ at $d\in A$ is $\sum_y R(d,y)$, by expanding both sides over ordered pairs of coordinates with difference projecting to $d$. It has values 1360 at zero, $-1036$ at $e$, and 148 at the other fourteen shifts. Across the four rows there are 32 unordered $e$-pairs, with entries $z_i,z_i'$. Therefore
$$
 \sum_i(z_i^2+(z_i')^2)=1360,\qquad
 \sum_i(z_i+z_i')^2=1360-1036=324.
 \tag{2}
$$
The $e$-correlation counts the two directed contributions to each unordered pair, exactly the cross term in the second identity.

Then compress over $N$ to $Q=G/N=A/\langle e\rangle$. Put $\beta_f(c)=\sum_{x\in c}a_f(x)$, $\gamma_f(c)=\beta_f(c)/2\in\mathbb Z$ (cosets have even size 74), and set $s_f=\sum_c\gamma_f(c)$, $t_f=\sum_c\gamma_f(c)^2$. The seven nontrivial characters of $Q$ give divided power sum 7, and the trivial one gives 599. Parseval on this quotient yields, irrespective of its abelian group type,
$$
 \sum_f s_f^2=599,\qquad
 \sum_f t_f=(599+7\cdot7)/8=81,\qquad
 \sum_f(8t_f-s_f^2)=49.
 \tag{3}
$$
Each row term $8t_f-s_f^2=\sum_{c<c'}(\gamma_f(c)-\gamma_f(c'))^2$ is nonnegative; the unordered-pair sum uses any ordering of $Q$.

## Scalar order exclusions

A scalar of order $h\mid36$ fixes one position in each $37$-point slice and has $36/h$ other orbits of length $h$. If the sign function is fixed, its slice sum is $z=c+h\,j$, $c\in\{\pm1\}$, with $j$ the sum of $36/h$ signs. The same allowable sums follow if, independently for every slice and every block, the affine center of the scalar is different; the argument requires only those row orbit sizes, not that translating the slices preserves the correlation. The single scalar subgroup of each order lies in the cyclic group $\mathbb F_{37}^{\times}$.

**All divisors $h≥6$ of 36.** The slice sign sum has $j\equiv36/h\pmod2$, so $(j+j')/2$ is integral for each pair of slices. Thus $\gamma=(z+z')/2$ always belongs to $h\mathbb Z+\{-1,0,1\}$. A light value lies in $\{-1,0,1\}$; every other value has absolute value at least $h-1\ge5$. A row of eight that mixes light and heavy has at least seven cross pairs, each differing by at least four, so its contribution $8t_f-s_f^2$ is at least 112, above the total 49. An entirely heavy row exceeds the total norm 81. All rows must then be light, putting the total norm at most 32, a contradiction. This single argument excludes 6, 9, 12, 18 and 36.

For an independent cross-check using only the first compression, the following two cases also apply.

**Order 9.** Here $z=c+18v$, $v\in\{-2,-1,0,1,2\}$. Take $m=\sum v^2$, $r=\sum cv$ over all 64 slice rows. By the first equation of (2), $9m+r=36$. Since $|r|\le\sum|v|\le m$, $m=4$ and $r=0$. A single term of absolute value 2 would give $r=\pm2$, impossible. Thus there are exactly four nonzero terms, all of absolute value 1. The number of $e$-pairs with exactly one such term is even, and each contributes at least $16^2$ to the second equation of (2), so there are none. Two equal nonzero terms in the same pair would contribute at least $34^2>324$; therefore both terms vanish or have opposite signs. Then every pair contributes at most $2^2$, for a total at most $128<324$.

**Order 12.** There are three length-12 orbits outside the center; their sign sum is odd, so $|z|\ge11$ on all 64 slices. This contradicts $64\cdot11^2>1360$. A subgroup containing either order 9 or order 12 is also excluded; this handles orders 18 and 36.

**Order 4.** There are nine length-4 orbits, so $z\equiv4\pm1\bmod8$. Pairing and then halving gives $\gamma\equiv-1,0,1\bmod4$. For such an integer, $\gamma^2\ge\pm\gamma$ and $\gamma^2\ge\pm4\gamma-3$: the only integer violating the latter positive-sign inequality lies strictly between 1 and 3, hence is 2, and negation gives the other case. Summing eight scalar inequalities gives
$$
 t_f\ge \max(|s_f|,4|s_f|-24).
$$
As $\sum_fs_f^2=599$, the modest complete list of 14 nonnegative sorted four-square representations of 599 suffices. The minimum of $\sum_f\max(|s_f|,4|s_f|-24)$ is 82, attained at $(3,5,6,23)$, contradicting 81. For transparency, the full list and resulting bounds are below. It is generated by trying $0\le a\le b\le c\le24$ and checking whether $599-a^2-b^2-c^2$ is a square $d^2$ with $d\ge c$.

| Sorted absolute values | Lower bound for $\sum t_f$ |
|---|---:|
| $(1,6,11,21)$ | 87 |
| $(1,7,15,18)$ | 92 |
| $(2,3,15,19)$ | 93 |
| $(2,9,15,17)$ | 94 |
| $(3,5,6,23)$ | 82 |
| $(3,5,9,22)$ | 84 |
| $(3,7,10,21)$ | 86 |
| $(3,13,14,15)$ | 99 |
| $(5,5,15,18)$ | 94 |
| $(5,9,13,18)$ | 93 |
| $(6,7,15,17)$ | 93 |
| $(6,9,11,19)$ | 90 |
| $(6,13,13,15)$ | 98 |
| $(7,10,15,15)$ | 95 |

These exclusions leave orders $1,2,3$. Under the pure action alone, the strongest remaining case $h=3$ has $1+36/3=13$ orbits per slice and thus $4\cdot16\cdot13=832$ sign-orbit positions before any other constraints or symmetries. The result does not bound all automorphisms of a core or matrix; an action moving $A$ or permuting blocks need not contain the required pure slice actions. A signed affine action with a constant sign on a slice cannot have sign $-1$, since the affine center is fixed and its entry is nonzero. No claim is made for arbitrary position-dependent sign switches.

## Exact replay and scope

    python3 -B verify.py | diff -u replay.txt -
    sha256sum -c SHA256SUMS

The original integer-only verifier reruns the four-square list, checks the pair norms by residual-budget recursion for orders 9, 12, 18, and 36, and verifies that the pair norms alone do **not** exclude 4 or 6. A small generated width-28 core on $C_8$ at target order 60 has all seven off-subgroup correlations equal to 4. A separate exact witness after compression to $A=C_2\times(\mathbb F_8,+)$, $e=(1,0)$, satisfies all 16 compressed correlations and the slice orbit-sum rules simultaneously for $h=1,2,3$. It retains Singer scalar invariance on $\mathbb F_8$; it is **not** a full-core lift.

Primary status: [Egan et al., Table 5](https://arxiv.org/html/2610.27175v2) and [Munemasa–Núñez Ponasso, introduction](https://arxiv.org/pdf/2611.27591) identify 2396 as the first unresolved order. The [Davis–Astra 2396 report](https://raw.githubusercontent.com/1337isnot1337/solved-conjectures/main/results/hadamard-2396/README.md) already proves a different global signed order-479 automorphism restriction, discusses a width-4 order-11 multiplier obstruction on length 598, width-12 conditional borders and multiplier bounds on length 596, and length-1197 Legendre-pair multiplier constraints. None is the length-592, pure-$\mathbb F_{37}$, width-28 bound proved here. The [Ramos report](https://github.com/Arthur742Ramos/hadamard-668-multiplier-obstructions) also studies related scalar restrictions over 37 for length-333 Legendre pairs toward order 668, not this four-row template. Fourier restriction, character compression and orbit congruences are standard; the present contribution is the uniform arithmetic gate for this particular 2396 route. No later matching result found in a bounded literature check; no priority or general nonexistence claim.
