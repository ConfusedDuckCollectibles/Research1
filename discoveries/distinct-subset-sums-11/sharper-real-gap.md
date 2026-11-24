# A one-fifth real gap above the subset-sum bandwidth bound

Let $f_{\mathbb R}(n)$ be the least possible largest element of $n$ positive real numbers whose $2^n$ subset sums are pairwise at least 1 apart. Write
$$
 H_n=\sum_{i=0}^{n-1}\binom{i}{\lfloor i/2\rfloor}.
$$
Then $f_{\mathbb R}(n)\ge H_n+1/5$ for every $n\ge6$, with $f_{\mathbb R}(9)\ge148+1/4$.

This addendum strengthens the fraction in the [earlier bandwidth-gap argument](bandwidth-gap.md), using its saturated cuts and eleven-coordinate transfer. The integer corollary $f_{\mathbb Z}(n)\ge H_n+1$ is unchanged. This is partial progress on [Notik's full-unit conjecture](https://arxiv.org/pdf/2611.17360), §7, item 4, not a proof of it, an exact value at $n\ge9$, or an asymptotic improvement. The bandwidth benchmark $H_n$ comes from [Costa and Della Fiore](https://arxiv.org/abs/2609.27941). This AI-assisted proof and finite replay are not peer reviewed; no priority claim is made.

## 1. Rank and rounding inputs

Order the actual subset sums $t_r$ by increasing value, starting at rank 1. Put $H=H_n$ and $a_n=H+\epsilon$, where $a_n$ is the largest weight. The bandwidth bound gives $\epsilon\ge0$. We will rule out $\epsilon<1/5$, or $\epsilon<1/4$ at $n=9$.

A *saturated cut* $p$ is a canonical simplicial prefix whose exterior vertex boundary in the cube has size $H$. Such prefixes end at sets not containing the largest coordinate, taking one of each pair $(n-2,n-1),(n-4,n-3),\ldots$, and optionally coordinate 1 if $n$ is even. By Harper, the corresponding actual rank prefix has at least $H$ exterior neighbors. All have sums in $(t_p,t_p+H+\epsilon]$, which contains at most $H$ sums above $t_p$ when $\epsilon<1$. Hence they are exactly the next $H$ ranks. If $L,U$ are the first and last saturated cuts, $U-L<H$. On the band $[L,U+H]$, write
$$
 t_r=u+r+\theta_r.
$$
The excess $\theta_r$ is nondecreasing and, choosing $\theta_L=0$, lies in $[0,2\epsilon]$: the extreme saturated intervals overlap and each has excess rise at most $\epsilon$. Within **any** saturated interval $[p,p+H]$, every pair's excess difference is at most $\epsilon$.

Use the exact support-prefix counts at low and shifted cuts from §§2–3 of the earlier proof. They follow from Harper equality at all radii via [Räty, Theorem 4](https://arxiv.org/pdf/1806.11061v2) and binomial inversion. For a cut $p$, if $C_j(p)$ counts subsets of $[j]$ among the first $p$ ranks, the relevant identity is
$$
 C_j(p+H)-C_{j-1}(p+H)=C_{j-1}(p).
$$
In the canonical order this follows by deleting $j$ from the sets of $[j]$ containing it; the rank profiles transfer it to the actual order. A cut pair $p<q$ with $C_{j-1}(q)>C_{j-1}(p)$ supplies an actual such subset of rank $r\in(p,q]$. Adding weight $a_j$ sends it to a rank $r'\in(p+H,q+H]$: the identity gives the lower endpoint, and $t_r+a_j\le t_{q+H}+\epsilon$, with $\epsilon<1$, gives the upper endpoint by separation. Thus $a_j=b_j+\delta_j$ with
$$
 H-(q-p-1)\le b_j=r'-r\le H,\qquad 0\le\delta_j\le2\epsilon.
$$
Since $2\epsilon<1$, $b_j$ is the floor. The largest floor is $H$, with error $\epsilon$. The floors strictly increase by singleton separation. The last five floor deficits from $H$ are therefore $(a,b,c,1,0)$, with $2\le c\le5$, $c<b\le20$, and $b<a\le76$; in dimension 9 use just $(b,c,1,0)$. The verifier finds the necessary cut pairs by brute force in each reference cube; for $n\ge11$ use the first 32 low cuts, keeping the leading coordinates fixed and transferring the lexicographic order in the last eleven coordinates. The difference conditions use only count changes inside those cut intervals and transfer with that order. Finally, since $b_{n-1}=H-1$, singleton separation gives $\delta_{n-1}\le\epsilon$.

## 2. A signed bound only at queried ranks

Let $k=4$ in dimension 9 and $k=5$ otherwise. For an upper mask $P$ let $e_P$ sum the integer floors and $\delta_P$ sum their errors. Assign weights $2,\ldots,2,1,0$ to the $k$ bit positions, denoting their sum over any mask by $W(P)$, and let $f=1_{\{k\in Q\}}-1_{\{k\in P\}}$. With
$$
 A=W(Q\setminus P)+f,\qquad B=W(P\setminus Q)-f,
$$
the signed error satisfies $-B\epsilon\le\delta_Q-\delta_P\le A\epsilon$, including $\epsilon=0$. Its absolute value is at most $(2k-2)\epsilon$.

Retain only the ranks in the windows below. If $r,s$ are two retained ranks and $d=s-r$, take $\beta=0$ when $d=0$, $\beta=1$ if both ranks share a saturated interval, and $\beta=2$ otherwise. Monotonicity bounds the absolute value of
$$
 (\delta_Q-\delta_P)-(\theta_s-\theta_r)
$$
by $\epsilon\max(A,B+\beta)$ for $d>0$, $\epsilon\max(A+\beta,B)$ for $d<0$, and $\epsilon\max(A,B)$ for $d=0$. We need the bound only for **near shifts** $d-(e_Q-e_P)\in\{-1,0,1\}$. Its coefficient there is at most 4 in dimension 9 and 5 otherwise.

For $n=9,10$ the replay checks this claim by enumerating all ordered mask pairs, including identical masks. For each absolute rank difference it uses the worst $\beta$ among all retained pairs of that difference. In dimension 10 it is enough to omit the first three low sites. For $n\ge11$ the two rank windows have offsets $[0,238]$ and $H+[0,175]$ relative to $L$. The low window is contained in $[L,L+H]$. If $p=L+175$ is the last cut of the local group, the high window is contained in $[p,p+H]$. Within a window $H\ge526$ and the total deficit at most 102 force the same cardinality for near-shift masks: unequal cardinalities give a difference of size at least $H-102>239$. Both exclusive parts have at most two bits, so $A,B\le4$. With $\beta\le1$ this suffices. Across windows orient $r$ low and $s$ high. Then $d\in[H-238,H+175]$ forces $|Q|-|P|=1$ at a near shift. If $Q$ gains at most two exclusive bits, $A\le4$ and the bound $B+2\le5$ follows because it loses at most one. In the remaining case $Q$ gains three and loses two, using all five exclusive bits. If it gains the last bit then $A\le4+1$, $B\le4-1$; if it loses that bit then $A\le6-1$, $B\le2+1$. This proves the coefficient bound, with reversed orientation handled by swapping masks and ranks.

## 3. Exact cover after rounding

Each retained rank $r$ comes from one lower subset of sum $x$ and one upper mask $P$. Give $x$ the integer tag $z=r-e_P$. If the same lower subset also participates at $s$ with mask $Q$, the difference of the two tags is
$$
 d-(e_Q-e_P)=(\delta_Q-\delta_P)-(\theta_s-\theta_r).
$$
The coarse bound $2k\epsilon<2$ puts it among $-1,0,1$. The near-shift bound now makes its magnitude strictly below 1, so the two tags agree. If two different lower subsets get the same tag, the exact-shift bound makes the difference of their actual sums below 1, also impossible. If a translate $z+e_Q$ lands at a retained rank $s$, the same exact-shift bound puts the translated actual sum less than 1 from $t_s$. It must be that actual rank and must obey every applicable mask ban. In particular different translates cannot hit the same retained rank. Translated collisions **outside** retained ranks are not prohibited by this rounding argument, and are not tested.

Two distinct upper masks cannot have equal $e_P$ in a genuine example. Since $H$ exceeds the sum of the deficits in every template, equality would mean the same cardinality, leaving at most two exclusive bits on each side. Their actual sums would differ by at most $4\epsilon<1$. This excludes the duplicate-floor tails before search.

For the remaining tails, try all tags $r-e_P$ obtained from retained sites $r$. Reject a tag if one of its translates lands on a retained site with a banned upper mask. The bans are the earlier exact rank-profile restrictions: on consecutive-cut intervals $(p,q]$, the count of masks supported on an initial coordinate prefix is prescribed, and counts of zero or the whole interval forbid the corresponding class. Complementation reverses ranks by $r\mapsto 2^n+1-r$. The cut interval $(2^n+1-q,2^n+1-p]$ therefore maps to the local interval $[p,q)$; the transformed mask predicate asks whether the mask contains the tail suffix. The checker regenerates both restrictions in the template, not any conjectured reflection of the actual tag set. The remaining tags are rows on their sets of retained-site incidences. A true instance supplies an exact cover. Exhaustive row choice at a least-supported uncovered site, deleting only rows sharing one of its sites, shows no such cover exists.

The verifier uses the native $H$ in dimensions 9 and 10 and $H_{11}=526$ in the universal template. For the latter every tag has form $vH+r$ with residue $0\le r\le340$. Translating by an upper floor sum leaves the residue in $[-102,340]$, a range of width 442. Therefore comparisons with a retained rank or another floor sum agree at every $H\ge526$ with those in dimension 11. No computational search in growing dimension is used.

| Template | $H$ in replay | Retained rank offsets from $L$ | Duplicate floors | Forced rejection | Branched rejection |
|:--|--:|:--|--:|--:|--:|
| 9 | 148 | $[0,112]\cup(H+[0,49])$ | 4 | 59 | 3 |
| 10 | 274 | $[3,238]\cup(H+[0,175])$ | 455 | 3,493 | 262 |
| 11 and above | 526 | $[0,238]\cup(H+[0,175])$ | 455 | 3,487 | 268 |

All 8,486 template-tail cases are excluded. The separate untrimmed dimension-10 arithmetic control reaches coefficient 6 on ten tails, at absolute shifts 447, 448 and 449; the checked coefficient is 5 after trimming, not before. [Notik, Theorem 5.1](https://arxiv.org/pdf/2611.17360) gives the exact real values 24, 44, 84 for $n=6,7,8$, respectively, versus $H_n=23,43,78$. These finish the all-dimension statement.

## Reproduce

```sh
c++ -O3 -std=c++17 -Wall -Wextra -Wpedantic replay-sharper.cpp -o /tmp/subset-gap-replay
/tmp/subset-gap-replay
```

The code needs no inputs and uses only exact integer arithmetic and the C++ standard library. The C++ replay was independently implemented to rederive and check the rounding and cover calculation.
