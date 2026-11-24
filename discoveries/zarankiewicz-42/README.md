# At most 283 ones without a rectangle

Take a 42 × 42 grid and mark some cells. You may not mark all four corners of any rectangle whose sides follow the grid. How many marks can you place? The previous upper bound we found in the literature was **287**. This note rules out **284 and above**, so **283** is the new ceiling. A published 43 × 43 example gives **279** marks after deleting one row and one column, leaving the exact value open between 279 and 283.

In extremal graph notation this is the unconditional bound **$z(42,42;2,2)\le 283$**. The new reduction is a counting proof. Its only computational dependency is the known bound $\mathrm{pa}(5;6)\le31$: at most 31 length-five words on a six-symbol alphabet can pairwise agree in at most one position. An independently written, self-contained checker for that dependency is included.

## Proof

Work with the equivalent incidence structure, with 42 points and 42 lines, distinct lines intersecting in at most one point. It suffices to exclude exactly $E=284$ incidences, since any larger example could be thinned to 284. For each line $M$, let $a_M=7-|M|$. The signed defect sum is $D=\sum_M a_M=294-E=10$. Let $L$ be the number of uncovered unordered point pairs, and set $X=\sum_M a_M(a_M-1)/2$. The polynomial convention matters for negative integers: every term is nonnegative, and a line of size eight contributes 1. From
$$
  \binom{7-a}{2}=21-6a+\frac{a(a-1)}2
$$
and the fact that no pair is covered twice,
$$
  L+X=\binom{42}{2}-42\cdot21+6D=39.
$$

**High line exclusion.** Fix a line $B$ of size $s$, let $w$ be the number of other lines disjoint from $B$, and let $L_{\rm cross}$ count uncovered pairs with one point in $B$. Each other line meeting $B$ covers exactly its size minus one cross pairs, and no pair is counted twice. Thus
$$
  L_{\rm cross}=s(43-s)-(E-41)+6w-\sum_{M\cap B=\varnothing}a_M,
$$
where the sum is over **other lines** (in this argument $s\ge8$). Also $a(a-1)/2\ge a-1$ for every integer $a$. Applying that inequality just to the disjoint lines, and retaining the term for $B$, yields
$$
  L+X\ge\ s(43-s)-(E-41)+\frac{(7-s)(6-s)}2+5w.
$$
At $E=284$, the right side without $5w$ is the concave quadratic $s(73-s)/2-222$; on $9\le s\le42$ its endpoint minimum is 66. Thus $s\le8$ for every line. If $s=8$, the bound is $38+5w\le39$, so $w=0$. Then $L_{\rm cross}=37$, $1\le X\le2$, and there is at most one other line of size eight. There is at most one uncovered pair among the 34 points outside $B$, since $L\le38$. Each such point is therefore collinear with at least 32 other outside points. All other lines meet $B$ exactly once. An outside point on at most six lines could see at most $6\cdot5+1=31$ other outside points, even allowing the one other size-eight line. Thus every outside point has degree at least seven. The eight points of $B$ have total degree $8+41=49$, hence their signed defect sums to 7. The outside points contribute nonpositive signed defect. This contradicts the total point-side signed defect 10. Transposing also excludes points of degree eight; both line sizes and point degrees are at most seven.

**Clean point and packing.** There are at least 32 points of degree seven, since the total nonnegative point defects sum to 10. For a point $p$ of degree seven, its number of uncovered pairs is
$$
  \delta(p)=41-\sum_{M\ni p}(|M|-1)=\sum_{M\ni p}a_M-1.
$$
If every degree-seven point had $\delta(p)\ge1$, there would be at least $2\cdot32=64$ weighted incidences with short lines. The total of all these weights over all points is at most $\sum_M (7-a_M)a_M\le6\sum_M a_M=60$, impossible. Some degree-seven point has no uncovered pair and lies on one six-point line and six seven-point lines. Removing that point from these seven lines partitions the other 41 points into one group of five and six groups of six. Every one of the 35 other lines meets each group in at most one point. The seven lines through the chosen point contribute $6+6\cdot7=48$ incidences; the 35 off-point lines contribute 236, so there are $35\cdot7-236=9$ missing group-line intersections. The group of five is met at most $5\cdot6=30$ times, since each of its points has degree at most seven and is already on a line through the chosen point. Therefore at most four missed intersections occur among the six groups of six. Discard the group of six with the most misses; the other five have at most three between them (if there are no misses this also holds). At least $35-3=32$ lines meet every one of the five retained groups. Recording their chosen symbol in each group gives 32 distinct words; two words cannot agree in two positions, since that would be two intersecting points on their lines. This contradicts the prior packing bound. ∎

## Reproduce the finite dependency

The proof above takes the hexary packing ceiling from [Sadhu, Theorem 6.4 and §12](https://arxiv.org/html/2608.01606v2). The supplied C++20 checker independently rules out 32 words, uses no external data or solver, and needs a compiler supporting unsigned 128-bit integers (GCC or Clang). Do not define `NDEBUG`.

```sh
g++ -O3 -std=c++20 check-pa56.cpp -o check-pa56
./check-pa56
```

It enumerates 393,120 possible second full classes, groups them into 56 orbits using simultaneous symbol permutations and coordinate permutations, and refutes every size profile of the four remaining classes. The five profiles exhaust nonincreasing tuples of four counts between 0 and 6 summing to 20. There must be at least two full classes in a 32-word candidate because one full class and five smaller ones give at most 31. Normalize the first full class by relabeling the four tail coordinates so its words have tails `(j,j,j,j)`. Every other tail then has four distinct symbols. Within a remaining class tails must disagree everywhere; across classes they can agree at most once. The checker generates completions subject to these constraints. Within one class it tries tail indices in increasing order; equal-sized classes can be ordered by their smallest tail index because those indices are distinct. The final line should be `PASS: no 32-word packing; pa(5,6)<=31`.

## Sources, lower bound, and scope

- The previous ceiling of 287 used here follows from [Damásdi, Héger and Szőnyi, Corollary 3.19, Eq. (3.3)](https://real.mtak.hu/13977/1/DHSzFIN.pdf): at $n=c=6$ the ceiling is 288, with equality forcing an embedding in a projective plane of order six. Such a plane does not exist. [Dybizbański, Dzido and Radziszowski, Theorem 3](https://cs.rit.edu/~spr/PUBL/zr14.pdf) gives the weaker ceiling of 288 directly.
- [Sadhu's §13](https://arxiv.org/html/2608.01606v2#S13) explicitly identifies improving the 42 × 42 bound as an open direction, using 288 as its reference. The packing-array input and clean-point strategy are credited to Sadhu. The improved bound here does not automatically give an improvement for 43 × 43.
- For the lower bound, use [Sadhu's Appendix A](https://arxiv.org/html/2608.01606v2#A1) and delete the first listed line and point 12 (1-based). Both have six incidences and their shared cell is marked; the result has $290-6-6+1=279$ marks. The original 290-one construction is credited in Sadhu's paper to MathOverflow user user48028. No third-party matrix, checker, or article text is bundled.

This note was prepared with AI assistance and independently checked within the project. It has not undergone external peer review.
