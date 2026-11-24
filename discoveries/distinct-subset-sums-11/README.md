# A lower bound of 527 for eleven distinct subset sums

Pick eleven distinct positive integers. If all 2,048 subsets must have different sums, **at least one of the integers must be 527 or larger**. The previous lower bound derived from the hypercube bandwidth was 526. A known construction reaches 594, so the exact minimum is still open.

In the notation of the Erdős distinct subset sums problem, this gives $527\le f_{\mathbb Z}(11)\le594$. The result is global over positive integers, not a test limited to modifying a known construction. It does not prove a bound of 527 for the real-valued relaxation or an asymptotic result.

## 1. What equality at 526 would require

Suppose the entries are $a_1<\cdots<a_{11}\le526$. Rank the 2,048 subsets by increasing sum, starting at 1. Let $F_m$ be the first $m$ subsets and $t_m$ their last sum. Let $I_m$ instead be the first $m$ subsets in the simplicial order, comparing first by cardinality and then lexicographically by the increasing index lists. On the 11-cube, let $N^d$ denote closed radius-$d$ neighborhoods.

[Harper's theorem](https://doi.org/10.1016/S0021-9800(66)80059-5) says $I_m$ minimizes the size of $N^1$ among families of size $m$. Counting the cube gives 32 cuts $m$ with $|N^1(I_m)|-m=526$. Their extremes are 674 and 849. An outside neighbor of $F_m$ arises by adding one entry, so its sum lies between $t_m+1$ and $t_m+a_{11}$. There must be at least 526 such neighbors with distinct integer sums. Thus $a_{11}=526$, the entire interval is filled, and $N^1(F_m)=F_{m+526}$. Overlap between the extreme intervals forces 702 consecutive values at ranks 674 through 1375. Write those values as $q+u$, with rank $q$ and integer offset $u$. Complementation reverses ranks; since $674+1375=2049$, the entries sum to $2u+2049$.

For a family $F$ defined by a threshold on positive increasing weights, let $h_j(F)$ count its subsets whose largest index is $j$, with index zero for the empty subset. For $0\le d\le11$,

$$
|N^d(F)|=\sum_{i<d}\binom{11}{i}+\sum_{j=0}^{11}h_j(F)\binom{11-j}{d}.
$$

To test whether a set belongs to the neighborhood, remove its $d$ largest indices, or all its indices if it has fewer. This produces the least-weight set within distance $d$. The sets of size less than $d$ give the first sum; for a set in $F$ with largest index $j$, the deletion has $\binom{11-j}{d}$ inverses among larger sets. We use this identity for nonempty families containing the empty set. It also applies to $I_m$, realized by the increasing synthetic weights $K-2^{11-i}$ with $K>2^{11}$. The binomial matrix becomes triangular when the order of the indices is reversed, so the neighborhood sizes determine all the $h_j$.

[Räty's Theorem 4](https://arxiv.org/pdf/1806.11061v2) says a family minimizing the one-step neighborhood also minimizes it at all radii. At the 32 saturated cuts, therefore, the histogram of $F_m$ is that of $I_m$. The same is true of $F_{m+526}$: taking a neighborhood of its one-step predecessor adds one to the radius, and $N^1(I_m)=I_{m+526}$, as checked by the replay. At all 64 cuts this determines the cumulative count $C_j(q)$ of subsets of the first $j$ entries with sum at most $q+u$.

The identity $C_j(p)-C_{j-1}(p)=C_{j-1}(p-a_j)$, together with monotonicity, gives $a_j\ge p-q+1$ whenever the left side is smaller than $C_{j-1}(q)$. Here are four comparisons fully regenerated from the cube:

| $j$ | $p$ | $q$ | Difference at $p$ | $C_{j-1}(q)$ | Lower bound |
|---:|---:|---:|---:|---:|---:|
| 7 | 1249 | 800 | 62 | 63 | 450 |
| 8 | 1340 | 835 | 118 | 119 | 506 |
| 9 | 1365 | 845 | 211 | 212 | 521 |
| 10 | 1372 | 848 | 354 | 355 | 525 |

There are exactly 4,210 increasing tails $H=(a,b,c,525,526)$ with $a\ge450$, $b\ge506$, and $c\ge521$.

## 2. Exhaustive obstruction to covering the forced sums

For each such $H$, form its 32 subset sums $E$. Any repeat excludes the tail. Let $X$ be the six smaller entries' subset sums translated by subtracting $u$, and put $\tau=2049-\sum H$. A genuine $X$ is invariant under reflection $x\mapsto\tau-x$ and contains no fixed point of that reflection: complementary subsets of the six entries are different. The six entries sum to at most $B=6a-21$, so every $x\in X$ has $|2x-\tau|\le B$. Every central rank 674 through 1375 must occur in $X+E$, and all representations must be unique, including outside the central interval.

For $6\le j\le11$, let $E_j$ be the subset sums using only the first $j-6$ positions of the tail. The known count difference $C_j(q)-C_j(p)$ at consecutive known cuts $p<q$ counts the values of $X+E_j$ in $(p,q]$. A zero difference bans sums in that interval with $e\in E_j$; a difference of $q-p$ bans sums with $e\notin E_j$, since the interval is already completely filled by the prefix. We include fixed zero and full-count cuts outside the allowed support.

Discard individual positions violating those bans and reflection pairs with internal collisions. For any as-yet uncovered central rank, select one of the live reflection pairs that can cover it. No remaining pair means contradiction, one is forced, and with several the search tries every one. Selecting a pair deletes every other pair whose translates collide with it. This is a relaxation: we do not require the surviving positions to factor into six entries' subset sums. Every genuine solution would nevertheless survive some branch at each step.

The replay excludes 455 tails by repeated upper sums, 3,639 with only forced choices, and the other 116 by exhaustive branching. There are 265 branch edges in total, with at most six for any tail. Every branch ultimately has a missing central rank. There is no surviving candidate, so the integer maximum cannot be 526 or less.

## Reproduce

With Python 3.10 or newer and no external input:

```sh
python3 occupancy-cover.py
python3 verify-occupancy.py occupancy-full.cert.jsonl
```

The second program reconstructs the cube and geometric constraints independently of the first, checks all certificate branches, and tests known positive cases in dimensions 3, 4 and 5. It rejects an intentionally false certificate in each control. Do not run Python with `-O`.

An independent C++17 search needs no certificate file or external data:

```sh
g++ -std=c++17 -O3 -Wall -Wextra independent.cpp -o independent
./independent
```

## Prior work and scope

[Costa and Della Fiore, Theorem 1.1](https://arxiv.org/abs/2609.27941) give the bandwidth lower bound, which evaluates to 526 at eleven. The upper endpoint 594 comes from the [Conway–Guy construction, proved by Bohman](https://doi.org/10.1090/S0002-9939-96-03653-2): `285, 433, 510, 550, 570, 581, 587, 590, 592, 593, 594`. [Dyson's exhaustive search](https://github.com/pwdyson/erdos_1) gives the exact integer value at ten; the eleven-entry case is still open. [Notik's Problem 4](https://arxiv.org/abs/2611.17360) concerns a related conjecture over the reals, to which this integer cover argument does not apply.

This is an AI-assisted research note with independent proof replays, not external peer review or a claim of exhaustive priority.
