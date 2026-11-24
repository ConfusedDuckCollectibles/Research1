# A uniform gap above the subset-sum bandwidth bound

Choose numbers so that every subset has a different sum. An earlier formula gave a lower bound for the largest number you need. We prove that the bound cannot be attained once you choose six or more numbers, even if you allow real numbers and require the sums to stay at least one unit apart. The argument for arbitrarily many numbers reduces to one fixed eleven-coordinate check. This extends the [earlier eleven-integer proof](README.md).

Let $f_{\mathbb R}(n)$ be the least possible largest element of an
$n$-element set of positive reals whose $2^n$ subset sums are pairwise
at least 1 apart. Let $f_{\mathbb Z}(n)$ use positive integers, and put
$H_n=\sum_{d=0}^{n-1}\binom d{\lfloor d/2\rfloor}$.

**Theorem.** For every $n\ge6$,
$f_{\mathbb R}(n)\ge H_n+1/20$ and $f_{\mathbb Z}(n)\ge H_n+1$.
A separate finite check gives the slightly larger real gap $1/16$
at $n=9$. The real results use an 11-coordinate suffix for all
$n\ge11$, finite covers for $n=9,10$, and known exact values for
$n=6,7,8$. In particular $f_{\mathbb Z}(11)\ge527$,
$f_{\mathbb Z}(12)\ge989$ and $f_{\mathbb Z}(13)\ge1913$.
No unsettled exact value or proposed real gap of a full unit is
proved. This AI-assisted work is not peer reviewed, and makes no
priority claim.

## 1. Saturated cuts are pair transversals

Let $r_n(S)$ denote the rank, starting at 1, of $S\subseteq[n]$ in the
simplicial order (cardinality, then lexicographic order of increasing
index lists). Write $I_m$ for its first $m$ sets. The smallest
$k$-subset of a $(k+1)$-set deletes its largest index. It follows
that if the last $k$-set $S$ of $I_m$ excludes $n$, then the closed
one-step neighborhood $N(I_m)$ is the initial segment ending at
$S\cup\{n\}$. The lower layers are already in $N(I_m)$, and in layer
$k+1$ the largest admissible extension appends $n$ to $S$.

For $S\subseteq[n-1]$ define
$D_n(S)=r_n(S\cup\{n\})-r_n(S)$. The within-layer lexicographic
order splits into the block containing index 1 followed by the block
omitting 1. In both cases deleting the first index from the ground
set gives
$D_n(S)=\binom{n-1}{|S|}+D_{n-1}(S\setminus\{1\})$, with
suffix indices relabeled. If $1\in S$, subtract the lower-dimensional
layer increment $\binom{n-1}{|S|-1}$ from $\binom n{|S|}$.
If $1\notin S$, the offsets of the preceding 1-blocks differ by
$\binom{n-1}{|S|}-\binom{n-1}{|S|-1}$, and subtract the lower
increment $\binom{n-1}{|S|}$. Both yield the displayed recursion.
With the one-coordinate base case this gives

$$
D_n(S)=\sum_{i=1}^n\binom{n-i}{|S\cap[i,n-1]|}\le H_n.
$$

If $S$ contains $n$, any last neighbor in layer $k+1$ has a lex-first
parent strictly earlier than $S$: that parent cannot contain $n$.
Appending $n$ to the parent only increases the neighbor's rank;
the bound just proved shows the resulting gap from $S$ is strictly
less than $H_n$. If no next-layer neighbor exists, the gap is also
strictly smaller. Thus $|N(I_m)|-m=H_n$ occurs exactly when every
term in the sum for $D_n(S)$ is central. These are the sets excluding
$n$ and containing exactly one of each pair
$(n-2,n-1),(n-4,n-3),\ldots$, with optional index 1 when $n$ is
even. Call such a cut saturated.

Let $L$ and $U$ be the smallest and largest saturated cuts. Changing
the choice in the $i$th pair from the right changes the rank by
$\binom{2i-1}{i-1}$. Indeed, with all other choices fixed, the
first-index lex blocks give
$\binom{2i}{i-1}-\binom{2i-1}{i-2}$, with a negative lower index
interpreted as zero, giving this value by Pascal. When $n=2r$,
adding the optional index 1 changes layer from $r-1$ to $r$, giving
$\binom n{r-1}-\binom{n-1}{r-2}=\binom{n-1}{r-1}$. Consequently

$$
U-L=W_n:=\sum_{i=1}^{\lfloor n/2\rfloor}\binom{2i-1}{i-1}<H_n.
$$

The saturated vertices are invariant under complement in $[n-1]$.
Full complementation reverses simplicial rank and exchanges the
set of saturated cuts with the same cuts shifted by $H_n$.
In particular $L+U+H_n=2^n+1$.

## 2. Real near-saturation yields the exact rank profiles

Suppose the weights are $a_1<\cdots<a_n$ and their subset sums are
1-separated. Let $t_r$ be the $r$th actual subset sum in increasing
order and let $F_r$ be the first $r$ subsets in that order. Harper's
vertex-isoperimetric theorem compares $N(F_r)$ to $N(I_r)$. At a
saturated cut $p$, each exterior neighbor is reached by adding a
weight and has sum in $(t_p,t_p+a_n]$. As there are at least $H_n$
neighbors separated by 1 and each is separated by 1 from $t_p$,
$a_n\ge H_n$.

Set $H=H_n$ and suppose $a_n=H+\epsilon$ for
$0\le\epsilon<1/20$. The same interval fits at most $H$ subset sums
above $t_p$, giving $N(F_p)=F_{p+H}$ and
$t_{p+H}-t_p\le H+\epsilon$. Thus every saturated cut attains
Harper equality. By [Räty, Theorem 4](https://arxiv.org/pdf/1806.11061v2),
a one-step minimizer is a minimizer at every neighborhood radius;
no assumption on the complement is needed. Since
$N^d(N(F_p))=N^{d+1}(F_p)$ and $N(I_p)=I_{p+H}$, the corresponding
equality also holds between $F_{p+H}$ and $I_{p+H}$ at every radius.

For a threshold family $F$ on the ordered coordinates, write $h_j(F)$
for the number of its vertices of largest index $j$, with the empty
set using index 0. Removing the $d$ largest elements of any vertex
minimizes its attainable sum within distance at most $d$. For a
fixed resulting set of maximum index $j$, there are
$\binom{n-j}{d}$ ways to restore $d$ larger indices, so

$$
|N^d(F)|=\sum_{i<d}\binom ni+\sum_{j=0}^n h_j(F)\binom{n-j}{d}.
$$

The formula also applies to $I_r$, realized by weights
$K-2^{n-i}$ with $K>2^n$. The binomial matrix is invertible.
Consequently at each of the stated cuts the exact number
$C_j(r)$ of actual subsets of $[j]$ at ranks at most $r$ equals
the corresponding count in $I_r$.

Write $t_r=r+u+\theta_r$ on the band $L\le r\le U+H$, with
$\theta_L=0$. Because of separation, $\theta_r$ is nondecreasing.
The two extreme rank intervals $[L,L+H]$ and $[U,U+H]$ overlap,
and each has total excess over unit spacing at most $\epsilon$.
Therefore $0\le\theta_r\le2\epsilon$ throughout the band.

## 3. Transfer from one fixed 11-coordinate suffix

For $n\ge11$, the first 32 saturated cuts choose all possibilities for
the five rightmost pairs. Their earlier coordinates are fixed; for
even $n$ the optional index 1 is absent throughout this group. These
cuts all lie in one fixed-cardinality layer. The cuts shifted by $H$
append $n$, and lie in the next layer with the same early prefix.

Between any two of these cuts in its own group, the early portion
and cardinality are fixed. Indeed the lex order on a layer groups
sets with identical membership on consecutively earliest coordinates.
Thus the intervening suffix masks, in order, are precisely those
between the respective cuts of the 11-cube, on the last 11 coordinates.
All within-group rank differences and count changes for support on
the first $j$ coordinates, $n-5\le j\le n$, are literally the
11-dimensional ones. At dimension 11, these two groups of cuts begin
at 674 and $674+526$. The range of offsets in either group is
$[0,175]$. The verifier constructs the 32 cut ranks directly from
the five pairs, checking the shifted ranks by appending coordinate 11.

For any of these first-group cuts $p$, and $j<n$, the **reference**
profiles satisfy

$$
C_j(p+H)-C_{j-1}(p+H)=C_{j-1}(p).
$$

The left side counts subsets of $[j]$ containing $j$ in $I_{p+H}$.
Deleting $j$ gives exactly the subsets of $[j-1]$ in $I_p$:
all lower layers are included, and in the top layer appending $j$
can be compared with the cut set with $n$ appended. A parent
equal to the cut set, if allowed in $[j-1]$, stays in the segment
because $j<n$. The measured actual profiles are the same. Writing
$B_{j-1}(x)$ for the number of actual subsets of $[j-1]$ of sum
at most $x$, the identity gives
$B_{j-1}(t_{p+H}-a_j)=B_{j-1}(t_p)$.

The following suffix ranks are given in the 11-dimensional
numbering; when transferring, preserve their offsets from the start
of the group. Each interval $(p,q]$ contains a set of maximum index
less than $j$ in the template, so $C_{j-1}(q)>C_{j-1}(p)$.

| $j$ in the 11-cube | $p$ | $q$ | $q-p-1$ |
|---:|---:|---:|---:|
| 7 | 723 | 800 | 76 |
| 8 | 814 | 835 | 20 |
| 9 | 839 | 845 | 5 |
| 10 | 846 | 848 | 1 |

Choose an actual subset of $[j-1]$ in $(p,q]$ of rank $s$.
The identity forces the sum after adjoining $j$ to exceed
$t_{p+H}$, so the new rank $s'$ is at least $p+H+1$. Both sums are
in the central band: $s\le q\le L+174\le U-1$, and separation
gives at least $H+1$ available between $t_s$ and $t_{U+H}$,
whereas $a_j<H+1$. We therefore have
$a_j=(s'-s)+(\theta_{s'}-\theta_s)$. With $b_j=\lfloor a_j\rfloor$,
this proves $b_j=s'-s\ge H-(q-p-1)$ and
$0\le a_j-b_j\le2\epsilon$, for the four largest weights below $a_n$.
For $j=n$, $b_n=H$ and $a_n-b_n=\epsilon$ by its assumed range.
Singleton sums are 1-separated, so the five integer floors are
strictly increasing. Their deficits from $H$ are therefore among the
same 4,210 decreasing tuples, bounded coordinatewise by
$(76,20,5,1,0)$, in every dimension $n\ge11$.

## 4. A robust, dimension-free exact cover

For a tail deficit tuple, let $e_M$ be the sum of the five integer
floors for mask $M\subseteq[5]$; write $\delta_M$ for the difference
between the actual upper sum and $e_M$. We have
$0\le\delta_M\le9\epsilon$, using $2\epsilon$ for each of four
weights and $\epsilon$ for the last one. Equal $e_M$ for two distinct
masks would contradict separation, so 455 tuples are immediately out.

The required full rank sites are **both** windows

$$
L+[0,238]\quad\hbox{and}\quad L+H+[0,175],
$$

all within the central band. Keep the lower subsets that take part in
at least one representation of a required rank and assign their lower
sum $x$ the nominal integer $z=r-e_M$ if the representation of
rank $r$ uses upper mask $M$. Then

$$
-9\epsilon\le x-u-z=\theta_r-\delta_M\le2\epsilon.
$$

Because the interval has width $11\epsilon<1$, the nominal integer
is uniquely defined, even if a lower subset occurs in more than one
representation. Distinct lower subsets have distinct tags by
separation. A full sum from nominal lower point $z$ and any upper
mask $M'$ has nominal coordinate $z+e_{M'}$ and error relative
to its nominal value plus $u$ in $[-9\epsilon,11\epsilon]$.
If the nominal coordinate is a central rank, its distance from
the corresponding actual $t_r$ is at most $11\epsilon<1$, so it
must be that actual subset sum. If two distinct representations
have the same nominal coordinate anywhere, their actual distance
is at most $20\epsilon<1$, again impossible.

To describe the finite problem, use integer coordinate pairs
$(k,r)$ to represent $kH+r$ relative to $L$. A tail translate is
$(|M|,-\sum_{i\in M}d_i)$ with the deficits $d_i$. Candidate lower
points come from subtracting one translate from one required point.
Their residues lie in $[0,340]$, and after adding any upper mask
their residues lie in $[-102,340]$. This whole interval has width
442, less than $H\ge526$. Thus coordinate equality, translation
collisions, and alignment with any required rank window are
independent of the actual $H$. The generator uses synthetic spacing
2048; the independent verifier uses literal coordinate pairs.

There are two types of bans on an upper mask at a nominal full sum,
using successive cuts $p<q$ within each group.

1. For each last-five prefix mask $g$, the count change for
   prefix-restricted subsets in $(p,q]$ equals the reference count
   of masks contained in $g$. A zero count forbids every such mask
   at these sites; a count of $q-p$ forbids every other mask.
2. The mirrored cuts $(2^n+1-q,2^n+1-p]$ are known cuts as well.
   Under full subset complementation this interval corresponds to
   **$[p,q)$**, including its left endpoint and excluding its right
   endpoint. Small-prefix membership of the complementary subset
   means the original mask contains all tail bits outside $g$.
   Counting such masks in the reference interval $[p,q)$ gives
   exactly the same zero/full bans for that predicate.

No whole-spectrum support bound, cardinality of the remaining low
set, symmetry of the selected lower points, or assumption on the
remaining values is used. The real configuration gives a set of
distinct nominal lower points satisfying all bans, with all its
translates disjoint and every required site covered.

## 5. Certificates and replay

The finite checker tests all 4,210 deficit tails. A permitted lower
point arises by subtracting one upper translate from one required
site, and is removed if any one of its full translates is banned.
One uncovered necessary site refutes the tail. A site with exactly
one candidate forces that point; every other point whose translates
would collide is removed. At a branch, the certificate covers
each candidate for its chosen site, in upper mask-index order.
Thus this checks a necessary exact cover, not a random search or
a numerical solver. No cutoff is treated as proof.

| Certificate type | Tails |
|---|---:|
| Repeated rounded upper subset sum | 455 |
| Contradiction with zero or more forced selections, no branch | 3,487 |
| Exhaustive branching with contradictory leaves | 268 |
| **Total** | **4,210** |

The checker verifies 738 branch edges, 4,139 contradiction leaves and
2,138 unit selections in total. The largest tree has 14 branch edges.
With Python 3.10 or newer
without `-O`, run:

```sh
python universal11-cover.py
python verify-local11.py universal11.cert.jsonl
```

Both use the standard library and compute the suffix template.
The verifier is separately written, using ordered tuples,
explicit mask sets and coordinate pairs rather than integer
bitsets. It checks a feasible control that must reject false
proofs, an infeasible single-row control, every tail and every
certificate leaf. The generated certificate need not be checked in. An independent C++17 checker derives the same obstruction without reading a certificate or using synthetic integer spacing:

```sh
g++ -std=c++17 -O3 replay-uniform.cpp -o /tmp/dss-uniform
/tmp/dss-uniform
```

All possible tails have contradictions, excluding
$a_n<H_n+1/20$. This proves the theorem for every $n\ge11$.

## 6. The remaining dimensions

The real cases $n=6,7,8$ have exact values $24,44,84$ by
[Notik, Theorem 5.1](https://arxiv.org/abs/2611.17360); the respective
bandwidth baselines are $23,43,78$. The same real rounding argument
applies to dimensions 9 and 10, with finite rank templates rather
than the uniform 11-coordinate suffix.

| Dimension | $H$ | Extreme saturated low cuts | Lower query offsets | Upper query offsets |
|---:|---:|:---:|:---:|:---:|
| 9 | 148 | 158–207 | $[0,112]$ | $H+[0,49]$ |
| 10 | 274 | 288–463 | $[0,238]$ | $H+[0,175]$ |

Use the absolute spacing $H$ in these two finite cases: their
residue windows need not have width less than $H$, so replacing
actual positions by formal coordinate pairs could miss collisions.
The grade-count witnesses, each in the low cuts, are:

| Dimension | Upper tail size | Final deficit bounds | Witness pairs $(p,q)$ |
|---:|---:|:---:|:---|
|9|4|$(20,5,1,0)$|$(172,193),(162,168),(159,161)$|
|10|5|$(76,20,5,1,0)$|$(337,414),(302,323),(292,298),(289,291)$|

For clarity, the rounding step at any such witness works even if
$q$ is the last saturated low cut. Adjoining $a_j$ to a prefix
subset of rank $s\in(p,q]$ gives a rank $s'>p+H$ by the cut identity.
Its value is at most $t_q+H+\epsilon\le t_{q+H}+\epsilon$.
Since distinct subset sums are 1-separated and $\epsilon<1$,
$s'\le q+H$. Thus both ranks lie in the band as required.
For a tail of size $k$, with last error $\epsilon$ and each of the
other $k-1$ errors at most $2\epsilon$, the total error is at most
$\Delta=(2k-1)\epsilon$. A nominal full sum then has error in
$[-\Delta,\Delta+2\epsilon]$, whose width is $4k\epsilon$.
Using $\epsilon<1/(4k)$ makes the same rounding and collision
arguments valid.

The finite replay counts follow:

| Dimension | Candidate tails | Repeated upper sum | Forced-only contradiction | Exhaustive branched contradiction | Branch edges |
|---:|---:|---:|---:|---:|---:|
|9|66|4|59|3|6|
|10|4,210|455|3,493|262|620|

Run `python real-small-replay.py` for the standalone exact-coordinate
checks and feasible/infeasible controls. In particular this proves
$f_{\mathbb R}(9)\ge148+1/16$ and $f_{\mathbb R}(10)\ge274+1/20$.
The real bound $H_n+1/20$ therefore holds for **every** $n\ge6$. A separate C++17 replay checks both exact-coordinate cases:

```sh
for n in 9 10; do
  g++ -std=c++17 -O3 -DDIM=$n replay-small.cpp -o /tmp/dss-small
  /tmp/dss-small
done
```

The previously reported finite upper witnesses still give
$f_{\mathbb Z}(12)\le1157$ and $f_{\mathbb Z}(13)\le2249$.

### Sources

* L. H. Harper, [Optimal numberings and isoperimetric problems on graphs](https://doi.org/10.1016/S0021-9800(66)80059-5), 1966.
* E. Räty, [Uniqueness in Harper's vertex-isoperimetric theorem](https://doi.org/10.1016/j.disc.2019.111696), 2020, [Theorem 4, arXiv v2](https://arxiv.org/pdf/1806.11061v2).
* S. Costa and S. Della Fiore, [A bandwidth refinement of the Erdős distinct subset sums bound](https://arxiv.org/abs/2609.27941), 2026, Theorem 1.1 gives $H_n$ (988 and 1912 in dimensions 12 and 13).
* E. Notik, [Positive dual certificates for real separated subset sums](https://arxiv.org/abs/2611.17360), 2026, §7, Problem 4 conjectures the real gap of at least 1 for all $n\ge6$; this note gives a smaller uniform positive gap, not that constant.
* [OEIS A276661, including explicit upper witnesses](https://oeis.org/A276661/internal); P. Dyson's [original n=10 computation](https://github.com/pwdyson/erdos_1); Borwein–Mossinghoff's [n=9 computation](https://doi.org/10.1090/S0025-5718-02-01460-6); [Lunnon's n≤8 computation](https://doi.org/10.1090/S0025-5718-1988-0917837-5).
