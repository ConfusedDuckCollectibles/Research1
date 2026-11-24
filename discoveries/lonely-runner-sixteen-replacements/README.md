# Arbitrary replacements of the standard sixteen-runner speed set

For a finite set V of positive integer moving speeds, write

    ML(V) = max_t min_{v in V} ||vt||,

where ||x|| is distance to the nearest integer. Sixteen total runners
means fifteen moving speeds and one stationary runner; its lonely-runner
threshold is 1/16.

**Theorem.** Delete any k-element subset of {1, …, 15} and add any k
distinct integers at least 16, where k is 2, 3, 4, or 5. The resulting set V
satisfies ML(V) > 1/16. There is **no upper bound on the added speeds**
in this statement. This is a restricted result, not a proof for all
sixteen-runner configurations and not a full classification of tight sets.

## Proof and exact certificate

Let n=16, α=2/n, and B_w={t: ||wt||≤1/n} be the closed bad set for speed
w. For any interval J of length L,

    |J ∩ B_w| ≤ α L + α(1−α)/w.                        (1)

In the scaled coordinate wt, the bad set is 1-periodic and has measure α
per period. Split the scaled interval into floor(wL) full periods plus
a remainder of length q in [0,1). Its bad measure is at most min(q,α):
modulo one, a remainder shorter than one period has no overlap with
itself except possibly at endpoints. The excess over αq is at most
α(1−α), as seen by splitting at q=α. Divide by w to obtain (1).

Fix a deleted subset and install the retained speeds. Insert the new
speeds in increasing order. At an intermediate prefix, suppose the
installed speeds all have gap strictly greater than 1/n throughout an
open interval J of length L>0. If r new speeds remain, with smallest u,
and the completed set had ML≤1/n, their bad sets would cover J.
Subadditivity in (1), using that all their speeds are at least u, gives

    L ≤ rα L + rα(1−α)/u,
    so u ≤ rα(1−α)/[(1−rα)L].                       (2)

Here n>2r in every use. This is a proved cutoff for *each* prefix,
not a selected numerical window. For each allowed integer u within
the inclusive floor, install it and proceed. An empty prefix interval
would stop the proof, not exclude a branch. At a complete branch, one
rational time with every gap strictly greater than 1/n suffices.
An alleged failing tuple must follow some enumerated branch at every
depth, by (2), but every terminal branch has a strict witness.

The C++ generator obtains one open interval by exact intersections
of rational good bands. The independent Python checker does not
recompute or intersect the good regions. For a proposed interval (a,b), it checks separately
for each installed speed w that there exists an integer j such that
n(wa−j)≥1 and n(wb−j)≤n−1. The function wt−j is strictly increasing, so
every interior point is in the strict good band. It computes (2) itself,
checks that all branches below each inclusive cap occur, and at leaf
time p/q checks n·min(wp mod q, q−(wp mod q))>q for each speed. Exact
endpoints need not be maximal for this proof, and no result from
the smaller-runner preprints is used by the checks.

| k | Deleted subsets | Recorded prefix counts, then leaf count | Largest caps at each depth |
|---|---:|---|---|
| 2 | 105 | 105, 1,188, 697 | 68, 50 |
| 3 | 455 | 455, 8,149, 24,965, 3,084 | 94, 81, 56 |
| 4 | 1,365 | 1,365, 38,237, 274,721, 303,607, 9,341 | 140, 176, 105, 60 |
| 5 | 3,003 | 3,003, 130,243, 1,843,194, 6,149,515, 2,881,258, 25,182 | 157, 189, 192, 132, 72 |

The caps are maxima over prefixes and need not be attained in full
tuples. The selected rational witnesses give exact gaps individually;
they are not a claimed global optimum or quantitative uniform margin.

## Replay without downloading certificates

Run `sh replay.sh` for k=2,3,4, or `sh replay.sh 5` for the
longer five-replacement check. Use GCC or Clang C++17 (including
128-bit integer support) and Python 3.10+. Everything works offline.
The script generates temporary text certificates, independently
checks them, then deletes them. No large data files or external code are bundled. Any timeout,
representation cap, invalid witness or unexpected truncation is not
a proof of an unfinished case.

## Scope and attribution

- [Yuhan Zhang, *Tight instances*, v3](https://arxiv.org/pdf/2608.13599v3)
  gives the closed-cover setup (Lemma 2.2), a related counting bound
  (Lemma 6.1), a complete two-replacement finite census through maximum
  added speed 8n (Theorem 6.10), and a deeper finite census for removals
  {2,3} through 40n for 8≤n≤18 (Theorem 6.4). His finite two-replacement bound 128
  at n=16 exceeds all terminal speeds retained in our k=2 proof.
  The difference here is the analytic unbounded tail rule and the
  unrestricted k=3,4,5 replacement patterns.
- [JenW1N’s October 2026 Lean proof and review](https://conjectures.io/results/4c651db1-3b53-48da-ac31-0b0020a8d87c)
  force matched multiplicative/gcd criteria in the upper-half
  two-removal case. The underlying proof already implies strictness
  for the 28 pairs from {8,…,15} at n=16: see the
  [full Lean source](https://conjectures.io/results/4c651db1-3b53-48da-ac31-0b0020a8d87c/solution/download),
  `target_of_remaining_window`, and the gcd criterion. In this use
  equality is used only to rule out n≤4, and the arguments work
  from ML≤1/n at n=16; the individual gcd conditions already
  fail at multiplier 2. This attribution is an inference from
  the proof, not an independent Lean rerun.
- [Allikvere, *Fourteen and fifteen lonely runners*, v2](https://arxiv.org/pdf/2609.02604v2)
  reports the largest total runner counts in the checked preprint
  literature. There is no claim here that the fifteen-total-runner
  result is peer reviewed, nor that the present certificate verifies it.

A separate general observation: if the conjecture is known through
n−1 total runners and r of k replacements remain, an installed prefix
of maximum speed M has a strict-good interval of length at least
2r/[n(n−r)M] by the lower-runner gap and the M-Lipschitz property.
Equation (2) gives next speed ≤(n−2)(n−r)M/(n−2r). For fixed k this
is O_k(n^(k+1)) after iteration from a core in {1,…,n−1}. The
weaker bound already in Zhang yields the same conditional asymptotic
order, with n+2 replacing n−2 in each recursion factor; this note
does not claim that order as novel or unconditional.

AI-assisted research note; not peer reviewed. No claim of priority.
