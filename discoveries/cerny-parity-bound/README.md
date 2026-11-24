# A parity-distance bound for the three outstanding binary quasi-Eulerian configurations

Status: independently checked research note, 24 November 2026. This
does not settle the Černý conjecture globally, or in all members of the
three configurations. It supplies an unconditional uniform bound with
leading coefficient 5/4 in these configurations and proves the conjectured
bound for a polynomially recognizable subfamily.

## Exact statement

Let Q be the n-state set (n≥4) of a strongly connected synchronizing
binary automaton, with letters a,b. Fix three distinct states x_1,x_2,s.
Write x=x_1, y=x_2. There are exactly two incoming arcs, counting both
labels with multiplicity, at every state outside {x,y,s}. Both letters
send each exceptional state x,y into {x,y,s}. Name E=Q\{x,y}. Exactly
two arcs from E have their target in {x,y,s} (by counting all arcs). The
following are the three cases, up to renamings:

|Case|(x.a,x.b,y.a,y.b)|Two boundary targets from E|w(x),w(y); ordinary w|
|---|---|---|---|
| A |(s,x,x,s)| y,s |1,1; 2|
| B |(s,s,x,x)| y,s |1,1; 2|
| C |(s,x,x,s)| x,y |3,1; 2|

These are exactly the exclusions in Theorem 5.9 of [Goel][goel]
(with the paper's X/E and boundary notation).

Initially assign binary labels 1 to x and 0 to y,s. A step on a flips
the assigned label, and a step on b preserves it in cases A,C. In case
B both steps flip it. Regard this as an ordinary reachability process
on (state,binary label), with the three labeled seeds. Let τ be the
minimal t such that some state has been reached at *both* labels along
paths of length at most t (paths may have different lengths). It is
a breadth-first-search quantity with 2n nodes. It is finite: otherwise
all reachable states, which comprise Q by strong connectivity, receive
one label each and the letters induce permutations of a nontrivial
two-state quotient. A reset word would be impossible. In our labeling
τ≥1.

Define D=min(n−1,τ+1). The following integer quantities bound the reset
threshold, with empty sums interpreted as 0. Let

                 E_A=E_B = ∑(j=1..n−2) max(j,n−1−j),
                 E_C = ∑(j=2..n−1) max(j−1,n−j−1).

Let L_A=L_B=n−2 and L_C=n−1. For h=1,2 put G(h)=L. For h≥3 put
```
G_AB(h) = max(0, min(n−2, n+1−2 ceil(h/2)))
G_C(h)  = max(0, min(n−1, n−2 floor(h/2))).
```

Then

                      rt(A) ≤ E + ∑(h=1..D) G(h).                 (1)

In particular, taking D=n−1 gives a uniform (5/4)n²+O(n) bound for
each of the three configurations. For D≥2 the odd sums can also be
written exactly as
```
O_AB(n,D) = nD − ceil(D²/2) − 2,
O_C(n,D)  = nD − floor(D²/2).
```
The exact test that the right-hand
side of (1) is at most (n−1)² proves the Černý bound for any qualifying
automaton, independently of its reset word. Asymptotically, for each
fixed η>0 it succeeds whenever τ≤(1−1/sqrt(2)−η)n and n is sufficiently
large. This slogan is not the finite test; use the exact inequality.
Examples of largest sufficient τ for A/B,C respectively:
n=13 (3,4), n=20 (5,6), n=30 (8,9), n=100 (29,29). These cutoffs
refer to (1), not to a necessary condition for the Černý bound.

## Proof

**Stationary extension costs.** Write W(K)=∑(q∈K)w(q), with w from the
table. At each ordinary target outside {x,y,s} the two incoming arcs
both come from ordinary states, hence their weights sum to 4. Directly
at x,y,s, the sums of incoming weights across both letters are 2w
as well: in A and B the sums are 2,2,4; in C they are 6,2,4.
Thus w is stationary for the average of the two transition matrices.

For each nonempty proper set K define m(K) as the least word length
whose preimage of K changes W(K). It is also the least word length
whose preimage raises W(K): averaging all words of any fixed length
returns W(K), so a decrease at that length implies there is an
increase at that same length. Existence with m(K)≤n−1 follows by the
standard span argument, supplied here to fix conventions. 
Consider T_i=span{w M_u−w: |u|≤i}. Then
T_{i+1}=T_1+∑(c=a,b)T_i M_c. If T_i=T_{i+1}, the sequence stabilizes;
all T_i lie in the (n−1)-dimensional space of rows of sum zero, since
each M_c is row stochastic. Thus T_{n−1} is the stable span.
If there had been no mass change through n−1, then 1_K would be
annihilated on this full span. A reset word followed by a path to a
state inside K changes W(K) to W(Q), a contradiction.

Start from a singleton of mass 2 (s) in A,B, and the singleton x of
mass 3 in C. Repeatedly take a preimage by a minimally long mass-raising
word, reaching all Q. The concatenation in reverse order resets into
the starting singleton. Each proper mass value is visited at most
once; separate their even and odd values.

**Parity observation.** For any word u, write c(q)=1_K(q). The parity
of W(Ku^{-1}) is c(x.u)+c(y.u) modulo two, because only the weights at
x and y are odd. Because x.a=s and y.a=x in all cases, the difference
modulo two between W(K(au)^{-1}) and W(Ku^{-1}) is c(s.u)+c(y.u).
Consequently if m(K)≥h, then for every |u|≤h−2,
```
c(s.u)=c(y.u),
c(x.u)+c(y.u) ≡ W(K) (mod 2).
```
The second relation holds through length h−1 as well. For h=1 the
first range is empty. These are consequences of *no change*, not merely
of having found no gain in a nonaveraged set of words.

**Even levels.** If W(K) is even and m(K)≥2, parity at u=empty forces
the whole triple {x,y,s} to lie on one side of K. Suppose the triple
is inside K and k=|K|. The first exit from K by a word starting at the
triple occurs from s (not x or y, whose first moves stay in the triple).
If its length is f, the joint reachable balls of lengths 0,...,f−1
strictly increase and remain in K by strong connectivity.
Starting with size 3 yields f≤k−2. For an exiting word u of length f,
x.u and y.u are still in K by first-step minimality, while s.u is
not. Then the parity observation shows W(K(au)^{-1}) and W(Ku^{-1})
have opposite parity. At least one differs from W(K), and
m(K)≤f+1≤k−1. For the triple outside K, apply the same argument to
the complement; m(K)=m(Q\K), so the bound becomes n−k−1.

In A,B an even proper visited mass 2j has 1≤j≤n−2. It is represented
by a K of size j if the triple is outside and j+1 if the triple is
inside. Its cost is bounded by max(1,j,n−j−1)=max(j,n−j−1).
In C the even proper levels have 2≤j≤n−1 due to the start at mass 3,
and k=j whether the triple is in or out. Their cost is bounded by
max(1,j−1,n−j−1)=max(j−1,n−j−1) for n≥4. This proves the E terms.

**Odd levels: identification with the prescribed parity cover.**
When W(K) is odd and m(K)≥h≥2, labels of x.u and y.u are opposite
for |u|≤h−1, and labels of y.u and s.u agree through |u|≤h−2.
Set d=c(s). For |u|≤h−2 this forces
```
c(s.u) = d XOR π(u)     c(y.u) = d XOR π(u)
c(x.u) = 1 XOR d XOR π(u).
```
Here π sums the flips of the letters mod 2. For s this follows by
inducting on the **first** letter of u: replace c(s.au') by c(y.au'),
then use y.a=x; in A,C replace c(s.bu') by c(y.bu')=c(s.u'),
and in B use y.b=x, another flip. The corresponding equations for
x,y then follow from the parity observation. This also proves m(K)≤τ+1:
if τ≤m(K)−2, the two conflicting paths at depth ≤τ would force
different indicators at one state. Thus every odd stage cost is at most D.

**Odd levels: balancing the visible colors.** Fix h≥3 and a set with
m(K)≥h and odd mass. Let R_i be the states reachable from any of x,y,s
within i steps. For i≤h−2 the preceding forced coloring partitions R_i
unambiguously into R_i^0,R_i^1 (assigning y,s color zero, x color one).
These sets are increasing in i. Every transition whose source is in
R_i for i≤h−3 obeys the designated flip and has its target in
R_{i+1}; an arc from the last layer is not assumed to obey it. By strong
connectivity |R_i|≥min(n,3+i), so |R_{h−3}|≥h. Let R^b=R_{h−2}^b.

In A each state has indegree 2 except y of indegree 1 and s of
indegree 3. Both are always in R^0, so the total indegree capacities
of each color are exactly twice their sizes. Each source in R_{h−3}
sends one arc into each color (a flips, b preserves), hence both
|R^b|≥ceil(h/2).

In B the indegree capacities are the same. Both arcs flip: each
source of color b in R_{h−3} uses two incoming slots of R^{1−b}.
Thus R^{1−b} has at least the size of the previous color b, and by
inclusion also at least the size of the previous color 1−b.
Both |R^b|≥ceil(h/2).

In C y has indegree 1, x indegree 3, and all other states indegree 2.
Thus the capacities of R^0,R^1 are respectively 2|R^0|−1 and
2|R^1|+1. Each source in R_{h−3} sends one arc to each color, giving
|R^0|≥ceil((h+1)/2) and |R^1|≥ceil((h−1)/2).

For A,B either color is included in K and the other is excluded.
It includes exactly one state of weight 1. So its odd mass and that
of its complement are each at least 2 ceil(h/2)−1. The full mass
is 2n−2, hence the number of possible odd masses with m≥h is at most
n+1−2 ceil(h/2), clipped below by zero and above by L. For C the
included color 0 gives weight ≥2 ceil((h+1)/2)−1, and included color
1 gives ≥2 ceil((h−1)/2)+1. These are equal. Complement has the same
lower bound, and full mass 2n, giving at most n−2 floor(h/2) odd
values. The total number of odd values in the starting-to-full
range is L for h=1,2. By summing the counts at each integral
height h, the sum of odd costs is at most ∑(h=1..D)G(h).
This proves (1).

Finally E=(3/4)n²+O(n), while ∑(h=1..D)G(h)=nD−D²/2+O(n) uniformly
for 2≤D≤n−1. At D≤n−1 the leading coefficient is at most 5/4.
For τ≤(1−1/sqrt(2)−η)n, the positive savings in the leading coefficient
beat all linear remainders, so the exact Černý bound holds eventually.

## Stronger bound using the actual reachable mass profile

The same argument retains full balls rather than just their minimum
sizes. For h≥2 let R=R_{h−2}, and c=n−|R|; all c states outside R
are ordinary and have weight 2. For even W(K) and m(K)≥h the
mod-2 equations make all three indicators equal for every word
of length ≤h−2. The same first-letter induction as above, now
with no flips of the indicators, shows that the entire ball R
is inside K or outside it. Thus its possible even masses lie in
`{2j, W(R)+2j: 0≤j≤c}`.
For odd W(K) and m(K)≥h the prescribed lift must be consistent on
R, and its two colored classes R^0,R^1 must be respectively fully
included or excluded (either orientation). Thus its possible odd
masses lie in `{W(R^0)+2j, W(R^1)+2j: 0≤j≤c}`;
if the lift is inconsistent there are none. In each parity keep
only levels between the starting mass (inclusive) and the full
mass (exclusive). At h=1 use all such levels (no tighter restriction
needed). Intersect each level set with the preceding sets, sum their
cardinalities over h=1,...,n−1. This is another upper bound on the
same extension trajectory, never worse than (1), computed in
polynomial time by `profile` in `bound.py`.

## Scope and possible next step

The argument is *not* a proof for late collisions. A future attempt
must improve the bound for large τ, not re-run small-state enumerations.
The precise potential obstruction to improving (1) by level counts
alone is that the estimate charges independently all aligned even
levels and the odd levels that rigidly agree with the same large
two-color reachable ball. A useful next test would prove that these
two expensive behaviors cannot alternate freely along the *same*
mass-raising trajectory, or exploit the nearly two-color quotient
present when τ is large. The first exit lemma and stationary
extension ideas are standard and appear in Goel. The new step to audit
is the prescribed parity cover and indegree-controlled count of slow
odd mass values. No priority or peer-review claim.

For comparison, besides the original degree-two estimate
4(n−1)², Berlinkov's general stationary-weight theorem directly
specializes to (2+o(1))n² for these weights, and the elementary
even-level accounting already gives (7/4)n²+O(n) when every odd
level is simply charged n−1. The parity improvement is from that
7/4 immediate baseline to 5/4; do not imply 4 is the sharp known
starting point for these exact configurations.

### A concrete unbounded limitation of both tests

For any odd n≥9 take the cyclic permutation in order 1,0,2,3,...,n−1
as a and let b fix all states except 1 maps to 2. Let u=n−5,v=n−3;
swap the images of u,v under a, and also under b. This is case A
with (x,y,s)=(0,1,2). Now a has a main odd cycle of size n−2 and
a 2-cycle containing v; b transposes u and v, bridging the cycles,
and identifies 1 with 2. The automaton is strongly connected.
The letter power b² identifies 1 with 2 and fixes every other state.
Acting with conjugates of b² by powers of a gives directed two-step
merges around the main cycle (1 and 2 are two steps apart), so
successively empty the main cycle into one state, leaving the small
cycle untouched. Rotate the main survivor away from u before each
b step; b then moves v to u without moving the survivor off the main
cycle. Rotate the small cycle if necessary so an unmerged small
state is at v. Repeating and then merging on the main cycle proves
synchronization.

The intended coloring 0 on state 1, 1 on state 0, and alternating
with color 0 at state 2 along the remaining natural order is obeyed
by all transitions except a from n−1 to 1. From the seed triple the
shortest route to this arc goes through s=2 and u=n−5, then jumps to
n−2 and advances twice (the last step is the violating arc).
Thus τ=n−4. Through depth i=n−7 the joint ball has i+3 states,
then adds two states in each of the next two steps. At levels
`h≤n−5` before these two steps the even and odd profiles have
exactly the universal numbers of possible mass levels: the ordinary
color counts differ by at most one, so the two odd ranges overlap
or abut on the odd lattice. At the tail h=n−4,n−3,n−2,n−1 the
universal and actual even counts are respectively (6,4,2,0) and
(4,0,0,0); the universal and actual odd counts are (4,4,2,2) and
(4,2,0,0). The universal AB total for odd n is
5(n−1)²/4−2. The actual profile loses 8 even and 6 odd levels,
so it still equals 5(n−1)²/4−16, above Černý for odd n≥11.
This is a limitation of the proof, **not a counterexample to the
conjecture**. For n=13 the exact reset threshold is only 56,
confirmed by two subset BFS methods in `late.py`.

[goel]: https://hexagonmath.org/pdf/2611.00635v1

## Reproduce the checks

The bound is proved above; finite tests are controls, not the proof.
A separately written C++ checker exhausts 15,696 prescribed labeled
completions at sizes 4–6 and tests actual extension costs and parity
layers. It also includes six witnesses
at sizes 13 and 17 outside the compared known structural classes; see
the [independent audit](audit.md) for the precise scope and overlap
checks. The Python scripts supply three further twenty-state examples,
a late-collision example and the bound arithmetic. No packages beyond
the Python standard library or C++17 standard library are required.

```sh
python3 bound.py
python3 examples.py
python3 late.py
python3 pilot.py 8 160
python3 pilot.py 10 240
c++ -O3 -std=c++17 -Wall -Wextra -Wpedantic replay.cpp -o /tmp/cerny-audit
/tmp/cerny-audit
```
