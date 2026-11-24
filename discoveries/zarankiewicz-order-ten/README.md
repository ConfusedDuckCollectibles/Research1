# Two rectangle-free grid bounds near order ten

How densely can you mark a square grid without marking all four corners of a rectangle? Two grid sizes, 110 and 111, are especially close to what would be possible if a finite projective plane of order ten existed. It does not. The proofs below show that several near-misses are impossible, too.

Write $z(N)=z(N,N;2,2)$, the most ones in an $N\times N$ binary matrix with no all-one $2\times2$. We prove
$$
 z(111)\le1211,\qquad z(110)\le1198 .
$$
These are upper bounds, not determinations of either extremal number. The proof uses established results about finite planes and does not use any assumption about three or more unknown mutually orthogonal Latin squares of order ten. No new computer search is required. The strongest comparison bounds identified here are 1214 and 1199 respectively; the exact answers remain open.

## A crossing identity

Regard the columns as points and rows as lines. Any two points lie on at most one common line, equivalently any two lines meet at most once. Put $E$ for the number of incidences, $D=11N-E$, $a_M=11-|M|$, and $X=\sum_M a_M(a_M-1)/2$. All summands in $X$ are nonnegative, including those from negative defects. If $L$ counts uncovered unordered point pairs, counting covered pairs gives
$$
 L+X=\binom N2-55N+10D
$$
because $\binom{11-a}{2}=55-10a+a(a-1)/2$.

For a line $B$ of size $s\ge1$, let $w$ be the number of other lines disjoint from $B$. Each other line that meets $B$ covers exactly its size minus one cross pairs (pairs with exactly one point of $B$). Thus the number of uncovered cross pairs is
$$
 L_{\mathrm{cross}}=s(N+1-s)-(E-(N-1))+10w-\sum_{M\ne B,\ M\cap B=\varnothing}a_M.
$$
Using the term for $B$ in $X$ and $a(a-1)/2\ge a-1$ on each disjoint line yields
$$
 L+X\ge g(s)+9w,\qquad
 g(s)=s(N+1-s)-(E-(N-1))+(11-s)(10-s)/2 . \tag{1}
$$
The function $g$ is concave. For $N=111,E=1212$, $L+X=90$, and $\min\{g(12),g(111)\}=\min\{99,4059\}=99$, so no line has size at least 12. For $N=110,E=1199$, $L+X=55$, and $\min\{g(12),g(110)\}=\min\{99,3970\}=99$, giving the same cap. Transposition caps point degrees at 11 as well.

## Excluding 1212 incidences on side 111

Assume exactly $E=1212$, hence total point defect and total line defect both equal 9. A line or point is **full** when its size or degree is 11; otherwise it is **short** or **low**, respectively. Let $\alpha$ count short lines. There are at least $111-9=102$ full lines. [Dow's embedding theorem](https://doi.org/10.1016/0012-365X(83)90036-5), also stated by [Stinson, Theorem 1.1](https://doi.org/10.1017/S1446788700021728), says a family of 11-point lines on 111 points that are pairwise intersecting extends to an order-ten projective plane when it contains more than $106-2\sqrt{13}$ lines. Such a plane does not exist by [Lam, Thiel and Swiercz](https://doi.org/10.4153/CJM-1989-049-4). So a disjoint pair of full lines must exist. This conclusion also holds in the transpose.

For a full line $A$ in such a pair, let $k_A$ count full lines disjoint from it and $l_A$ count its low points. Two simple bounds are useful:

* $k_A\ge 11-\alpha$, as in Stinson's Lemma 2.1. Within the full-line packing, a point of $A$ cannot be on eleven full lines, because they would cover all other 110 points. In particular it would then need eleven *other* lines through it to meet all points of a full line disjoint from $A$. Thus each point of $A$ is on at most nine full lines other than $A$. Of the $110-\alpha$ full lines other than $A$, at most $11\cdot9=99$ meet $A$.
* $l_A\ge 11-\alpha$. Any full point $p\in A$ misses a point of each full line disjoint from $A$: it has only ten lines besides $A$ with which to reach eleven points on that line. Its number $\delta(p)$ of uncovered point pairs through $p$ is $\sum_{M\ni p}a_M$ when $p$ has degree 11. So it lies on a short line. Each short line hits $A$ at most once; at most $\alpha$ full points of $A$ can have such an incidence. Also, for any $A$,
  $\sum_{p\in A,\ d(p)=11}\delta(p)\le9$, by counting short-line weights.

If some full point $p\in A$ has $\delta(p)=1$, all full lines disjoint from $A$ pass through the unique point $x$ not covered with $p$. If $k_A\ge2$, that point is uniquely determined by the neighbors, so all full points with $\delta=1$ on $A$ miss this same $x$.

**Both defect profiles are unit.** Take disjoint full lines $A,B$. They have at most nine low points between them, so at least thirteen full points, with total $\delta$ at most eighteen. As all these full points have $\delta\ge1$, some has $\delta=1$; call its line $A$ after swapping the names if needed. If $\alpha\le8$, put $r=11-\alpha\ge3$. There are at least $r$ full lines disjoint from $A$, all passing through one point and otherwise pairwise disjoint. Each has at least $r$ low points. They contain at least $1+r(r-1)$ distinct low points outside $A$ (and at least $r^2$ if the common point is not low). Together with the at least $r$ low points on $A$, this demands at least $r^2+1\ge10$ low points. Contradiction. Thus $\alpha=9$, and transposition gives exactly nine low points as well. All positive defects are one.

**The disjoint full lines are impossible.** In this balanced profile every nonspanning full line $A$ (one with at least one disjoint full line) has $k_A,l_A\ge2$, and is disjoint from exactly $l_A$ lines in all: each of its points of degree $d$ accounts for $d-1$ other lines meeting $A$, without overlap. Also $l_A\le6$. Otherwise a disjoint full neighbor $B$ has $l_B\ge2$, so $l_A=7,l_B=2$ and all low points lie on $A\cup B$. A second full line disjoint from $B$ exists, and meets $A$ at most once. It would have at most one low point, a contradiction.

A nonspanning full line $A$ has at least five full points with $\delta\ge1$ and total $\delta\le9$. Some has $\delta=1$; let $x_A$ be the common point of its disjoint full neighbors. This hub $x_A$ is low. If it were full, the full point detecting it would be disjoint in the dual full packing from another full point, and so would have at least two uncovered full points by the transposed $k_A\ge11-\alpha$ bound, not one. Each of the $k_A$ disjoint full neighbors contains $x_A$ and at least one other low point, distinct across these neighbors, so
$$
 l_A+k_A\le8. \tag{2}
$$
Exactly $l_A-k_A$ short lines miss $A$, hence the full points on $A$ have total $\delta\le9-(l_A-k_A)$. Of its $11-l_A$ full points, the number with $\delta=1$ is therefore at least
$$
 2(11-l_A)-[9-(l_A-k_A)]=13-l_A-k_A\ge5.
$$
All these points miss $x_A$. The hub has degree ten; each of its lines meeting $A$ reaches one different point of $A$, and at most $l_A$ of its lines are disjoint from $A$. Thus $x_A$ misses at most $11-(10-l_A)=1+l_A$ points of $A$. Consequently $l_A\ge4$ for every nonspanning full line. Its $k_A\ge2$ disjoint full neighbors now contain at least $1+3k_A\ge7$ low points outside $A$, whereas at most $9-4=5$ remain. Contradiction.

Hence 1212 incidences are impossible. Deleting incidences from any larger example would yield exactly 1212, proving $z(111)\le1211$.

## Excluding 1199 incidences on side 110

Assume exactly $E=1199$, so $D=11$. The cap gives at least $110-11=99$ points of degree 11. For such a point $p$, the uncovered-pair count is
$$
 \delta(p)=109-\sum_{M\ni p}(|M|-1)=\sum_{M\ni p}a_M-1 .
$$
The total weighted incidence with short lines is $\sum_M |M|a_M\le10D=110$. Since $2\cdot99>110$, some degree-11 point has $\delta=0$. Its lines consist of one of size 10 and ten of size 11. Their points other than $p$ partition the other 109 points into groups of sizes 9 and $10^{10}$.

There are 99 lines not through $p$, meeting each group at most once. The lines through $p$ contribute 120 incidences, so the others contribute $1199-120=1079$. Thus there are $99\cdot11-1079=10$ missed group-line pairs. The group of size nine is met at most $9\cdot10=90$ times, since each of its points already lies on one line through $p$ and has degree at most 11. Consequently at least nine of the ten size-ten groups are met by **every** off-point line.

Using just those nine groups yields 99 words on ten symbols, length nine, with no two agreeing twice. In every two-coordinate projection, exactly one pair of symbols is absent. In each coordinate all symbols occur ten times except one symbol occurring nine times: no symbol can occur more than ten times (compare with any other coordinate). Thus the absent pair in every projection is exactly the pair of exceptional symbols. Appending their nine-coordinate word gives a complete strength-two orthogonal array of 100 words. Taking two coordinates as row and column indices, the other seven are seven mutually orthogonal Latin squares of order ten, which would extend to a complete set by [Shrikhande’s theorem, §2](https://www.renyi.hu/~mweiner/shrikhande.pdf). A complete set gives a projective plane of order ten, contrary to Lam, Thiel and Swiercz. The question of **three** such squares is not needed.

This excludes 1199 incidences, and deletion gives $z(110)\le1198$.

## Scope and sources

For 110, [Damásdi, Héger and Szőnyi, Corollary 3.19(3.3)](https://real.mtak.hu/13977/1/DHSzFIN.pdf) gives 1200 and requires a projective embedding in its equality case (order $n=10$, $c=10$), hence 1199 combined with Lam et al. For 111, one earlier ceiling is 1214 directly from the cap (1) at 1215 and the classical [Stinson uniform packing bound, Corollary 2.9 and Lemma 3.1](https://doi.org/10.1017/S1446788700021728) plus Dow and Lam et al.: ≥105 full lines must either intersect pairwise and extend or contradict the strict pseudo-partial upper bound of fewer than 105 lines. We make no exhaustive priority claim.

Run `python3 check-arithmetic.py` to check the numerical inputs to the proofs and the comparison bounds. This is an arithmetic check, not a proof checker for the cited theorems.

The crossing-and-leave approach is motivated by [Sadhu's order-six analysis and order-ten question](https://arxiv.org/html/2608.01606v2), especially Sections 3 and 13. The order-ten bounds are independent of any order-six packing computation. This work is AI-assisted and not peer reviewed.
