# A rectangle-free bound next to a projective-plane order

Write $z(M,N;2,2)$ for the largest number of ones in an $M\times N$ binary matrix with no all-one $2\times2$. This note extends the incidence argument in [the Research1 order-ten bounds](README.md) to an unequal pair of sides.

**Theorem.** If $q\ge9$ and there is no projective plane of order $q$, then
$$
z(q^2+q,q^2+q+1;2,2)
\le (q+1)(q^2+q)-\lfloor q/2\rfloor-2.
$$
In particular, [Lam, Thiel and Swiercz's order-ten result](https://doi.org/10.4153/CJM-1989-049-4) gives $z(110,111;2,2)\le1203$. This is an upper bound, not an exact value.

## 1. The signed crossing estimate

For any linear incidence structure with $P$ points, $B$ lines and $E$ incidences, distinct points share at most one line and distinct lines meet at most once. Write $a_M=q+1-|M|$ for the signed line defect, $X=\sum_M a_M(a_M-1)/2$, and $L$ for the number of uncovered point pairs. The identity
$$
L+X=\binom P2-B\binom{q+1}{2}+q\big((q+1)B-E\big)
$$
follows by counting covered pairs. Here and below the polynomial $a(a-1)/2$ is nonnegative also when $a<0$.

Fix a line $A$ of size $s$, with $w$ other lines disjoint from it. Each other line that meets it covers exactly $|M|-1$ pairs with exactly one end on $A$. Thus the number left uncovered across $A$ is
$$
s(P+1-s)-(E-(B-1))+\sum_{M\ne A,\ M\cap A=\varnothing}(|M|-1).
$$
After adding the contribution from $A$ to $X$, each disjoint term combines with its summand in $X$ to give $q-1+(a_M-1)(a_M-2)/2\ge q-1$. Therefore
$$
L+X\ge g(s)+(q-1)w,\qquad
g(s)=s(P+1-s)-(E-(B-1))+(q+1-s)(q-s)/2. \tag{1}
$$
The function $g$ is a concave quadratic; specifically $g(P)-g(t)=(P-t)(P-2q+1-t)/2$.

## 2. Both degrees are at most $q+1$

Put $N=q^2+q$, $V=N+1$, $j=\lfloor q/2\rfloor+1$, and assume the rectangle has exactly $E=(q+1)N-j$ marks. Any larger example could be thinned to this size. For $P=V,B=N$, equation (1) gives $L+X=N/2+qj$, $g(q+2)=q^2+j$, and slack $K=(q-1)(j-q/2)\le q-1$. The two endpoints of $[q+3,V]$ satisfy $g(V)>g(q+3)>L+X$: the first follows from the endpoint formula and the second since $g(q+3)-g(q+2)=q^2-q-1>K$. Concavity excludes all sizes above $q+2$.

Suppose $A$ has size $q+2$. There is at most one other line disjoint from $A$, and $q^2-1$ points outside it.

- If $w=1$, equality holds in (1): the disjoint line has defect 1 or 2, every other line except $A$ has size at most $q+1$, and all outside pairs are covered. Each line through an outside point covers at most $q-1$ other outside points, whether it meets $A$ or is the one short disjoint line. Because $q(q-1)<q^2-2$, every outside point has degree at least $q+1$. The degrees on $A$ sum to $(q+2)+(N-2)$, so the signed point defect there is $q+2$. But the total signed point defect is $q+1+j>q+2$, incompatible with all outside defects being nonpositive.
- If $w=0$, let $h$ count the other $(q+2)$-lines and $L_{\rm out}$ count uncovered outside pairs. Each such oversized line contributes one to $X$; the exact crossing count gives $h+L_{\rm out}\le K\le q-1$. An outside point of degree at most $q-1$ could cover at most $(q-1)^2+h$ other outside points and miss at most $q-1-h$, too few to account for $q^2-2$ others. Thus every outside degree is at least $q$. The signed point defect on $A$ is now $q+1$ (all $N-1$ other lines meet it), leaving an outside signed sum $j$. A positive outside defect can only be one, so there are at least $j$ degree-$q$ outside points even before the point-side cap is known.

  Choose $j$ such points. For a chosen $x$, let $h_x$ count its oversized lines and $\delta_x$ its missing outside partners. Comparing $q(q-1)+h_x$ with $q^2-2-\delta_x$ gives $h_x+\delta_x\ge q-2$. Thus
  $$
  \binom h2\ge\sum_x\binom{h_x}{2}
  \ge\sum_x h_x-j
  \ge(q-3)j-2(q-1-h).
  $$
  The first inequality uses uniqueness of the meeting point of two lines; the last uses $\sum_x\delta_x\le2L_{\rm out}$. Subtract $\binom h2$ from the final bound. This concave expression is positive at both endpoints of $[0,q-1]$: using $j\ge(q+1)/2$ its endpoint values are at least $(q^2-6q+1)/2$ and $(q-5)/2$. Both are positive for $q\ge9$, a contradiction.

In the transpose $P=N,B=V$, one has $L+X=qj$ and $g(q+2)=q^2-q+j-1$. Its excess over $qj$ is $(q-1)(q-j)-1>0$, and the endpoint formula gives $g(N)>g(q+2)$. This caps all point degrees at $q+1$ as well.

## 3. Full lines force too many low points

Call a line full at size $q+1$, and a point low at degree at most $q$. There are at least $N-j$ full lines. If every two of them met, [Dow's theorem](https://doi.org/10.1016/0012-365X(83)90036-5), as restated in [Stinson, Theorem 1.1](https://doi.org/10.1017/S1446788700021728), would extend them to a projective plane of order $q$: $N-j>q^2-2\sqrt{q+3}+6$. Thus there is a full line $A$ disjoint from another full line. Call any full line with such a neighbor nonspanning.

Every point on a nonspanning full line $A$ lies on at most $q$ full lines. Otherwise its $q+1$ full lines would cover all other $q(q+1)$ points, yet only $q$ of those lines could meet the $q+1$ points of a full line disjoint from $A$. So $A$ meets at most $(q+1)(q-1)=q^2-1$ other full lines and has at least $t=q-j$ disjoint full neighbors.

For a degree-$(q+1)$ point $x$, the number of points not collinear with it is $\sum_{M\ni x}a_M$. On a nonspanning full line this is at least one, since the $q$ other lines through $x$ cannot cover a disjoint full line. The sum of these missing-partner counts over full-degree points on $A$ is at most $j$: each short line hits $A$ at most once. Thus each nonspanning full line has at least $q+1-j=t+1$ low points.

If a full-degree point on $A$ misses exactly one point, every disjoint full neighbor contains that same missed point. The neighbors intersect one another exactly there, hence together with $A$ use at least $(t+1)+1+t^2$ distinct low points (the lower bound lets the shared point be low). Since $t\ge q/2-1$, this is at least $q^2/4-q/2+2$, while the total point defect is at most $3q/2+2$. Their difference is strictly positive for $q\ge9$.

Therefore all full-degree points on every nonspanning full line miss at least two points. Each such line consequently has at least $q+1-\lfloor j/2\rfloor$ low points. Taking one nonspanning full line and two of its disjoint full neighbors gives at least $3(q+1-\lfloor j/2\rfloor)-1$ distinct lows. This exceeds the total point defect $q+1+j$, since the difference is at least $3q/4-3/2>0$. This final contradiction excludes $E$.

## Earlier bounds and scope

- [Damásdi, Héger and Szőnyi, Corollary 3.19, Eq. (3.3)](https://real.mtak.hu/13977/1/DHSzFIN.pdf) with $c=q+1$ gives $(q+1)(q^2+q)$; equality would embed in a plane of order $q$. Thus the earlier nonplane comparison is 1209 at $q=10$.
- A closer comparison follows from [Stinson, Corollary 2.9](https://doi.org/10.1017/S1446788700021728) and the crossing estimate. Set $j_0=\lfloor q/2\rfloor$. At this deficit, the slack at $q+2$ is negative for odd $q$ and zero for even $q$. In the latter case, a size-$(q+2)$ line would meet all other lines, be the only oversized line, and leave every outside pair covered. Then every outside degree is at least $q+1$ as above, although the signed outside defect is $j_0>0$. Thus at least $N-j_0\ge q^2+q/2$ lines are full, contradicting Stinson's strict packing bound. This elementary calculation gives 1204 at $q=10$; the theorem goes one beyond it. This is a derived comparator, not a claim that 1204 was printed in Stinson.
- The existing [110-square bound 1197](README.md) alone gives only 1207 by deleting a column of minimum degree, or averaging the 111 deletions. The new rectangle result is not a consequence of that square count. The order-ten exact maximum remains open. No exhaustive priority claim.
- The signed crossing and incidence proof is AI-assisted and independently checked, not peer reviewed. The optional replay uses only standard Python and confirms arithmetic and synthetic crossing controls; the proof does not rely on computation.

Run `python3 replay-rectangle.py` for the optional controls.
