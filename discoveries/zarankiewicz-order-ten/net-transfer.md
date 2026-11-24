# From dense Zarankiewicz matrices to finite nets

An addendum to [the Research1 order-ten note](README.md). Write $z(N)=z(N,N;2,2)$ for the maximum weight of a binary $N\times N$ matrix containing no all-one $2\times2$ submatrix. The previous note gives $z(110)\le1198$ and $z(111)\le1211$. Here we improve the first bound to
$$
z(110)\le1197
$$
and give an order-dependent upper bound at $N=q^2+q$ when no projective plane of order $q$ exists. The side-111 bound is unchanged.

## Net transfer

**Theorem.** Let $q\ge4$ and $1\le j\le q/2$ be integers, and set $N=q^2+q$. If an $N\times N$ rectangle-free matrix has at least
$$
E=(q+1)N-q-j
$$
ones, there is a strength-two orthogonal array on $q^2$ rows and $q+1-j$ columns with $q$ symbols and index one. Equivalently, there is a net of order $q$ and deficiency $j$, or a set of $q-1-j$ mutually orthogonal Latin squares of order $q$.

**Proof.** Delete ones if necessary so the weight is exactly $E$. View rows as lines and columns as points, with ones as incidences. Set $D=q+j$. For a line $M$ write $a_M=q+1-|M|$. Let $L$ count pairs of points not on a common line and $X=\sum_M a_M(a_M-1)/2$. Because pairs appear at most once,
$$
L+X=\binom N2-N\binom{q+1}{2}+qD
=\frac{q^2+(2j-1)q}{2}.
$$

For a line $B$ of size $s$, each other line intersecting $B$ covers exactly its size minus one pairs with exactly one point on $B$. The count of such pairs left uncovered is therefore
$$
s(N+1-s)-(E-(N-1))+
\sum_{M\ne B,\,M\cap B=\varnothing}(|M|-1).
$$
Adding the term for $B$ in $X$, and combining each disjoint line's sum term with its contribution to $X$, gives
$$
L+X\ge g(s),\qquad
g(s)=s(N+1-s)-(E-(N-1))+
\frac{(q+1-s)(q-s)}2.
$$
Indeed, each disjoint line contributes $q-1+(a_M-1)(a_M-2)/2\ge0$, even when $a_M$ is negative. The quadratic $g$ is concave. Its minimum on $[q+2,N]$ is attained at $q+2$, since
$$
g(N)-g(q+2)=\frac{(N-q-2)(N-3q-1)}2>0
\quad(q\ge4).
$$
But $g(q+2)=q^2+j-2>L+X$: twice the gap is $q^2-(2j-1)q+2j-4\ge q-2>0$. Thus all line sizes are at most $q+1$, and transposition gives the same bound on point degrees.

For a point $p$ set $h(p)=\sum_{M\ni p} a_M$. A point of degree $q+1$ has exactly $h(p)-1$ other points not on a common line with it. There are at least $N-D=q^2-j$ such full-degree points, and
$$
\sum_p h(p)=\sum_M |M|a_M\le qD
$$
since a line of positive defect has size at most $q$. As $2(q^2-j)>qD$ (the difference is $q^2-qj-2j\ge q^2/2-q>0$), some full-degree $p$ has $h(p)=1$. Its pencil consists of $q$ lines of size $q+1$ and one of size $q$. Removing $p$ from these lines partitions all other points into $q$ groups of size $q$ and one small group of size $q-1$.

The $q^2-1$ lines not through $p$ meet each group at most once. Record their incidences as rows of a table with one coordinate per group, leaving a blank when the line misses a group. Two rows cannot agree at two filled coordinates: the two corresponding lines would share two points. The table has precisely $D-1=q+j-1$ blanks: the pencil accounts for $(q+1)^2-1$ incidences, leaving $E-((q+1)^2-1)$ filled cells out of $(q^2-1)(q+1)$. Each small-group point is incident with at most $q$ lines outside the pencil, so at least $q-1$ blanks lie in the small-group column.

There are at least $q-j\ge2$ complete large-group columns. In each such column each symbol occurs at most $q$ times, by comparing with another complete column. Exactly one symbol occurs $q-1$ times, and all others occur $q$ times. In any pair of complete columns the sole missing ordered pair consists of the deficient symbol in each. Append one row carrying the deficient symbol in every complete column. All their two-column projections are now complete.

If at most $j-1$ blanks lie in the large groups, at least $q+1-j$ large columns are already complete, as needed. Otherwise there are exactly $j$ large-group blanks and the small group attains its incidence bound. Each of its $q-1$ symbols occurs in exactly $q$ old rows. Against any complete large column, a small-group symbol cannot occur twice with the same large-group symbol, so its $q$ occurrences pair once with each symbol in that column. After appending the new row, the $q$ rows whose small-group cells are blank therefore see each large-group symbol precisely once in every complete column. Fill all these $q$ small-group cells with one new symbol. The small-group column is now orthogonal to every completed large column. This again gives at least $q+1-j$ coordinates, and we retain exactly that many. The symbol classes in each coordinate form a parallel class, giving a net of degree $q+1-j$; alternatively, select two coordinates as the row and column indices and read the rest as Latin squares. $\square$

## Consequences

- **Side 110.** Set $q=10,j=2$. At 1198 ones the theorem gives nine orthogonal coordinates, hence seven MOLS of order ten. [Shrikhande's deficiency-two theorem](https://www.jstor.org/stable/25049140) ([open scan](https://www.renyi.hu/~mweiner/shrikhande.pdf), §2) would complete them to an affine plane. Adjoining the line at infinity gives a projective plane, impossible by [Lam, Thiel and Swiercz](https://doi.org/10.4153/CJM-1989-049-4). Thus $z(110)\le1197$.
- **Nonplane orders.** [Metsch's 1991 net-completion theorem](https://doi.org/10.1007/BF00157614) extends an order-$q$, deficiency-$j$ net to an affine plane if $3q>8j^3-18j^2+8j+4$. Whenever no projective plane of order $q$ exists, this inequality and the transfer hypotheses yield
  $$
  z(q^2+q)\le (q+1)(q^2+q)-q-j-1.
  $$
  One admissible choice is $j=\lfloor(3q/8)^{1/3}\rfloor$ for every $q\ge4$: it lies between 1 and $q/2$, and the correction to $8j^3$ is strictly negative. The full polynomial condition can allow larger $j$. For example, $q=30,j=3$ satisfies it and gives $z(930)\le28796$. Order 30 is ruled out by [Bruck–Ryser, Theorem 1](https://doi.org/10.4153/CJM-1949-009-2), as is every positive $q\equiv6\pmod{36}$: these orders are 2 modulo 4 and their square-free parts contain the prime 3.
- **Comparison.** [Damásdi, Héger and Szőnyi, Corollary 3.19(3.3)](https://real.mtak.hu/13977/1/DHSzFIN.pdf) with $c=q$ gives $(q+1)(q^2+q)-q$, with equality requiring an order-$q$ projective plane. On nonplane orders this gives the ceiling one below that value. The new ceiling is $j$ still lower, hence an additive improvement of order $q^{1/3}=\Theta(N^{1/6})$ over that specific comparator along infinitely many proven nonplane orders. This does not change the leading asymptotics or claim an improvement over all other bounds at each order. In particular, the [existing side-42 note](../zarankiewicz-42/) is stronger at $q=6$. Exact values remain open here.

The pencil method extends the earlier Research1 order-ten argument, itself motivated by [Sadhu's leave and packing approach](https://arxiv.org/html/2608.01606v2). The extra step is to complete the small coordinate after restoring the missing row, not to repair the punctured large coordinates. The completion theorems and plane nonexistence results are classical cited ingredients. This is an AI-assisted, independently rechecked research note, not peer reviewed; no exhaustive priority claim is made. Run `python3 replay-net-transfer.py` for exact arithmetic and constructive finite-field controls. The universal argument above does not depend on computation.
