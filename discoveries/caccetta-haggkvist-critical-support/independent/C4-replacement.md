# Four-cycle support cannot extend (without a numerical triangle bound)

Suppose \(n=3r\), every outdegree is at least \(r\), every indegree at least \(r-1\), and the orientation has no directed triangle. Suppose a nonzero weighting \(p\) is supported precisely on four critical vertices \(S=\{s_0,s_1,s_2,s_3\}\). Assume \(p(N^-(v))=p(U(v))\) for every vertex \(v\), with \(U(v)\) the strict nonneighbors of \(v\). The directed induced cycle on \(S\) has edges \(s_i\to s_{i+1}\), modulo four. As the internal balance at \(s_i\) identifies its predecessor's weight with its opposite vertex's weight, the four weights coincide.

For any outsider \(v\), let \(I(v)\) be the support vertices sending to it, \(O(v)\) the support vertices receiving from it, and \(M(v)\) the missing support vertices. The balance gives \(|I(v)|=|M(v)|\). Thus the possible types are:

- \(A\): \(I=M=\varnothing\); these vertices send to all of \(S\).
- \(B_i\): \(I=\{s_i\}\), \(M=\{s_{i-1}\}\), \(O=\{s_{i+1},s_{i+2}\}\). Indeed \(s_{i-1}\) could not belong to \(O\) without closing a triangle.
- \(T\): \(|I|=|M|=2\), \(O=\varnothing\).

There is no edge into \(A\) from \(S\), nor from another type: if \(v\notin S\cup A\) sent to \(a\in A\), choose \(s_i\in I(v)\), obtaining \(a\to s_i\to v\to a\). But each \(s_i\) has indegree \(r-1\), already including an edge from its predecessor and all vertices of \(A\); hence \(|A|\le r-2\). If \(A\) were nonempty, each vertex in \(A\) would consequently have indegree at most \(|A|-1<r-1\), a contradiction. So \(A\) is empty.

Write \(b_i=|B_i|\). The total number of outsiders and the total number of arcs out of \(S\) to them give
\[
\sum_i b_i+|T|=3r-4,\qquad \sum_i b_i+2|T|=4(r-1).
\]
So \(|T|=r\). A vertex \(v\in T\) cannot send to \(S\) and has at most \(r-1\) possible targets within \(T\). As it has at least \(r\) outneighbors, it must send to some \(B_i\). If \(v\to b\in B_i\), then \(I(v)\) misses \(O(b)\), or else \(v\to b\to s\to v\) is a triangle. This forces \(I(v)=\{s_{i-1},s_i\}\), and that same \(B_i\) is the only type outside \(T\) to which \(v\) can send.

Let \(y_i\) be the number of \(v\in T\) with that consecutive input pair. The critical indegrees and outdegrees of \(S\) yield
\[
 b_{i-1}+b_{i-2}=r-2,\qquad b_i+y_i+y_{i+1}=r-1.
\]
The first equation gives \(b_i=b_{i+2}\); subtracting opposite instances of the second gives \(y_i=y_{i+2}\). Since \(\sum_i y_i=r\), it follows that \(y_i+y_{i+1}=r/2\) and consequently \(b_i=r/2-1\) for each \(i\). (An odd \(r\) would already be impossible.) Each \(v\in T\) therefore has at least \(r-b_i=r/2+1\) outneighbors in \(T\). An oriented graph on \(r\) vertices has average outdegree at most \((r-1)/2\), a contradiction. This excludes the four-vertex support using only triangle-freeness and the assumed degree conditions.
