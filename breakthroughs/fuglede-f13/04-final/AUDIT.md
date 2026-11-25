# Independent audit of linked-layer size 78

Private, 25 November 2026. The audited step is the last linked-layer obstruction; it takes the previously independently checked perpendicular linkage, cap five, and finite abscissa atlas as inputs. In conjunction with those results it eliminates the remaining layered size-78 case. The proposed integer profile kernel search was not run after the direct geometric finish emerged.

## Coordinate and completeness audit

In the established linked coordinates \(E=(x,f_s(x),s)\), \(A=(u,t,g_t(u))\), orthogonality implies \(\sum_t\zeta^{wt}C_t(r)=0\) for every normalized interlayer slope \(w\in I_r\). This follows by applying the same inverse interlayer-height Galois scaling to the original orthogonality relation for each generating pair, so the column coefficients are independent of the pair. If any column is active, Chebotarev bounds its entire union \(I_r\) by one less than the number of active columns, at most five. Every cross plane \(x-rs=x_0\) gives six selected points, one for each label; line cap five and the six-point direction gate force at least five slopes. Thus on each active \(r\), every such cross plane has exactly the same palette \(I_r\) of five slopes. The direct closure needs two active \(r\) only (the earlier argument gives at least ten), so stronger common graph-direction counts and chart-size bounds are not used here.

The six-point shadow gate sends the layer support \(S\) to a fixed affine normalization from the five 13-bit masks \(111,187,215,311,427\). This means one global coordinate change, not a relabeling per \(r\) or \(x_0\). It preserves the paired nature of the axes after the inverse change on the spectrum. Each canonical support contains its first two labels 0 and 1.

The independent audit source directly enumerates all residual height assignments after fixing \(h(0)=h(1)=0\), using actual attained slopes, not looping over prescribed direction palettes or selecting heights from \(sI\). With initial slope zero and four further heights, the search covers all \(13^4\) assignments on each support before exact pruning at six attained directions. There is one collinear zero residual, exactly twelve nonaffine residuals and no intermediate slope count. Undo the shear by \(h(s)\mapsto h(s)+as\) for every \(a\in\mathbb F_{13}\), leaving \(h(0)=0\). This recovers all 156 checked arrays on each support and determines their true palettes. Then compute all fifteen pair slopes directly to recheck the recovered palettes. The root and independent programs have the same conclusion; the independent count itself uses no root finite outputs as input other than the S list.

| Support mask | Palettes with models | Models per palette | Pairs invariant for *every fixed* palette |
|---:|---:|---:|---|
| 111 | 156 | 1 | All pairs |
| 187 | 78 | 2 | (1,7), (3,5) |
| 215 | 156 | 1 | All pairs |
| 311 | 156 | 1 | All pairs |
| 427 | 78 | 2 | (0,8), (1,7), (3,5) |

There is a useful simplification of the audit. Each support's twelve nonzero residuals form one multiplicative orbit, because a nonzero height vector has trivial scalar stabilizer. Thus all 156 zero-based profiles come from one representative \(h\) under \(a h(s)+b s\), \(a\ne0\). For 187 and 427 the representative palette is \(\{0,2,7,9,11\}\), with sole nontrivial affine symmetry \(w\mapsto9-w\). The two representative profiles give \(h(5)-h(3)=9\), unchanged by \(h\mapsto-h+9s\). Under any common palette transform this preserves equality, directly explaining the fixed pair. The other three palettes have trivial affine stabilizers. These observations are compatible with, but not needed for, the audit's exhaustive difference checks.

## Contradiction for a fixed pair

Fix a valid pair \(s\ne s'\) for the chosen support. For each active \(r\), every cross-plane height \(s\mapsto f_s(x_0+rs)\), normalized by the constant \(f_0(x_0)\), has palette \(I_r\). Its difference at \(s,s'\) therefore equals a value \(c_r\) independent of \(x_0\). Put \(z=x_0+rs\):
\[
 f_s(z)-f_{s'}\bigl(z+r(s'-s)\bigr)=c_r.
\]
Comparing two active \(r,r'\) produces a constant finite derivative of \(f_{s'}\) at nonzero displacement \((r-r')(s'-s)\). In the prime field this displacement generates the additive group, so the function is affine; its layer graph then has thirteen collinear points, violating cap five. The conclusion is unchanged if the derivative constant is zero. No direction or layer labels were selected depending on \(r\).

Replay the two adjacent files after restoring the .js and .expected suffixes. Default V8 output is the five exact JSON rows. No numerical Fourier computation or integer-kernel conjecture is used.
