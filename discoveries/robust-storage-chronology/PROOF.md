# Exact reserve and infeasibility certificates from coarse storage chronologies

## Frozen setup and result

This is an abstract single-bus, single-store energy balance, not a model of the reliability or operation of a specific electric grid. A bounded measurable net power profile \(p(t)\) repeats with period \(T\). Positive means surplus generation after serving load. Only its chronological bin totals \(A_j=\int_{a_j}^{b_j}p(t)\,dt\) and pointwise envelope \(-P_-\le p(t)\le P_+\) are known, with \(P_\pm>0\). Consecutive bins partition \([0,T]\) and \(\Delta_j=b_j-a_j>0\). Totals are feasible in the envelope and the same completed sub-bin profile repeats every period. There is no ramp, statistical correlation, network flow, alternate generation, or loss-of-load allowance.

The store has nonnegative energy state \(X\le E\), standing-loss rate \(\rho\ge0\), charge/discharge efficiencies \(0<\eta_c,\eta_d\le1\), optional charging power cap \(P_c>0\), and enough discharge power to cover \(P_-\). Define the maximum available internal net rate \(f\) and rejected internal charge \(z\) by
\[
 f(p)=\eta_c\min(p_+,P_c)-p_-/\eta_d,\qquad
 \dot X=f(p)-\rho X-z,\qquad 0\le z\le f(p)_+,
\]
where \(p_\pm=\max(\pm p,0)\) in \(f\). The store may reject surplus but may not reject load; discharge is only to cover the realized deficit. Let \(c=f(P_+)>0\), \(d=-f(-P_-)>0\). Every admissible completed trace must be supportable indefinitely with this one proposed energy capacity \(E\); the policy may use instantaneous realized surplus/deficit. Periodic operation may choose its own initial state; starting full, the same greedy charging/rejection policy is also sufficient for all traces.

**The exact robust result.** For an unwrapped arc \([t,u]\), \(0\le t\le T,\ t\le u\le t+T\), let
\[
 {\cal M}(t,u)=\min_{p\ \mathrm{admissible}}
       \int_t^u e^{\rho(s-t)} f(p(s))\,ds . \tag{1}
\]
This minimum has an explicit envelope-only expression. For each bin \(j\) set
\[
 h_j=\frac{A_j+P_-\Delta_j}{P_++P_-},\qquad \ell_j=\Delta_j-h_j.
\]
In unwrapped order within the arc, let \(S_j\) be the *latest* \(\min\{\ell_j,|[t,u]\cap j|\}\) time units of bin \(j\) encountered by the arc. If the same physical bin appears in two pieces due to wrapping, treat the pieces together and fill the later piece first. Then
\[
 {\cal M}(t,u)=c\int_t^u e^{\rho(s-t)}ds-(c+d)\sum_j\int_{S_j} e^{\rho(s-t)}ds . \tag{2}
\]
At \(\rho=0\), “latest” can be any equal-length in-arc subset. This is an attained minimum, not a relaxation substituting independently chosen signs for the same phase.

* If any full-cycle cut has \({\cal M}(t,t+T)<0\), **no finite storage capacity** robustly supports the coarse records; some compatible repeat trace is physically impossible indefinitely even from an arbitrary finite initial charge.
* Otherwise the least robust energy capacity is exactly
\[
 E_{\rm robust}=\max_{0\le t\le T,\ t\le u\le t+T}[-{\cal M}(t,u)]. \tag{3}
\]
Both extrema reduce to an explicit finite set of at most \(O(J^2)\) candidate cuts for \(J\) bins, including cut points *inside* observed bins. There is no fine resolution parameter or numerical optimization over chronologies. The maximizing deficit cut has a matching bang-bang witness and all histories attain the claimed feasible guarantee via the upper-cap greedy policy (initialized full).

**Known baselines.** The lossless interval-deficit/sequent-peak method and the equivalence of deterministic single-reservoir storage sizing with linear programming are classical [Gidley 1986](https://doi.org/10.1016/0022-1694(86)90146-0); [Lele 1987, improved reservoir sequent peak](https://doi.org/10.1029/WR023i010p01819) and other reservoir literature already account for storage-dependent losses. [Storage requirements to mitigate intermittent renewable energy sources, 2023](https://doi.org/10.3389/fenvs.2023.1076830) use loss-aware modified cumulative deficits and a forward fill-and-spill check on known power histories; [self-discharge parameterization, 2026](https://doi.org/10.1186/s42162-026-00676-z) optimizes profiles for individual known load events using a discretized proportional loss and without cyclicity. [Merrick et al., chronological representation, 2024](https://www.sciencedirect.com/science/article/pii/S0140988324003839), [representative-period methods, 2026](https://doi.org/10.1016/j.apenergy.2026.128989), and [event-enriched aggregation, Nov. 2026](https://doi.org/10.3390/en19225423) document the consequences of losing chronological information using full-resolution traces or holdouts. [Exact temporal aggregation under perfect information, arXiv:2510.14451](https://arxiv.org/html/2510.14451v1#S3) uses active sets and duals obtained from the full trace, not a tight guarantee using only bin totals. [Zou–Lin–Aliprantis, robust multi-stage grid operation](https://staff.ie.cuhk.edu.hk/~xjlin/paper/infocom18-robust-storage.pdf) develops a different robust contract with generators, a network and intertemporal ramp bounds. The prospective contribution here is the exact robust threshold and constructive finite certificate for the explicitly impoverished coarse-bin/loss contract, not recognition that chronology matters, invention of cumulative deficits, or a forecast for a grid.

## 1. Exact deterministic capacity after the chronology is specified

Let \(f(t)=f(p(t))\) be any one bounded periodic completion. If \(\rho>0\), the unique periodic no-rejection solution on the unconstrained real line is
\[
 Y(t)=\frac{\int_{t-T}^{t} e^{-\rho(t-s)} f(s)\,ds}{1-e^{-\rho T}} .
\]
A feasible periodic solution cannot exceed \(Y\), by integrating the nonnegative rejection over one period. More generally, starting at the same phase in any finite initial nonnegative state, the state after \(n\) periods is bounded above by \(e^{-\rho nT}X(t)+(1-e^{-\rho nT})Y(t)\). If \(Y\) is negative at any phase, no indefinitely feasible finite initial state exists. The condition \(Y\ge0\) is exactly \(\int_t^{t+T}e^{\rho(s-t)}f(s)\,ds\ge0\) for every \(t\), by periodicity and positive rescaling.

Assume \(Y\ge0\). Any interval \([t,u]\) needs starting state at least
\[
 D(t,u)=-\int_t^u e^{\rho(s-t)} f(s)\,ds =Y(t)-e^{\rho(u-t)}Y(u),
\]
if this number is positive; the integrating factor gives the bound even with any rejection. This also shows that intervals longer than \(T\) cannot raise the maximum: remove any full periods from their end while holding the initial phase, and \(Y(u)\) is unchanged and nonnegative, while its multiplying factor only decreases. To prove sufficiency at capacity \(E\ge\max D(t,u)\), use the upper-barrier trajectory
\[
 X_E(t)=Y(t)-\sup_{s\le t} e^{-\rho(t-s)}[Y(s)-E]_+
       =\min\!\left(Y(t),\inf_{s\le t}\left[
          e^{-\rho(t-s)}E+\int_s^t e^{-\rho(t-v)}f(v)\,dv\right]\right).
\]
All terms are nonnegative by the interval condition and \(Y\ge0\), and taking \(s=t\) gives \(X_E\le E\). The function is periodic. For the differential law, multiply the running maximum by \(e^{\rho t}\); it becomes the usual nondecreasing running maximum of the absolutely continuous \(e^{\rho s}[Y(s)-E]_+\), hence is locally absolutely continuous and increases only at a freshly touching maximum with \(X_E=E\). Consequently \(X_E'=f-\rho X_E-z\), with \(z\ge0\) supported at the upper cap. Since the derivative of an absolutely continuous function is zero a.e. on a level set, at those times \(z=f-\rho E\le f_+\), so no forbidden discharge or load rejection occurs. This is the greedy cap policy. Starting at \(E\) instead of its periodic initial phase yields a trajectory no smaller than \(X_E\) by comparison with the scalar cap dynamics (or the same running maximum beginning at the initial time). It is feasible without knowing future \(p\).

If \(\rho=0\), the unbounded periodic prerequisite is simply \(\int_0^T f\ge0\); negative mean defeats any finite initial charge. In the feasible case let \(Y_0(t)=\int_0^t f\) on all real times. The periodic reflected state \(X_E=E+Y_0(t)-\sup_{s\le t}Y_0(s)\) lies in \([0,E]\) exactly under the interval deficit condition. The supremum is finite and uses at most the previous period because full periods have nonnegative net. It has the same feasible cap reflection law and comparison from full. This recovers the ordinary sequent-peak criterion.

## 2. Exact least favorable chronology for a given cut

The map \(f(p)\) is concave on \([-P_-,P_+]\), since its ordered slopes from left to right are \(1/\eta_d,\eta_c,0\) (the last slope only if charge clipping is active in the envelope). Hence \(f(p)\) is at least its secant between the endpoint rates \(c,-d\). Writing \(\theta=(P_+-p)/(P_++P_-)\) gives \(f(p)\ge c-(c+d)\theta\), where \(0\le\theta\le1\), \(\int_j\theta=\ell_j\). Equality holds for either endpoint raw power. The arc's weight is nonnegative on the physical phase bin (zero outside the arc). Maximizing the weighted integral of \(\theta\) at given mass \(\ell_j\) fills its highest weights first: replace any lower-weight occupied mass by the same higher-weight unused mass, or subtract the cutoff weight times the unchanged total mass. Since the weight is monotone nondecreasing in unwrapped time on included portions, this selects \(S_j\) exactly; fill any remaining deficit duration in the bin's excluded region of zero weight. Taking \(p=-P_-\) there and \(p=P_+\) elsewhere realizes the optimum simultaneously for all bins and respects each exact total. In particular wrapping does not allow an independent choice of the same bin's deficit budget at both ends. This proves (2). A different cut can have a different worst chronology; the feasibility statement quantifies over every chronology and cut, so order of taking the two minima is legitimate.

A negative \({\cal M}\) on a full-cycle cut constructs an incompatible-with-indefinite-operation trace by the deterministic theorem. If all are nonnegative, every possible trace has \(Y\ge0\) (or nonnegative total when \(\rho=0\)), and the least capacity is the worst interval bound (3) by the same theorem. On periodic initialization this is a robust sizing statement, not stochastic reliability. The full-start feedback result comes from one greedy policy, not perfect-foresight control.

## 3. Finite candidate cuts without temporal discretization

Use the unwrapped time domain \(0\le t\le T\), \(0\le u-t\le T\). Form an arrangement of the following straight lines and domain edges:

* \(t\) or \(u\) at bin boundaries, using the copy of the boundaries shifted by \(T\) for \(u\);
* \(t=b_j-\ell_j\) and \(u=a_j+\ell_j\) (also the \(T\)-shifted copies for \(u\));
* \(u-t=\ell_j\) or \(u-t=T-\Delta_j+\ell_j\).

Discard portions outside the domain. Include every vertex of the induced arrangement, including on \(u=t+T\); using redundant lines causes no problem. This takes \(O(J^2)\) cuts.

To see why, on each cell the start/end bins and changes of every selected subset have fixed membership. For a single included piece \([A,B]\) of a phase bin, the selected interval is \([\max(A,B-\ell_j),B]\); switch lines are exactly those where \(B-A=\ell_j\). For a split wrapping bin, its later piece is \([a_j+T,u]\), its earlier piece \([t,b_j]\). Use the later piece first. Regime changes are \(u-a_j-T=\ell_j\) and \((u-t)+\Delta_j-T=\ell_j\). If selection reaches but does not fill the early piece, its lower endpoint is \(b_j-\ell_j+u-a_j-T=u+\Delta_j-\ell_j-T\). All integration endpoints within a cell thus have the form constant, \(t\) plus constant, or \(u\) plus constant. When \(\rho>0\), (2) on the cell is affine in \(P=e^{-\rho t}\) and \(Q=e^{\rho(u-t)}\). Each original cell is an intersection of halfplanes imposing bounds on \(t,u,u-t\); in \((P,Q)\) these become linear halfplanes (\(t=\text{const}\) fixes \(P\), \(u-t=\text{const}\) fixes \(Q\), \(u=\text{const}\) gives \(Q=e^{\rho u}P\)), so the transformed domain is again a compact polygon. The affine objective takes its extrema at vertices. When \(\rho=0\), the cell objective is affine in the original coordinates and the conclusion is immediate. The same proof restricted to \(u=t+T\) proves the finite full-cycle test. Formulas extend continuously to cell boundaries and zero-length pieces. No numerical quadrature is needed for the exponential integrals.

## 4. Sharp remedy using separate surplus and deficit totals

Suppose the meter preserves \(G_j^+=\int_j p_+\) and \(G_j^-=\int_j p_-\) separately. Let \(h_j^+=G_j^+/P_+\) and \(h_j^-=G_j^-/P_-\). Compatible metering implies \(h_j^++h_j^-\le\Delta_j\), and the remainder may be zero raw net power. Let \(L_j(q;t,u)\) be the weighted integral on the latest \(\min(q,|[t,u]\cap j|)\) in-arc duration, again with wrap handled in unwrapped chronology. Then the exact replacement of (2) is
\[
 {\cal M}_{\mathrm{gross}}(t,u)=c\int_t^u e^{\rho(s-t)}ds-
     \sum_j \big[c L_j(\Delta_j-h_j^+;t,u)+d L_j(h_j^-;t,u)\big]. \tag{4}
\]
On the positive side \(f(p)\ge c p_+/P_+\) by the positive-side secant even with clipping; the negative contribution is exactly \(-d p_-/P_-\). The two extremal occupancies minimize the positive term by filling the smallest weights and maximize the subtracted negative term by filling the largest. Since their durations sum to at most \(\Delta_j\), these selections can be disjoint even after resolving ties of zero weight, and are realized by \(p\in\{P_+,0,-P_-\}\). The earliest positive duration is the complement of the latest \(\Delta_j-h_j^+\) across the whole bin, giving (4) after excluding zero-weight times. The deterministic capacity/feasibility tests remain exact, and the same finite-cut list works with both duration parameters in place of \(\ell_j\).

An exact example shows that recording gross totals can change the conclusion from infinite robust capacity to finite. Take one bin of duration 8, raw bounds \(\pm1\), ideal conversion, \(\rho=\log(6/5)\), and net total 2. Net-only metering allows duration 5 at \(+1\) and 3 at \(-1\), for which \(2(6/5)^5-(6/5)^8-1<0\). But if the gross totals are separately \(G^+=3,G^-=1\), the worst full-cycle schedule has three surplus units, four idle units, and one deficit unit. Its full-cycle weighted integral multiplied by \(\rho\) is \((6/5)^3-1-(6/5)^8+(6/5)^7>0\). The exact robust capacity now has \(\rho E=1296/3125\), attained by starting at the four idle units before the deficit. The offline replay independently enumerates every integer unit-slot trinary chronology for these data; the polygon proof in §3 makes this integer-vertex check valid for the continuous completion set.

## 5. Closed one-bin threshold and additional controls

With **one** bin covering the complete period let \(h=h_1\), \(\ell=T-h\). The worst full-period chronology charges for duration \(h\) before a contiguous maximum deficit of duration \(\ell\); for any full cut the entire phase bin is freely reorderable. For \(\rho>0\),
\[
 {\cal M}(t,t+T)=\{(c+d)e^{\rho h}-c-d e^{\rho T}\}/\rho,\qquad
 h_{\rm crit}=\rho^{-1}\log\frac{c+d e^{\rho T}}{c+d}.
\]
If \(h<h_{\rm crit}\) the robust capacity is infinite. At and above the threshold, the largest weighted deficit is exactly \(d(e^{\rho\ell}-1)/\rho\). Indeed an arc of duration \(q\le\ell\) is filled by deficit and its required reserve increases in \(q\). For \(q>\ell\), the worst integral gives the required start reserve \(\{c+e^{\rho(q-\ell)}[d(e^{\rho\ell}-1)-c]\}/\rho\). Its bracket coefficient is negative in the nontrivial feasible case (\(\ell>0\), compare \(q=T\) where reserve is nonpositive), so the maximum is at \(q=\ell\). At \(\rho=0\) take the continuous limit: \(h_{\rm crit}=dT/(c+d)\), finite reserve \(d\ell\). This is an engineering-model threshold rather than a probability estimate.

With symmetric internal rates \(c=d\), \(h_{\rm crit}=T/2+\rho^{-1}\log\cosh(\rho T/2)\). Thus even arbitrarily small positive standing losses create a strictly positive net-surplus requirement in the least favorable coarse-compatible profile. For example with symmetric raw bounds \(\pm1\), ideal conversion, \(T=200\), \(h=101\), and \(\rho=\log(1001/1000)\) per normalized time unit, the mean-positive constant-profile replacement needs **zero** storage, yet an admissible repeat pattern cannot be supported by *any* finite store. The exact sign test is \(2(1001/1000)^{101}-(1001/1000)^{200}-1<0\).

For two bins of duration three, \(h=(2,2)\), \(c=d=1\), \(\rho=\log(11/10)\), every compatible profile is indefinitely feasible with sufficient capacity and the exact robust energy is \(21/(100\rho)\). The least favorable full-cycle test begins **inside** a bin; its minimum multiplied by \(\rho\) is \(156639/10^6\), whereas both boundary-only full-cycle cuts give \(207459/10^6\). A fine chronology with neighboring deficits across the boundary attains the finite storage maximum. Lowering only the charge efficiency to \(3/5\), and keeping all these raw totals, bounds and standing losses, makes the boundary-only robust signs **misclassify feasibility**: both boundary cuts still give \(11655/10^6>0\) after multiplication by \(\rho\), but the interior full-cycle minimum is \(-29001/10^6<0\). No finite capacity then guarantees operation. The exact offline replay independently enumerates compatible unit-slot chronologies and simulates the greedy store with exact one-step decay; the continuous polygon proof is what makes the finite replay certifying for every measurable completion under these particular integer-cut controls.

These examples are synthetic exact model checks, not estimates of storage requirements, outage vulnerability, or actual metering sufficiency. Fixed converter bounds are honored in the model; fixed ramps, leakage depending on time/state, multiunit interaction, forecast/probability models, and nonrepeating missing chronology would require new work.
