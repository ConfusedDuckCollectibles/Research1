# A three-setting paired-count witness for a two-mode class-D reflection sign

## Statement and scope

Let a coherent, elastic, quadratic Bogoliubov–de Gennes (BdG) scatterer couple to two normal leads and a superconducting condensate. One monitored lead has exactly two resolved incoming and outgoing physical electron modes and their hole partners. Both normal leads may carry transmitted quasiparticles. In the particle-hole-related Nambu bases, the monitored reflection block at the zero-energy limit satisfies
\[
 r=\begin{pmatrix}A&\overline B\\ B&\overline A\end{pmatrix},\qquad A,B\in\mathbb C^{2\times2}.
\]
It need not be unitary. The zero-energy finite-device reflection sign \(\operatorname{sgn}\det r\) is real, and is undefined as a binary classification at a zero determinant. Obtain calibrated access to the two outward orbitals with a number-conserving mixer \(V\) on electrons and its conjugate \(\overline V\) on holes. The analyzer changes the measurement basis, not the scattering boundary conditions of the device. A reflectionless normal mixer may rotate the inputs as well: any unitary on the fully occupied two-dimensional input subspace introduces only its Slater determinant phase. Independent normal input/output basis changes have even real orientation (\(\det V\det\overline V=1\)); hence the sign requires no absolute condensate phase or coherent mixing of an electron with a hole. It does require the correctly paired basis and calibrated **relative orbital control**, two distinct notions.

Prepare one electron excitation in each of the two occupied orthogonal incoming modes, with the same spectral wave packet, no other input excitations and Pauli-antisymmetric two-particle state. As described below, all zero-energy expressions refer to the continuous on-shell limit from positive excitation energy, not to treating a state at a Fermi-surface singularity as two independent electron/hole source modes.

Select uniformly among three normal output bases, the \(x,y,z\) axes of the two orbital modes. Score
- \(+1\) if both excitations are reflected into the monitored lead, **except** the next case;
- \(-2\) if they are a reflected electron and hole in *different* output orbitals;
- \(0\) if fewer than two excitations are reflected into that lead.

**The exact unconditional expectation of this score is \(\det r\)**, including when there is appreciable transmission out of the monitored lead. No conditioning on the trials with two returns is permitted for this equality. Two normal settings, even with the whole local outcome distribution in each and with exact data, do not suffice in general under this same fixed-input experiment. All individually injected, channel-resolved electron/hole probabilities in both leads can also be identical for opposite signs; this is proved with a separate no-hidden-channel pair.

This measures the sign of the **finite device reflection block**. As distinguished by [early scattering theory](https://arxiv.org/abs/1101.1749) and [reflection-zero theory for finite devices, November 2026](https://arxiv.org/abs/2611.17001v1), it does not itself establish topological protection, spatially separated bulk-associated Majorana modes, or that every trivial Andreev description is excluded.

## 1. Determinant from a two-electron occupation experiment

Use order \(e_1,e_2,h_1,h_2\) for the local outgoing channels. Let \(U=\binom A B\) be the two columns from the two occupied incoming electron states. Let \(p_{ij}\) for rows \(i<j\) be the \(2\times2\) minor of \(U\). Fermionic antisymmetry gives the outgoing pair amplitude \(p_{ij}\). For the full unitary two-lead scatterer the output is a normalized Slater determinant; retaining the minors of the local rows means that each \(|p_{ij}|^2\) is an **absolute** local two-return probability. It is not divided by return rate. Events with no local return are still counted as trials.

Set \(P_{ee}=|p_{12}|^2\), \(P_{hh}=|p_{34}|^2\), \(M_{ij}=p_{i,2+j}\) for \(i,j\in\{1,2\}\), and \(P_{eh}=\|M\|_F^2\). If \(C\) exchanges electron and hole rows, then \(r=[U,C\overline U]\). Laplace expansion through the two columns of \(U\) gives
\[
 \begin{split}
 \det r&=|p_{12}|^2+|p_{34}|^2-|p_{14}|^2-|p_{23}|^2
                   +2\operatorname{Re}(p_{13}\overline{p_{24}})\\
       &=P_{ee}+P_{hh}-P_{eh}+|\operatorname{tr}M|^2.                \tag{1}
 \end{split}
\]
In the Laplace sum a mixed complement has its two rows in the opposite order after toggling electron/hole, giving the extra minus sign; the same-charge complements have none. Formula (1) derives the sign without assuming \(r\) is near unitary.

The same-charge probabilities and \(P_{eh}\) are invariant under any lossless normal output orbital rotation. In a rotated basis \(M\) transforms as \(V M\overline V^T=V M V^\dagger\). For setting \(a\in\{x,y,z\}\), choose \(V_a^\dagger\sigma_z V_a=\sigma_a\), where \(\sigma_a\) are the Pauli matrices in the orbital space. The absolute probability \(S_a\) of a local electron-hole pair in the *same* resulting orbital is
\[
 S_a=\tfrac12\big(|\operatorname{tr}M|^2+|\operatorname{tr}(\sigma_aM)|^2\big).
\]
The Hilbert–Schmidt orthogonal basis \(I,\sigma_x,\sigma_y,\sigma_z\) gives
\[
 S_x+S_y+S_z=\|M\|_F^2+|\operatorname{tr}M|^2.
\]
Consequently the pair benchmark is
\[
 \boxed{\ \det r=P_{ee}+P_{hh}-2P_{eh}+S_x+S_y+S_z
       =P_{\rm 2\,returns}-3\,\mathbb E_{a}[P_{\rm mixed,different}(a)].\ }        \tag{2}
\]
This is the proposed mean score. At a given setting the second formula's rare mixed/different event is counted unconditionally. The first term constitutes a local two-excitation **number** probability, not a mean net charge. One cannot replace it by a net-charge readout that confuses an electron and hole in the same orbital with no return. Alternatives to sign-and-mode-resolved counting would need extra calibrated observables sufficient to recover that number probability (for example an appropriate energy-count observable in an energy-resolved idealization). The work does not assert that ordinary DC current measurement alone implements (2).

**If a charge counter cannot resolve a same-orbital neutral pair.** Keep separate net charges \(q_i=n_{ei}-n_{hi}\) for the two measured output orbitals in the paired trials. Then \(P_1=\Pr(|q_1+q_2|=1)\) is exactly the probability of a single local return, and \(P_{\rm mixed,different}=\Pr(q_1q_2=-1)\). Separately inject one electron at a time in input \(e_1\) and input \(e_2\) against the same vacuum reference and same calibrated scatterer. Their individual absolute return rates sum to \(\mu=E[N_{\rm reflected}]\) in the paired preparation: the Slater one-body density is the sum of their two one-body projectors. Thus \(P_{\rm 2\,returns}=(\mu-P_1)/2\); use this in (2), averaging \(P_1\) across settings if desired. With at most two excitations there is no ambiguity between zero and two in the *single-injection* controls. This alternative requires calibrated separate preparations and synchronized mode-resolved charge counting; it is not just a conductance measurement, and its own finite-sample error must include the uncertainty in the single-injection calibrations.

## 2. Exactly what lesser measurements fail to determine

### A. Full two-lead single-particle probability ambiguity, with no missing mode

Use two electron modes in each lead, \(c=4/5\), \(s=3/5\), and \(P=\begin{psmallmatrix}0&1\\1&0\end{psmallmatrix}\). In a real Majorana basis with two quadrature sectors, take
\[
 K=\begin{pmatrix}cI&sI\\-sI&cI\end{pmatrix},\qquad
 J_0=\begin{pmatrix}cP&sP\\-sP&cP\end{pmatrix},\qquad
 J_1=\operatorname{diag}(-1,1,-1,1) J_0 .
\]
The full scattering matrix in real quadrature order is \(\operatorname{diag}(K,J_i)\). Each block is real orthogonal; the full determinants are \(+1\) for both choices. Change to all-electron-before-all-hole Nambu ordering: \(A_{\rm full}=(K+J_i)/2\), \(B_{\rm full}=(K-J_i)/2\); the sign convention for the second formula is immaterial and fixed in both models. The supports of \(K\) and each \(J_i\) are disjoint entrywise, so both devices give *exactly* the same squared modulus of **every** entry of the full Nambu scattering matrix. In particular every ordinary fixed-basis electron or hole injection has the same per-channel outcome probabilities in either lead, all aggregate charge/linear conductance matrices match, and the zero-energy transmission probabilities entering elastic low-temperature heat response are indistinguishable. This is at zero energy; it is not a statement about all finite-bias spectra or higher-count data in the original basis.

For left reflection, in the real quadratures the blocks are \(cI\) and \(cP\) for device 0, and \(cI\) and \(c\,\operatorname{diag}(-1,1)P\) for device 1. Their determinants are respectively \(-c^4\) and \(+c^4\), i.e. \(-256/625\) versus \(+256/625\). Each has single-particle transmission probability \(s^2=9/25\) from any monitored normal mode, and (2) resolves the sign. There is no unobserved-mode substitution in these equal-data claims.

### B. Even full local coincidence data in *any two fixed* analyzer bases can fail

A normal analyzer basis corresponds to one Bloch axis of the orbital space (its output phases are irrelevant to probabilities). For any two settings, pick an axis normal to both and orthonormal spinors \(u_+,u_-\) oriented along its two signs. In both chosen bases, their two orbital amplitudes have equal modulus \(1/\sqrt2\). Choose \(a=2/5,b=7/12\) and define two reflections with the **same** \(A=aI\). Let
\[
 M_+=ab\,|u_+\rangle\langle u_+|,\qquad
 M_-=ab\,|u_+\rangle\langle u_-|.
\]
They are the actual mixed Plücker-amplitude matrices for \(B_\pm\) defined by
\[
 B_\pm^T=-J M_\pm/a,\qquad J=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]
Indeed \(M=a J B^T\) when \(A=aI\). \(B_\pm\) are rank one with operator norm \(b\). The associated \(r_\pm\) are contractions since \(\|r_\pm\|\le a+b=59/60<1\); they are completed to physical two-lead unitaries below.

In either observed setting, all mixed pair probabilities are \(a^2 b^2/4\), and \(P_{ee}=a^4,P_{hh}=0\) for both. The single-output marginal occupation probability is the corresponding diagonal of the rotated \(U U^\dagger\). It is \(a^2\) in each electron channel, and \(b^2/2\) in each hole channel because \(B_\pm B_\pm^\dagger=b^2|\overline{u_\pm}\rangle\langle\overline{u_\pm}|\). Therefore the entire unconditional local distribution with zero, one or two returns is identical, not just the three pair statistics. Explicitly each probability of exactly one local detector is its marginal minus the sum of its pair probabilities, and the probability of zero returns is the remaining mass. The two models continue to match all local electron- or hole-injection single-particle probabilities measured in those two settings: the differing rank-one output spinors have equal projection intensities in both settings; the input covector of \(B\) is the same. No equality of their right-lead or arbitrary three-setting data is asserted here.

By (1) the determinants are \(a^4=+16/625\) and \(a^4-a^2b^2=-649/22500\). Thus three is the smallest universally sufficient number of fixed normal output bases for the specified two-electron input and **local** detection family. The same obstruction applies when a deterministic rule chooses its second basis from the exact first local distribution: that first distribution is the same for every normal axis perpendicular to the first setting, so select a normal to both axes after identifying the rule's second choice. This does not bound experiments with coherent input and full process tomography, multiple lead/count arrangements, or phase references. For example coherent electronic inputs with arbitrary output tomography can determine each one-particle map \(A,B\) up to its overall phase; those two phases can be absorbed by independent uniform normal input/output phases, leaving \(\det r\) unchanged.

### Physical scattering realization

The constructions are admissible zero-energy class-D scattering, rather than merely real arrays of probabilities. Transform a particle-hole-symmetric contraction \(r\) to its real Majorana version \(r_M\). In the two leads of four real channels each, an admissible completion is
\[
 S_{\!M}=\begin{pmatrix}r_M&\sqrt{I-r_Mr_M^T}\\
                 -\sqrt{I-r_M^Tr_M}&r_M^T\end{pmatrix}.
\]
The roots are the positive square roots. The singular value decomposition of \(r_M\) shows orthogonality by splitting into four real orthogonal \(2\times2\) pairs. The full determinant is \(+1\), also by continuous deformation of \(r_M\) to zero. This includes contractions of norm 1; the strict examples have nonzero total transmission in every mode. No ancillary open detector modes outside the two stated leads are required. For strict \(r_M\), the dilation has no \(+1\) eigenvalue: writing \(D_*=\sqrt{I-r_Mr_M^T}\), \(D=\sqrt{I-r_M^Tr_M}\), a fixed point \((x,y)\) would have \(D_*y=(I-r_M)x\) and \(Dx=-(I-r_M^T)y\). Pairing with \(y,x\) gives \(y^TD_*y+x^TDx=0\), which forces \(x=y=0\). The explicit matrices in A also avoid \(+1\), since their lead rotation has eigenvalues \(c\pm is\) and their orbital factors have eigenvalues in \(\{1,-1,i,-i\}\).

For an explicit class-D causal energy-dependent model for any such \(S_{\!M}\in SO(8)\), omit its \(+1\) eigenspace (which can simply bypass) and let
\[
 D_{\rm int}=\gamma(S_{\!M}-I)^{-1}(S_{\!M}+I)
\]
on the orthogonal complement, where \(\gamma>0\); this is real antisymmetric. Couple a finite Majorana dot Hamiltonian \(H=iD_{\rm int}\) to the real on-shell lead channels of this complement with \(W=\sqrt{\gamma}I\) in the convention
\(S(E)=I-2iW^T(E-H+iWW^T)^{-1}W\), with identity on bypass. This is the finite-dimensional [Mahaux–Weidenmüller model in class D](https://arxiv.org/abs/1506.05995); the model has physical particle-hole symmetry \(S(-E)=\overline{S(E)}\) in the Majorana basis, is unitary at real \(E\), and at zero gives precisely \(S_{\!M}\). The internal dimension is even because every non-real conjugate pair is even-dimensional and the \(-1\) eigenspace of an \(SO(8)\) element has even dimension. With full-rank coupling on the retained subspace, there is no decoupled pole obstructing the near-zero limits. This shows no continuity or parity obstruction to these example matrices. It is an effective wide-band finite Hamiltonian realization, not evidence that a particular local materials sample realizes it. In particular a finite dot model with no bulk limit can have \(Q=-1\).

## 3. Finite-energy and calibration margin

At exactly mirrored nonzero excitation energies the two matrices are related by particle-hole symmetry, not individually real; use the one-sided on-shell limit \(S(0)\) in the formulas. A physical source can in principle supply two orthogonal orbital excitations with indistinguishable positive-energy spectral envelopes. Alternatively a formal cold input slice has precisely the two occupied incoming electron modes, but ordinary energy/time-integrated bias data are **not** automatically a sample of the required coincidence law. Assume the detector channels are the actual outgoing particle/hole excitations in a normal lead at positive energy with the corresponding offset sign in physical energy. A normal mixer must act on outgoing orbitals without changing \(r\) via reflection, contact loading, or setting-dependent interference. The \(e/h\) analyzers are conjugates at zero energy and must be calibrated across the actual energy window. A condensate phase reference is not used.

For the finite-bandwidth comparison, assume an independently justified uniform operator-norm bound \(\delta_s\) on \(S(E)-S(0)\) over the two inputs' supported energy window when regarded as operators on the full on-shell modes. Include an analogous \(\delta_v\) for the actual combined electron/hole normal analyzer relative to its target. The single-particle composition error is at most \(\delta_s+\delta_v\) by unitarity. For two identical-envelope inputs in orthogonal orbitals, the coherent evolution on wave packets is the exterior square of that single-particle operator; its operator norm changes by at most \(2(\delta_s+\delta_v)\) by expanding the difference on the tensor square and restricting. The ideal energy-independent comparison factors into the orbital problem and the common envelope, so counting without resolving energy leaves formula (2) valid for that comparison. This argument also accommodates an incoherent mixture of common narrow slice data obeying the same bound. It does **not** replace calibration of indistinguishability or require two particles to have an identical classical sampled energy in a finite wave packet.

If source plus finite-temperature input deviations from the ideal two-excitation sector have trace distance at most \(\eta\), and detector/reporting total variation is at most \(\zeta\), the mean score bias is bounded by
\[
 b_{\rm sys}\le 6(\delta_s+\delta_v)+3\eta+3\zeta.
\]
This uses score range \([-2,1]\), trace-distance/TV contraction and the fact that trace distance of pure states is at most their Euclidean vector distance. In a normal Fermi reservoir nonzero temperature changes the occupied/empty input sector near \(E=0\); small \(\eta\) must be certified by actual prepared states/operating window, not simply assumed from nominal zero bias. Distinguishable, interacting or incoherently mixed input orbitals are not in the exact contract, and if needed a separate source characterization must account for spectator excitations. For the explicit examples with full \(W=\sqrt\gamma I\), the resolvent has norm at most \(1/\gamma\) on the real axis so \(\|\partial_E S\|\le2/\gamma\). Thus a compact wave-packet support up to \(E_{\max}\) in that illustrative model gives \(\delta_s\le2E_{\max}/\gamma\), without estimating \(\gamma\) for any material. The operator bound models analyzer error and energy dependence where controlled, not changes in the number of physical modes, inelastic currents, missing channels or detector crosstalk without a TV calibration.

For \(n\) independent randomized-setting trials with a fixed expectation model, Hoeffding gives a two-sided score uncertainty of \(\sqrt{(9/2n)\log(2/\alpha)}\) with failure probability at most \(\alpha\). With a confidence interval for the actual mean, also subtract the systematic margin from the distance of the interval to zero. For illustration only, for one of the \(|\det r|=256/625\) examples, assuming \(\delta_s+\delta_v=\eta=\zeta=0.01\) and \(n=3000\) ideal trials at \(\alpha=0.05\), the worst-case total radius is \(<0.195\) (statistical \(<0.075\), systematic \(0.12\)). For an observed estimate, certification requires \(|\text{estimate}|>0.195\) with these conservative nuisance bounds. Since \(2(0.195)<256/625\), for this model and those fixed bounds the confidence event also guarantees that the threshold is met. No realized throughput, indistinguishability, temperature, device or detector errors are claimed. Near \(\det r=0\) the sign is inherently unstable. The example's exact det values are at a single fixed zero-energy contact tuning, not a junction phase diagram.
