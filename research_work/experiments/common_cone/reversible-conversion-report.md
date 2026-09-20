# RC-1: reversible conversion and environmental response

20 September 2026. A new effective interaction test in a nonexpanding setting, with no dark-matter component. Judged by its own Hamiltonian and conservation identities, not agreement with an older gravity formula. No observational fit or ordinary spin-2 identification is claimed.

## Direct answer to the proposed mechanism

Reversible conversion can be implemented without creating energy. In this local toy model, changing the relative response of the two channels also changes the environment's motion through an explicitly reciprocal force. Simply converting between identically responding channels does not give stronger gravity. Conversion, temporary binding, and subsequent travel on different paths are distinct steps; only local conversion and its environmental reaction are tested here.

Use the Hermitian matrix M(Q)=[[1+delta(Q),kappa(Q)],[kappa(Q),1-delta(Q)]], amplitudes psi=(a,b), and

H = P^2/2 + .2^2 Q^2/2 + psi^dagger M(Q) psi.

i psi_dot=M psi; Qdot=P; Pdot=-.2^2 Q-psi^dagger M'(Q) psi.

Here delta=.0 or .3 times tanh Q and kappa=kappa0(1+A tanh Q). The same interaction supplies conversion in both directions and the environmental force. A local state can influence the rate; no desired orbit or lens map is provided to it. The environment is one oscillator, not a completed spatial matter/field system. Its generalized force must not be advertised as a verified 3D momentum or angular-momentum budget.

## Executed evidence

- 100 finite-difference source-force/positive-matrix controls pass.
- First archive: 23/24 main trajectories meet every declared gate; run-19 fails backward-state accuracy (5.24614e-7 versus 1e-7). The failure is preserved.
- With unchanged equations and gates, stricter tolerances give 24/24 passes; all seven analytic/no-conversion controls and twelve deliberately missing-reaction controls pass their intended tests.
- Maximum relative energy drift 6.16013e-12; norm drift 4.38338e-12; backward-state error 3.06757e-09.
- Highest sampled receiving-mode population 0.999999876089. Constant-rate tests recover sin^2(kappa t), including return conversion. This is normalized mode population, not extra energy created or a guaranteed microscopic conversion probability.
- Identical-channel cases change the environmental coordinate from their no-conversion controls by at most 3.81817e-12. Their interaction-force expectation is zero for the declared pure initial state.
- Different-channel cases change that coordinate by as much as 6.64201 model units. These are response changes, not stable orbits or improved fits. Positive local mode energy does not prove full nonlinear stability.
- Omitting the environmental reaction produces relative energy errors from 0.207653 to 0.690519. Counting only the two mode energies and ignoring signed interaction energy is also invalid.
- Independent raw-state reconstruction: 320 audit checks pass; final-energy reconstruction difference <= 2.22045e-16.

## Enough gravitational response?

For a bookkeeping diagnostic with energy fraction f and derived active coefficients chi_gamma,chi_g:

M_active = [(1-f)chi_gamma+f chi_g] E_stored/c^2,

before any additional interaction source is included. In SR-1's weak branch, both coefficients are 1+s. The photon Hamiltonian derivative and the averaged wave-energy derivative therefore give no leading gain from conversion alone. This is a conditional result about SR-1, not a prohibition on a different invented interaction.

For the prior illustrative 200 km/s, R=6e20 m, M=2e41 kg budget, keeping only epsilon=.001 of source rest energy available requires an effective coefficient about 1780 in that same equivalent-energy interpretation, compared with SR-1's 2. The new theory could instead derive a different constitutive/interaction force; then its own force and energy relation must replace this budget. Choosing 1780 by hand does not explain observations. The budget's f is an energy fraction and must not be confused with the mode population when diagonal energies and interaction energy differ.

## Attachment and separate outgoing paths

The owner's extension is retained as the next spatial mechanism: a companion may couple to matter/radiation, travel in a bound or coherent state, and detach with a different outgoing direction. A complete implementation must derive binding energy, capture/release dynamics, the outgoing momenta and field/source reaction from one interaction. Deterministic or coherent release can in principle yield repeatable bending; random scatter is not assumed just because detachment occurs. Image-width and spectral predictions must decide which behavior the law actually produces.

Within any retained standard massless energy-momentum kinematics, an isolated massless particle cannot split into two non-collinear positive-energy massless particles without another participant: (p1+p2)^2=2 E1 E2(1-cos theta)/c^2. Matter, a finite-mass bound state or the surrounding field can supply the required balance. This is a stated kinematic assumption, not an old gravity law used as a veto. If the candidate changes it, specify and test the replacement conservation/propagation rules. No such bound state or spatial splitting dynamics is demonstrated by RC-1.

The next model must derive the field sourced by matter, spatial following and release, common signal propagation, an affordable long-term budget and both stellar and photon trajectories. Then test direct observations with shared constants. The current oscillator calculation supplies none of the missing galaxy/cluster profiles.

## Attribution and reproduction

Hermitian mode mixing and Hamiltonian reaction terms are established mathematical tools. Photon/graviton conversion has Gertsenshtein precedent; [Palessandro and Rothman](https://arxiv.org/abs/2301.02072) provide a derivation. Their magnetic-field mechanism and rates are not imported here. The local Q-dependent specialization is our declared effective candidate, not a claim of historical uniqueness. Existing CWC-1 radiation-to-scalar tests remain separate. [Weinberg's universal-coupling result](https://doi.org/10.1103/PhysRev.135.B1049) applies under stated framework assumptions; it is attribution/context, not a rejection gate for a different invented framework.

Use a separate checkout at b58ce4e, before the evidence archives are committed, and run reversible_conversion.py followed by refine_reversible_conversion.py. The first command intentionally exits unsuccessfully after preserving its complete failed first archive; proceed to the declared refinement. Run audit_reversible_conversion.py from the completed campaign against either reproduced or checked-in evidence. Source hashes and commits are preserved in both manifests.
