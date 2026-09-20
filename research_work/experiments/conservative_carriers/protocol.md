# SV-1: reciprocal moving-carrier attraction, before calculation
Declared 19 September 2026 after OP-1 failed its galaxy-transfer criterion.

## Purpose and boundary
Test the owner's moving-carrier attraction directly through a bounded conservative mechanical proxy. It is deliberately an equal-mass, low-speed particle Hamiltonian, not a massless graviton theory, photon lensing prediction, newly sourced cosmological field or galaxy fit. Its instantaneous interaction needs a causal mediator before relativistic use. No dark matter or expansion is introduced.

A prepared packet is released from an ordinary-source reservoir. The release/preparation cost is calculated against a specified assembled reference, including source recoil and compensating spin. Carrier creation, a real matter fuel budget and observed galaxy normalization remain unmodeled. Rest masses are unchanged between the reference and prepared state and cancel in this energy difference. No success may be described as complete photon-to-field production or proof of gravity.

## Hamiltonian and derived interaction
Dimensionless units: each carrier inertial mass m=1; source remnant M=100; I_source=100; softening epsilon=0.2; N=12.
For positions q_i, canonical momenta p_i, source R,P:
w_ij=exp[-|q_i-q_j|^2/(2 ell^2)], w_ii=0.
k=eta/(N-1).
H = (1/2) sum_i |p_i|^2 - k sum_(i<j) w_ij p_i dot p_j
    + |P|^2/(2M) - mu sum_i 1/sqrt(|q_i-R|^2+epsilon^2).
Source spin is a constant supplementary reservoir: S=-sum_i q_i cross p_i initially, E_spin=|S|^2/(2I). Initial P=-sum_i p_i. Global initial angular and linear momentum thus vanish (R=0).

Hamilton equations:
qdot_i = p_i - k sum_j w_ij p_j
pdot_i = -mu (q_i-R)/(|q_i-R|^2+epsilon^2)^(3/2)
         - (k/ell^2) sum_j w_ij (p_i dot p_j)(q_i-q_j)
Rdot=P/M; Pdot=+mu sum_i (q_i-R)/(... )^(3/2).

Positive eta produces an attractive pair contribution when canonical momenta are parallel; negative eta reverses it. Canonical momentum is NOT physical velocity, so force signs alone do not establish attraction of trajectories. Measure the actual separation acceleration and separation histories. Angular momentum in the conserved ledger uses canonical momenta. Do not confuse source-imposed rotation with spontaneously created total angular momentum.

Since each row sum of w is <=N-1 and |eta|<=0.8, the carrier kinetic matrix has eigenvalues >=0.2. H is bounded below by -N*mu/epsilon. The normalization by N-1 is a finite-packet mathematical choice; a local, particle-number-independent continuum law is NOT derived. The interaction selects a frame and has no Lorentz-invariance claim.

## Fixed scan and initialization
252 runs: eta={-0.8,-0.4,0,0.2,0.4,0.6,0.8}; ell={0.25,1,3}; mu={0,1}; six packets:
- radial ring: r=1, p=0.8 e_r;
- rotating ring: p=0.8 e_r+0.6 e_phi;
- reverse ring: p=0.8 e_r-0.6 e_phi;
- counter-rotating ring: alternating +/-0.6 e_phi plus 0.8 e_r;
- hot sphere: 12 Fibonacci-sphere positions at r=1, isotropic seeded momentum directions of magnitude 0.9;
- stream: q_x=-3+0.12 i, q_y=0.7, q_z=0; all p=(0.9,0,0).
Deterministic seed 20260919, no chosen noise added to the symmetric packets. The source recoils, including in the mu=0 control. eta=0 is the matched ordinary-potential control for each geometry; mu=0 isolates the carrier interaction without ordinary binding.

Use DOP853 to t=20, rtol=1e-8, atol=1e-10, max_step=0.1, 101 stored times. Refine eta={-0.8,0,0.8}, ell=1, both mu and all six packets (36 runs), with rtol=1e-10, atol=1e-12 and max_step=0.05. No case chosen from favorable outcomes.

## Controls, outputs, criteria
Before trajectories: finite-difference Hamilton gradients for q,p,R,P (relative <=1e-6), pair reciprocity, translation/rotation invariance, zero-coupling ordinary force, kinetic lower bound, sign and heading controls. Separately check two-carrier physical separation acceleration using the RHS derivative, with no ordinary source.
Save every run and integration failure. Record energy drift normalized by initial |kinetic|+|potential|+|recoil| or 1; total canonical momentum and angular momentum; minimum kinetic eigenvalue at stored times; initial preparation cost H_initial-H_assembled+E_spin (assembled positions coincide with source, momenta zero). This cost is a required input, not evidence that real matter supplies it.

Metrics at every stored time: RMS source-relative radius, RMS centroid-relative packet width, fraction within source radius 3, normalized mechanical-velocity circulation, local weighted heading correlation and its total neighbor weight. Report low neighbor weights; no-neighbor correlation is set to zero and cannot establish alignment. Primary descriptive comparisons use the final 20 stored times against matched eta=0; they are synthetic proxies, not observation scores.

Numerical gates: scaled energy drift <=1e-6; total momentum drift <=1e-7; angular-momentum drift <=1e-6; minimum kinetic eigenvalue >=0.2-1e-12. Refinement compares whole trajectories by maximum absolute position difference normalized by max(1, maximum refined position norm), threshold 1e-3. Retain any chaotic/nonconverged cases, without changing tolerances post hoc. A positive mechanical finding requires a resolved nonzero-coupling effect against its ordinary-potential control; it does not establish useful gravity or clean lensing.

## Attribution and next dependency
Hamiltonian mechanics, Gaussian kernels, reciprocal pair interactions, matrix positivity and adaptive integration are established tools. Velocity/momentum-dependent attraction also has precedents; this finite-packet combination is a project test, not a claim of historical originality. Record specific literature in the report. A favorable case would justify deriving a local field mediator and matter/light coupling, not attaching an arbitrary optical multiplier. OP-1's observed-data failures and the existing CMF constraints remain unchanged.
