# Finite-response temporal energy: what improved and what failed

## Result

A finite response time separates the companion energy from the energy that changes the clock field. This avoids the preceding algebraic time-matrix singularity. A slower, but still traveling, temporal sector also avoids the tested equal-speed degeneracy. However, neither version is a complete candidate: the proposed exchange rule can require negative companion energy in a spatially varying field. A safeguard that fixes uniform-cell accounting does not fix that spatial failure.

In ordinary language, a passing total-energy sum is insufficient if one account can go below zero. We need to specify which physical sector pays for radiation gaining energy as well as which receives radiation losing energy. This is an operational rule we can hypothesize; it does not require explaining the ultimate origin of time.

## Explicit new postulates

R is combined electromagnetic and gravitational-wave reference-energy density; T is traveling companion energy; M is traveling temporal-sector energy; D is deposited energy. K is a positive energy-density scale, tau a positive response time, chi a nonnegative dimensionless release coefficient and Gamma a nonnegative capture rate. Reference time and length are nondimensionalized with c0=1.

**Project postulates; no originality claim:**

\[
n=1+M/K,\quad v=1/n,\quad d\tau_{\rm clock}=dt/n,\quad b=R/(Kn),
\]

\[
\begin{aligned}
R_t+(vR)_x&=-bM_t,\\
T_t+(vT)_x&=bM_t-Q-\Gamma T,\\
M_t+(\sigma vM)_x&=Q,\\
D_t&=\Gamma T.
\end{aligned}
\]

The response timescale tau is distinct from clock proper time tau_clock. The original exchange is Q=(T-chi M)/tau. The separately specified availability variant is

\[
Q=\frac{T-\chi M\max(0,1-b)}{\tau}.
\]

This availability factor regulates an energy exchange, not a spatial screening of gravity or redshift. Its microscopic origin is not derived. M moves with flux sigma vM; it is not free stationary memory. Photons, GWs and T share local c0. Sigma=1 retains that speed for M's energy flux; sigma=0.5 is an explicit alternative with slower temporal-sector energy, not an adopted change to the user's preferred companion speed.

**Conditional derivation, ordinary conservation algebra:**

\[
\partial_t(R+T+M+D)+\partial_x[v(R+T)+\sigma vM]=0.
\]

This remains a reference-energy ledger, not a complete physical stress-energy or momentum theory.

## Local propagation

Eliminating M_t gives a standard first-order principal matrix in variables (R,T,M):

\[
A=\begin{pmatrix}
1/n&0&-R/(Kn^2)-R\sigma/(Kn^3)\\
0&1/n&-T/(Kn^2)+R\sigma/(Kn^3)\\
0&0&\sigma/n^2
\end{pmatrix}.
\]

**Known characteristic mathematics applied to these postulates:** its coordinate speeds are 1/n, 1/n and sigma/n squared. The corresponding clock-measured speeds relative to c0 are 1, 1 and sigma/n. The former singular denominator 1-R/(Kn) is absent.

For sigma=1 at M=0, the repeated-speed matrix is generically defective. Finite response does not repair this principal high-frequency problem. An independent Fourier matrix-exponential calculation shows normalized gains increasing from 1.397 at wavenumber 1 to 975.09 at wavenumber 1000 for the declared background. These are infinitesimal linear gains, not finite physical amplitudes.

For 0<sigma<1 and nonnegative M, the M characteristic is strictly separated from the repeated radiation/T speed. The displayed matrix then has a complete eigenbasis. All 27 sigma=0.5 grid states verify this; their largest normalized eigenvector condition number is 126.25, so diagonalizability does not mean small transient amplification. This is a local principal-part result, not proof of nonlinear stability or a successful redshift prediction.

## Uniform cells and the empty-reservoir test

At T=0 in a uniform cell, the raw exchange gives T_t=(1-b)chi M/tau, which is negative when b>1. Nine of 36 raw runs hit the declared negative-T stopping threshold. All 36 availability-regulated runs remain nonnegative and conserve the reference sum through reference time 20. Their homogeneous Rn invariant also passes. Three exact no-seed controls remain inactive: this rule does not initiate conversion from radiation alone.

But at T=0 with a spatial gradient in M, the regulated equations instead give

\[
\boxed{T_t=-\frac{b\sigma}{n^2}M_x+(b-1)Q.}
\]

The first term can overwhelm the positive local release. At K=0.5, R=0.22, M=0.1, T=0, sigma=0.5 and chi=tau=Gamma=1:

| M spatial gradient | Predicted T time derivative |
|---:|---:|
| 0 | +0.040111111 |
| 0.1 | +0.027379630 |
| 1 | **-0.087203704** |

A smooth positive M profile can have this local gradient while T is identically zero nearby. The last row therefore exits the allowed nonnegative-energy domain immediately. It is an analytical counterexample, independently evaluated from the original and reduced equations; no numerical clipping or long simulation is involved. A model-specific restriction on gradients would need a physical justification and a proof that evolution preserves it.

Thus the local response safeguard is insufficient. Spatial advection changes n, which changes the radiation reference energy; assigning the entire balancing exchange to T can demand energy from an empty T sector. Adding a response timescale does not by itself resolve that assignment.

## Verification, prediction limits and next step

The [protocol](protocol.md) distinguishes original tests from the later spatial-boundary check. [run.py](run.py) generates all [54 characteristic states](characteristics.csv), [72 cell runs](cells.csv) and [results](results.json). Run:

```text
python research_work/results/finite-response-transport/run.py
```

Maximum reference-energy error is 4.67e-15; scaled homogeneous Rn error is 3.45e-13. Fourier checks agree within 2.85e-13. Both independent spatial-boundary calculations agree to floating-point precision. Failed raw runs remain archived.

Homogeneous time-field evolution still cancels in the stipulated source/observer clock comparison: its measured spectral and infinitesimal event-stretch ratios are one. Changing reference photon energy alone is not an observed redshift. No spatial redshift solution, observation fit, holdout, momentum completion, gravitational force or lensing prediction is established by these tests.

The next necessary derivation is a spatially valid allocation of radiation/field work to actual energy-carrying sectors, with nonnegative budgets at empty-sector boundaries. Recheck its characteristic speeds before any redshift fit. Simply clipping T to zero would break the ledger, and prescribing an independent positive redshift coefficient would bypass the selected physical mechanism. The complete nine-goal program remains open.
