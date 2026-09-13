# Finite local formation histories

13 September 2026. Constant-field diagnostic with empty initial storage; no physical age imposed.

**Normalization correction (13 September 2026).** In the density, capacity and rate equations in this report, C denotes the pre-retention amplitude A=2 C0=9.457178483e7 Msun/kpc^3, where C0=4.728589242e7 is the stored fit parameter. C/original C multipliers are unchanged because the factor of two cancels. See [normalization audit](capacity-normalization-report.md).

## Finding

Finite buildup reduces local filling, but the tested shared durations do not outperform the identically amplitude-adjusted reference across the existing galaxy partitions. Training selection among six finite durations and equilibrium chooses u=10. Its validation/test velocity RMSE is 29.41/25.93 km/s versus 30.82/21.41 for the matched reference. The validation advantage is accompanied by a worse test result and a worse training score. This is not an adopted replacement or a new blind prediction.

This experiment isolates local formation in the original full-opacity field. It is not a time-dependent solution of the previous coupled transport system. It also differs from the earlier threshold-history calculation, which applied one filling factor per galaxy rather than allowing the input to vary spatially as XJ(r).

## Postulates and known solution

Retain the hypothetical capacity density Cg(r), g=[1+(r/a)^2]^-2, and hold the incident field at the full-opacity angular average J0(r) from kappa=k0 g. Use the prior normalized capacity-weighted threshold distribution on 1e-6<=t<=1e6. Its positive weights are the finite-cutoff version of the inverse-chosen exact-third mixture, not a first-principles derivation of one-third.

Set Y(r)=XJ0(r). Begin with empty threshold states. For dimensionless duration u=lambda Delta t, the known local rate equation and solution are

\[
\frac{df_t}{du}=Y(1-f_t)-t f_t,\qquad
f_t(Y,u)=\frac{Y}{Y+t}\left[1-e^{-(Y+t)u}\right].
\]

The deposited density is

\[
\rho_d(r,u)=Cg(r)\int f_t[Y(r),u]\rho(t)dt.
\]

Application to companions, the capacity Cg, the threshold weights and loading proxy are hypotheses. The exponential solution is established linear-rate mathematics. A common u assumes the same microscopic lambda and the same exposure duration in these controlled cases; it is not a claim that real galaxies share formation histories. Neither lambda nor a cosmic age is fixed.

As u tends to zero at fixed finite threshold cutoffs, f_t/u tends to Y. This initially recovers linear dependence on illumination, but with an amplitude proportional to X u rather than the reference eta(X). A single common amplitude cannot generally remove that different dependence on each galaxy's source proxy. At large u the density approaches the finite-cutoff local-equilibrium branch. It need not reproduce the exact infinite-support local formula perfectly; that distinction is retained.

## Controlled comparison

Predeclare the diagnostic grid u=0.01,0.1,1,10,100,1000 and equilibrium. Freeze the existing k0, a/R_d, observed inputs and finite threshold range. For each duration, fit only the common C multiplier by equal-galaxy velocity MSE on the 89 training galaxies. Apply the same amplitude adjustment to the reference. The amplitude search expands its upper bound when required, rather than capping early-time solutions at an arbitrary low value; every optimum is interior. There are no per-galaxy fitted durations or amplitudes.

| Duration u | C/original C | Training RMSE | Validation RMSE | Test RMSE |
|---|---:|---:|---:|---:|
| Reference control | 0.747109 | 26.87 | 30.82 | 21.41 |
| 0.01 | 10.091649 | 31.30 | 39.92 | 29.78 |
| 0.1 | 1.505178 | 30.43 | 37.98 | 29.07 |
| 1 | 0.505047 | 29.38 | 33.69 | 27.01 |
| 10 (training-selected grid duration) | 0.372963 | 28.55 | 29.41 | 25.93 |
| 100 | 0.354411 | 28.74 | 28.66 | 25.32 |
| 1000 | 0.353414 | 28.71 | 28.65 | 25.21 |
| Equilibrium, finite thresholds | 0.353413 | 28.71 | 28.65 | 25.15 |

Scores are km/s with equal weight per galaxy before taking the square root. Existing partitions contain 89/29/31 galaxies and have already been repeatedly exposed. The selected duration is a shared explored model choice in addition to the fitted amplitude; it is not a parameter-free prediction. Observational and baryonic uncertainties are not propagated. Frozen-amplitude scores and per-galaxy adjusted residuals are recorded in the results JSON as well.

The lowest training score among these duration branches occurs at u=10, but it remains above the reference's training error. Selecting equilibrium or another duration afterward because of test results would change the selection rule. No duration is adopted. This finite grid does not prove that every possible duration or history fails.

## Relation to coupled transport

For the empty-site opacity kappa=k0 g(1-f), opacity never exceeds the original k0 g. Under a maintained isotropic boundary and an initial field at least J0, reducing opacity or adding a nonnegative return source cannot make the incident field smaller than J0 along otherwise unchanged rays. Local capture/release dynamics are increasing in the supplied Y while 0<=f<=1. Thus the fixed-field filling provides a lower-field comparison for that class of consistent histories.

This ordering concerns occupancy with unchanged capacities, rates and source histories. It does not imply ordered velocity residuals, nor does it remain a direct fitted-mass bound after changing C or the physical rate mapping. Actual illumination changes as sites fill, and finite propagation delays must be accounted for in a time-dependent model. A constant J0 prescription is therefore used only as a discriminator, not quietly combined with the previously solved stationary feedback field.

## Energy accounting

For a threshold k=Y+t and equilibrium fraction f_eq=Y/k, the integrated local capture and release energies per capacity are

\[
C_t=Y\left[(1-f_{eq})u+\frac{f_{eq}(1-e^{-ku})}{k}\right],
\]

\[
R_t=t f_{eq}\left[u-\frac{1-e^{-ku}}{k}\right],\qquad C_t-R_t=f_t(Y,u).
\]

These known rate integrals preserve local energy for the prescribed field. They do not supply a global incident-power calculation or a physical release channel. Large early-time C adjustments represent larger capacities, not free energy; connecting them to supply and lambda remains necessary. No claim of permanent retention follows from having a finite formation duration.

## Consequence

The short-time linear-illumination limit is a useful successful mathematical limit, but its galaxy-to-galaxy amplitude scaling differs from the reference. The tested shared histories do not repair that difference with one amplitude. A measured or independently derived history dependence could change the outcome; a freely chosen duration for every galaxy would merely add flexibility. The reference, prior energy conclusions and unresolved propagation/lensing/support requirements remain intact.

## Reproduction

Run `python research_work/results/companion-extensions/local-formation.py`. It verifies raw-data and model hashes, all 149 accepted galaxy arrays and 3150 radii/velocities, baseline scores and monotonic local filling. It compares radial/angular quadrature orders 64/128 and threshold orders 128/256 for every duration and galaxy; the largest frozen-amplitude speed difference is 1.27e-9 km/s. These are integration checks, not independent physical validation. Results include both fixed and adjusted scores, per-galaxy adjusted errors, and the duration chosen using training alone.

Related: [local equilibrium capacity](local-capacity-report.md), [coupled stationary transport](occupancy-transport-report.md), [earlier global-factor histories](threshold-history-report.md).
