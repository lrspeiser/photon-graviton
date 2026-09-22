# JR-9 protocol: same-object cross-prediction and one derived response kernel

**Declared 22 September 2026, before any cross-prediction score was computed.**
Baseline commit for the frozen response: `e276202` (JR-8 head). The R10 forward
model is reproduced bit-for-bit from `research_work/results/joint-response-iteration`
before anything else runs; that reproduction check is part of the executed output.

This is a development experiment on historically examined objects. It is not a
blind test, not a new observation, and not a completed microscopic derivation.

## Why this experiment exists

JR-1's R10 is the project's best joint fit, but its four optimizer-fitted lenses
had their stellar-mass offsets and orbital anisotropies adjusted against stellar
kinematics **and** the Einstein angle at the same time. That cannot answer the
question the hypothesis actually needs answered: if one companion state produces
both observables, then constraining it with one observable must predict the other
without further adjustment.

JR-9 replaces the joint adjustment with two one-directional predictions.

## A. Frozen universal response, per-lens ordinary nuisance only

The ten universal field/source coefficients and the two shared stellar-population
scales stay at their published R10 values and are never re-fitted here:

    logA, logc, p, dA, dc, q, logt, logsph, core_sph, q_sph, logu, logusph

No object receives an individual companion amplitude and no lens receives a
lensing gain. The only per-lens freedom is the ordinary astrophysical nuisance
pair JR-1 already declared:

| Symbol | Meaning | Prior | Bounds |
|---|---|---|---|
| `dm` | stellar log-mass offset, dex | N(0, published `log10_mass_error`) | +/- 5 sigma |
| `beta` | constant orbital anisotropy | N(0, 0.30) | (-1.00, 0.35) |

`dm` enters the baryonic force, the companion source proxy and both observables
consistently, exactly as in JR-1. The Einstein angle is independent of `beta`;
that is a property of the deflection integral, not an assumption added here.

All six systems are treated identically. JR-1's split, where only the four fitted
lenses carried individual offsets, is not reproduced.

## B. The two declared directions

**Direction K -> L (motions predict the lens).** For each lens, form the posterior
over `(dm, beta)` from the resolved KCWI V_rms bins with their full released
covariance, plus the two priors above. The Einstein angle is never in this
likelihood. Map the resulting `dm` posterior through the deflection integral to
obtain a predicted Einstein angle with a credible interval, and report the
fractional error against the catalog image-model angle.

**Direction L -> K (the lens predicts the motions).** For each lens, solve for the
single `dm` that reproduces the catalog Einstein angle exactly. The resolved
kinematics are never in this solve. `beta` is not constrained by lensing, so the
predicted V_rms profile is reported both marginalized over its prior and at its
best value, and scored by chi-square against the measured bins and covariance.

**Consistency statistic.** The physical claim under test is that one state
explains both observables with compatible galaxy properties. The statistic is

    Delta_dm = dm(lensing) - dm(kinematics)

reported in dex and in units of the published stellar-mass uncertainty. A model
that needs a systematically heavier stellar population to explain the lens than
to explain the stars has failed this test, whatever its fitted residuals look
like.

**Control.** The identical two-direction protocol is run with the companion term
removed (`baryons_only`), with the same nuisance freedom and priors, so that any
cross-prediction success or failure can be attributed.

We declare in advance that we will report every system, including the three whose
R10 angles are already known to be accurate, and that we will not drop a lens
whose cross-prediction fails.

## C. Sensitivity before interpretation

The catalog SIE angles are image-model summaries with no released uncertainty.
JR-1's 5 percent was an optimization scale, not a measurement. JR-9 therefore
does not adopt a single lensing uncertainty. Every tension statement is reported
across an explicit scan of assumed fractional angle uncertainty, and the value at
which each discrepancy would fall below 1 sigma is stated.

## D. Radial response kernels and the joint template diagnostic

Before proposing any spatial correction we compute, for each lens, the two
functional derivatives of the observables with respect to companion mass added in
a thin shell at 3D radius r:

    K_lens(r)  = d(theta_E) / d(shell mass at r)
    K_kin(r,i) = d(V_rms in bin i) / d(shell mass at r)

Their ratio decides whether a redistribution that repairs an under-bent lens
without disturbing an already accurate V_rms profile exists at all. If the two
kernels have the same radial support over the measured range, the proposed
"spatial redistribution" direction is not available in these systems and we will
say so.

We then run the declared screening statistic from the starter package on the real
40-bin plus six-angle data vector:

    A_hat = (t' W r) / (t' W t),   sigma_A = 1 / sqrt(t' W t)

with r = observed minus R10 baseline, t = candidate minus baseline, and W the
inverse covariance after marginalizing the local `(dm, beta)` nuisance modes of
every lens. Separate motion and lensing amplitudes are computed as a diagnostic
only; the physical model is not granted two independent amplitudes.

## E. One interaction, not two fitted rate laws

JR-5 through JR-8 chose `k_plus ~ exp(2h)` and `k_minus ~ exp(4h^2)`
independently. JR-9 replaces both with a single matter-assisted event

    gamma + M  <->  gamma' + chi + M'

reduced to two effective channels with coupling `u`, detuning `Delta` and
environmental dephasing `Gamma`, giving

    k = 2 u^2 Gamma / (Gamma^2 + Delta^2)

The deliverable is that `u`, `Delta` and `Gamma` are computed from material
quantities with units, not assigned per galaxy:

* `Gamma(n, T)` from the collision rate of the medium;
* `Delta(q, v) = Delta_0 - q . v` from the advected material response
  `S_v(q, omega) = S_0(q, omega - q . v)`, which is where direction enters;
* `u` from the transition matrix element and the local density.

Forward and return rates come from the same matrix element with their own
occupation factors, so their ratio is derived, not chosen.

We declare the falsifiable consequence in advance: because `dk/dGamma` has the
sign of `Delta^2 - Gamma^2`, one kernel predicts that the orientation response
**changes sign** at a computable critical material condition. Sources on opposite
sides of that condition must respond to gas alignment in opposite directions.
That is a prediction about source physics that does not mention gravity, and the
next resolved gas sample can contradict it.

We also state in advance what this section does not do: it does not derive the
matrix element from a specified companion identity, it does not close the energy
ledger, and it does not produce a gravitational source. A rate law is not a force
law.

## F. What would count as failure

* `Delta_dm` larger than the published stellar-mass uncertainty in a consistent
  direction across systems.
* K -> L and L -> K requiring incompatible galaxy properties.
* `sqrt(t' W t)` below about 1, meaning the proposed correction is undetectable in
  this data and no amplitude can be measured.
* The two response kernels being proportional over the measured range, meaning no
  redistribution can separate the observables.

Any of these is recorded as the result. None of them is grounds for enlarging the
error bars or re-fitting the universal coefficients.
