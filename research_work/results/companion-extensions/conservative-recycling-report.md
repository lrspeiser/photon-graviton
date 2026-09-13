# Conservative recycling: a stationary return-channel endpoint

13 September 2026. Exact limiting solution and comparison on the exposed galaxy sample.

## Finding

If captured energy is released isotropically into the same interacting companion channel without energy loss, a stationary isotropic boundary bath admits the exact interior solution i=J=1. Re-emission replaces captured intensity. With the proposed illumination-independent capacity, the deposited density becomes Cg eta(X), losing the reference's radial attenuation factor.

This endpoint worsens the rotation comparison. After an identical training-only adjustment of the shared amplitude, test velocity RMSE is 30.88 km/s versus 26.70 for the previous escaping-release branch and 21.41 for the reference. Validation remains better than the reference but worse than the escaping-release branch. Neither endpoint establishes a better joint fit.

Recycling is not ruled out generally. The calculation constrains conservative, isotropic, frequency-independent re-emission with the current capacity and equilibrium response. It does not calculate retention after the external bath is removed.

## Assumptions and formulas

Keep capacity density Cg, g(r)=[1+(r/a)^2]^-2, exact-third filling f=eta(XJ), and capture opacity kappa=k0 g(1-f). Assign a fraction omega of released energy to isotropic re-emission in the same channel and the remainder to an escaping channel. The stationary frequency-integrated transfer equation is

\[
\frac{di}{ds}=-\kappa i+\omega\kappa J,\qquad
J=\frac{1}{4\pi}\int i\,d\Omega,\qquad i_{\rm incoming}=1.
\]

This is known absorption/re-emission transfer mathematics. Applying it to companions, choosing isotropic re-emission with the same effective interaction law, and the capacity interpretation are hypotheses. The model is gray: it does not derive a spectrum or frequency redistribution. Release into a different spectrum with a different opacity is not covered by this equation.

The previous occupancy-transport calculation is omega=0. This calculation uses omega=1. Substitution gives the exact stationary solution

\[
i=J=1,\qquad f=\eta(X),\qquad \rho_d=Cg\eta(X).
\]

In a finite-optical-depth medium with an isotropic boundary and no independent internal power source, this constant solution is also the usual conservative-transfer equilibrium. The present analysis uses that solution; it does not claim a dynamic stability proof for the coupled storage system. More generally, passive transfer places J at or below the boundary maximum, whereas adding a positive return source can increase illumination above the non-returning case. These field limits do not bound fitted velocity residuals: re-adjusting parameters or changing spatial effects can make intermediate omega behave differently. Intermediate values have not been fitted or selected here.

## Gravity calculation

For R/a=x, the known spherical mass integral is

\[
M_d(<R)=2\pi C\eta(X)a^3\left[\arctan x-\frac{x}{1+x^2}\right],
\qquad M_d(\infty)=\pi^2 C\eta(X)a^3.
\]

The exact-third response and all source proxies remain fixed. Each observed-radius prediction uses the ratio of this mass to the reference attenuated mass:

\[
v_{\rm recycle}^2=v_b^2+
\frac{M_{\rm recycle}(<R)}{M_{\rm ref}(<R)}(v_0^2-v_b^2).
\]

The formula for spherical gravity is established mathematics, not a new relativistic lensing law. Since both profiles contain the same eta(X), it cancels from the enclosed-mass ratio. That does not remove one-third from the original amplitude or imply a derivation of it. Deposits are still interpreted through the hypothetical g profile, not through a directly observed storage capacity.

## Comparison

The sample contains the existing 149 galaxies and 3150 accepted rotation points. Scores are equal-galaxy velocity RMSE in km/s, with 89/29/31 already-exposed training/validation/test partitions. Observational covariance and baryonic uncertainties are not propagated in this diagnostic.

| Branch | C/original C | Train | Validation | Test |
|---|---:|---:|---:|---:|
| Frozen reference | 1 | 29.03 | 32.49 | 23.59 |
| Frozen full recycling | 1 | 72.07 | 81.45 | 94.41 |
| Training-adjusted reference | 0.747109 | 26.87 | 30.82 | 21.41 |
| Training-adjusted escaping release (previous calculation) | 0.317027 | 29.44 | 28.64 | 26.70 |
| Training-adjusted full recycling | 0.277800 | 30.77 | 28.82 | 30.88 |

Each adjusted case re-estimates only C on training velocity MSE using common multiplier bounds 0..4; the fitted solutions are interior. k0 and a/R_d stay frozen. There are no per-galaxy adjusted parameters. The same adjustment is applied to the reference, so the comparison does not credit recycling for a benefit that comes from changing the amplitude objective. None of the adjusted parameters or branches replaces the reference. The existing partitions do not constitute an untouched-system test.

## Energy interpretation

At omega=1, stationary local capture and re-emission have equal power. The active radiation field has no net stationary sink, and the integrated radiation flux entering the object equals the flux leaving it. The store also has no net energy change. This does not create a net luminosity excess above the isotropic bath, although capture/release turnover is nonzero.

There is no independent steady drain into an invisible channel in this branch. Initial formation still requires supplying the stored energy, and removing the boundary illumination is a different, time-dependent problem in which energy can escape to infinity. Consequently, this endpoint does not demonstrate permanent memory, free initial energy, or retention in an isolated galaxy. The actual source history and normalization remain unknown. A uniform stationary bath must itself be reconciled with global propagation and gravity; it is not a supplied cosmological solution.

The returning channel is a frequency-integrated ansatz. Recoil, microscopic transitions, angular redistribution, detailed balance across frequencies, gravitational support and stability have not been derived. Vanishing net radiative energy transfer at equilibrium does not establish these missing requirements.

## What we learned

The two simple release endpoints now have explicit meanings. Escaping release attenuates the bath but costs continuing net converted power; conservative same-channel release restores the undimmed stationary bath and increases occupancy. Neither improves the test score over the matched reference with the current capacity. This narrows the next physical question to the actual capacity and reciprocal transition rules, rather than assuming that an unspecified return flow fixes the profile.

Intermediate return fractions, anisotropic emission, frequency-changing escape and state-dependent reverse reactions are separate choices requiring defined laws. They are not excluded by these endpoint results, and no interpolation is claimed to solve the observed discrepancies. The broader photon-transfer, lensing, source and support goals remain open.

## Verification

Run `python research_work/results/companion-extensions/conservative-recycling.py`. It verifies original input hashes and accepted data arrays, records model-input hashes, reproduces baseline scores, checks the analytic unattenuated mass against numerical radial quadrature, and compares reference-mass angular/radial quadrature orders 64 and 128. The largest corresponding speed difference is 1.55e-9 km/s. The constant transfer solution and its energy cancellation are analytic statements; numerical mass checks do not independently validate a physical companion channel. Results are in `conservative-recycling-results.json`.

Related: [escaping-release transport](occupancy-transport-report.md), [fixed-field capacity](local-capacity-report.md).
