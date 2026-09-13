# Stationary capture with occupancy-dependent transport

13 September 2026. Coupled spherical transfer/occupancy diagnostic; reference retained.

**Normalization correction (13 September 2026).** In the density, capacity and rate equations in this report, C denotes the pre-retention amplitude A=2 C0=9.457178483e7 Msun/kpc^3, where C0=4.728589242e7 is the stored fit parameter. C/original C multipliers are unchanged because the factor of two cancels. See [normalization audit](capacity-normalization-report.md).

## Result

Allowing filled sites to stop capturing makes the system more transparent and increases illumination and deposited density relative to the fixed-field local-capacity calculation. This feedback does not repair its rotation discrepancy. With the shared amplitude adjusted on training galaxies, validation/test RMSE is 28.64/26.70 km/s, compared with 28.63/25.19 for the fixed-field local rule and 30.82/21.41 for the identically adjusted reference. The validation advantage over the reference remains, but the test disadvantage grows. No new preferred law is selected.

This experiment solves the specific stationary energy-transfer closure below across all 149 existing galaxy inputs. It does not derive the microscopic storage states, complete momentum balance or gravitational support, measure external source power, or validate a full theory.

## Explicit physical branch

Assume spherical capacity density Cg(r), where g=[1+(r/a)^2]^-2 is the existing hypothetical capture profile, not the measured ordinary-matter density. The boundary illumination is isotropic, with the existing source proxy X used as its dimensionless loading normalization. Let i(r,n) be directional intensity relative to the incident intensity and J its angular average. Keep the exact-third equilibrium response and make the capture opacity proportional to the empty fraction:

\[
f(r)=\eta[XJ(r)],\qquad \eta(y)=\frac{y^{1/3}}{1+y^{1/3}},
\]

\[
\kappa_{\rm cap}(r)=k_0g(r)[1-f(r)],\qquad
\frac{di}{ds}=-\kappa_{\rm cap}i,\qquad i_{\rm incoming}=1,
\]

\[
J(r)=\frac{1}{4\pi}\int \exp\left[-\int_{\rm incoming\ ray}\kappa_{\rm cap}\,ds\right]d\Omega,
\qquad \rho_d(r)=Cg(r)f(r).
\]

The ray attenuation and angular averaging are known transfer mathematics. Empty-site opacity, capacity Cg, the use of X, and the companion interpretation are hypotheses. The one-third response remains the reference mathematical curve, not a newly derived exponent.

Crucially, all stored energy released at equilibrium is assigned to an escaping channel which is not recaptured in this calculation. This is an explicit alternative branch, not a claim that the incoming companions themselves have been shown to become noninteracting. No spectrum, microscopic selection rule or escape lifetime is derived. Re-emission into the incoming interacting channel would require a source term in the transfer equation and is not included here.

## Rate and energy consistency

For the previously stipulated threshold mixture, local capture per capacity per dimensionless time is XJ(1-f); the averaged release equals it at equilibrium. Let E_cap density be Cg c^2 in consistent SI units. Matching the capture opacity to the kinetic loading scale requires

\[
\lambda X=\frac{k_0u_\infty}{Cc},
\]

where u_infinity is the incoming companion energy density. This is a required mapping, not a measurement of u_infinity or lambda. In particular, using different catalogue X values stipulates different boundary loading strengths; this experiment is not the earlier common-background model. Rescaling C at fixed X leaves an undetermined supply/rate normalization to adjust consistently.

The same absorbed power can be computed from complete-ray attenuation or from volume absorption:

\[
P_{\rm abs}=cu_\infty\,2\pi\int_0^\infty b[1-e^{-\tau_{\rm chord}(b)}]db
=cu_\infty\int \kappa_{\rm cap}J\,dV.
\]

The assigned escaping release carries this power, leaving stationary stored energy unchanged. These expressions use the isotropic absorption cross section and do not multiply projected area by a redundant factor of four. Numerical agreement of the independent chord and volume evaluations is within 6.17e-5 relative across the sample. This closes the stated stationary energy ledger, not the physical mechanism for escape. Capture recoil, radiation forces, support and a stationary matter/field stress balance remain to be supplied; spherical cancellation of total vector momentum does not establish local support.

## Why the effect has this sign

For fixed J, occupancy reduces opacity below k0 g. Reduced opacity increases J, and eta is increasing, so the feedback further increases occupancy. The transfer map is monotone. Iterating it from the full-opacity field below and unattenuated J=1 above gives closing numerical brackets, rather than assuming one initial guess is the solution. The lower field is at least the previous fixed full-opacity J; hence the stationary density cannot be smaller under these assumptions. This ordering is a mathematical consequence of this particular closure, not a universal claim about companion interactions.

## Rotation comparison

Integrate the resulting spherical density within each of the 3150 accepted observed radii, and replace the reference extra-gravity component by the new enclosed-mass ratio. Preserve measured baryonic inputs, the reference scale a and k0, and the exact-third response. Scores use equal-galaxy mean squared velocity errors, followed by a square root, in km/s.

| Model | C/original C | Train RMSE | Validation RMSE | Test RMSE |
|---|---:|---:|---:|---:|
| Frozen reference | 1 | 29.03 | 32.49 | 23.59 |
| Frozen occupancy-feedback model | 1 | 64.02 | 74.81 | 76.04 |
| Training-adjusted reference | 0.747109 | 26.87 | 30.82 | 21.41 |
| Training-adjusted fixed-field local model (previous result) | 0.352688 | 28.75 | 28.63 | 25.19 |
| Training-adjusted occupancy-feedback model | 0.317027 | 29.44 | 28.64 | 26.70 |

One common amplitude is re-estimated on the 89 training galaxies in each adjusted branch, with the same 0..4 multiplier bounds and an interior solution. Two other reference constants remain frozen. The no-background reference receives the same adjustment and objective. No per-galaxy parameters are fitted. These are already-exposed partitions with 29 validation and 31 test galaxies; this is a diagnostic comparison, not a blind test or likelihood analysis with observational/baryonic uncertainty. The fitted C value is not adopted into the reference.

## Numerical verification and its limits

Run `python research_work/results/companion-extensions/occupancy-transport.py`. Input/model hashes, per-radius predictions, central fields, power checks and convergence information are saved in `occupancy-transport-results.json`.

The solver compactifies each infinite incoming ray with z=sqrt(1+b^2) tan(alpha) in units of a. It interpolates vacancy on a log(1+r/a) grid, using a sparse positive ray-integration operator. The constant-vacancy operator is checked against the analytic opacity integral; zero opacity yields J=1. Initial grids have 193/385 radial points, 40/80 angles and 80/160 path quadrature points. Two cases require 769 radial points; the last grid comparison changes a predicted speed by at most 0.0571 km/s. These differences are discretization checks, not rigorous error bounds. Far-grid vacancy is held constant outside r/a=1000; the declining opacity makes the ray contribution small, and this approximation is included in the numerical scope rather than a finite physical halo boundary.

Stationary brackets close to less than 1e-9 in J in 4-11 iterations. This brackets the discretized static transfer problem; it is not a proof of time-dependent stability or uniqueness for a more general field theory. Independent integrated absorption checks are reported above. Baseline scores and all accepted observed radii/velocities are verified without altering data.

## What this resolves and what remains

The fixed-field loophole has been tested for empty-site capture plus an escaping release channel. It increases, rather than reduces, the local-capacity model's density and test residual. Altering the transport consistently therefore does not rescue this branch by itself.

Remaining alternatives require different physics: a capacity distribution derived from actual constituents, return/reverse interactions that change occupancy, or a release law whose occupied states and energy budget are specified. These are not equivalent to freely changing an attenuation factor. Supply, equilibrium preparation, support, emission spectra and joint lensing must still be assessed. The result does not close other photon-transfer or storage branches and does not establish a unified model.

Related: [fixed-field local capacity](local-capacity-report.md), [capture/capacity audit](capture-capacity-report.md).
