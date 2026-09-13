# Accumulated force feedback does not repair galaxy predictions

The preceding depth-capture run made progress by testing a common deposition rule on actual rotation curves and exposing its poor fit and neglected feedback. This run tests a distinct force-response law and evolves the deposited density as its gravity grows. It does not test total-potential feedback in the old depth law.

## Formula and physical assumptions

Candidate postulate, not a derived microscopic interaction or originality claim:

    d rho_d(r,s)/ds = C [g_total(r,s)/(g_total(r,s)+g0)]^2
    rho_d(r,0)=0, 0 <= s <= 1

Known spherical Newtonian equations:

    M_d(<r,s)=4 pi integral_0^r rho_d(x,s) x^2 dx
    g_total(r,s)=g_b(r)+G M_d(<r,s)/r^2
    v_pred(r)^2=v_b(r)^2+G M_d(<r,1)/r

s is normalized accumulated exposure, not an imposed cosmic age. C=(kappa_max/c) integral u_c(t)dt remains a fitted, shared bath exposure. The companion energy is captured rather than destroyed. The density is accumulated throughout exposure; it is not assigned from the final force alone. A comparison branch uses g_b instead of g_total in the rate and is fitted separately with the same two free parameters.

The force magnitude is a new environmental assumption for this candidate. A covariant or preferred-environment definition is absent. Capture responsive to the total depth of a well would be a different law. An optically thin isotropic bath and retained, spatially stationary deposits are assumed; depletion, focusing, support/orbits and nonspherical gravity are not derived. Interior spherical force is independent of exterior shells, so predictions need no invented outer edge; this does not prove global finite deposited mass or a realizable formation history.

## Fit and frozen evaluation

All 89 training galaxies share C and g0. The 29 and 31 evaluation galaxies receive no parameter changes. These are exposed historical partitions, not blind confirmation. The fitting objective gives each galaxy equal weight in squared log10 speed residual.

| Group | Empirical baseline RMS km/s | Force capture without feedback | Force capture with feedback |
|---|---:|---:|---:|
| Training | 20.005 | 62.559 | 62.560 |
| Validation | 27.377 | 52.111 | 52.111 |
| Test | 17.200 | 44.465 | 44.464 |

Final numbers use the finer evaluation at frozen fitted parameters. The respective log10 RMS values without/with feedback are 0.175990/0.175991 (training), 0.146535/0.146534 (validation), and 0.127062/0.127062 (test). Neither fits comparably to the empirical baseline. The tiny differences are not evidence favoring a mechanism.

Without feedback C=4.06668e6 Msun/kpc^3; with feedback C=4.06645e6 Msun/kpc^3. Both fits push g0 to the lower bound, 0.01 (km/s)^2/kpc, for every optimization start. This is not a measurement of g0: the optimizer is moving toward saturated capture. In that regime the rate approaches C, so changing gravity hardly changes it. This explains why feedback has little impact in this fitted branch. It does not establish that feedback is unimportant for every companion model.

## Meaning for the next physical model

Uniform exposure with these saturating capture rules tends toward broadly uniform deposits over the measured region. Known gravity then gives an extra speed proportional to radius in the uniform-density limit, not automatically the required variety of rotation-curve shapes. Adding total-force feedback to this rule does not solve that mismatch.

The tested branch should not replace the empirical relation. A next transport model needs a specified spatial source/exposure or redistribution mechanism, with shared parameters and actual transfer evaluation; arbitrary individual galaxy amplitudes would conceal the problem. Total-depth feedback remains untested across this galaxy sample, as do attenuation and orbital redistribution. Photon supply, joint lensing and the time/brightness requirements remain open.

Input hashes match the prior sample. Deposited density remains between zero and C during integration. Doubling radial intervals and exposure steps changes group speed RMS by less than 0.000185 km/s at the fitted parameters, below the declared 0.1 km/s threshold. This validates the reported descriptive score to that tolerance, not parameter uncertainties or observational agreement. Distances, stellar mass/light, spherical geometry and prior data-use caveats remain in force.

Reproduce with `python research_work/results/luminosity-gravity-transfer/force-feedback.py`. force-feedback-results.json retains all starts, coarse/finer scores and parameters; force-feedback-predictions.json records every actual speed and finer-grid prediction for both branches. All six goals remain open.
