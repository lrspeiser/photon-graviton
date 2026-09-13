# Same deposits predict lensing and motion: no joint success yet

The preceding Milky Way run made concrete progress through frozen capture transfer. Here the galaxy-trained capture source is applied to six exposed SLACS systems and used for both resolved stellar dynamics and lensing. Inner stellar bins set ordinary stellar mass and constant orbital anisotropy only; outer stellar bins and lens angles do not enter the fit.

## Conditional morphology extension

These early-type galaxies do not supply the exponential disk scale used in SPARC. We explicitly postulate an equivalent scale R_equiv=R_half_light/1.67834699 using the archived combined equal-area half-light radius. The divisor is a known exponential-disk identity; applying it to these systems is a new morphology assumption, not a measured exponential disk. Measured Sersic component profiles remain the ordinary-mass and stellar-tracer shapes.

Capture parameters C,k0,s are unchanged from SPARC. The density depends on observed size and the fixed bath model, not the fitted stellar mass. The extra stellar-dynamical term consequently carries no empirical mass^p factor. Fits preserve the prior stellar-mass/anisotropy bounds, inner covariance, PSF and four starts. The source geometry uses the archived conditional nonexpanding distances, not independent lens distance measurements.

The lens prediction uses known equal-potential weak-field bending of the ordinary plus deposited mass:

    deflection(b)=4G/c^2 integral_0^(pi/2) M_total(<b/cos t)/(b/cos t) dt,
    theta=(D_ls/D_s) deflection(D_l theta).

The capture profile is a project postulate using known absorption mathematics. No new field equation, microscopic capture derivation or literature-priority claim follows. A different stress/metric response would require another explicit physical theory; none is fitted to these lensing outcomes.

## Actual comparison

| Galaxy | Catalog SIE angle, arcsec | Interception prediction, arcsec |
|---|---:|---:|
| J0037-0942 | 1.530 | 1.552 |
| J1112+0826 | 1.490 | 1.210 |
| J1204+0358 | 1.310 | 1.333 |
| J1402+6321 | 1.350 | 1.408 |
| J1621+3931 | 1.290 | 1.150 |
| J1630+4520 | 1.780 | 1.415 |

Catalog SIE radii are model-derived summaries of lens imaging, not exact spherical Einstein rings. No lensing error distribution is fitted or inferred from this table. The last, second and fifth systems remain substantially underpredicted; the other three are close or modestly overpredicted.

| Diagnostic, six systems | Previous empirical extra source | Transparent capture | Intercepted capture |
|---|---:|---:|---:|
| Inner stellar chi-square sum | 41.490 | 38.012 | 47.614 |
| Outer conditional residual-square sum | 49.394 | 39.135 | 57.978 |
| Lens fractional RMS | 12.036% | 12.454% | 12.352% |
| Anisotropy-boundary fits | 0 | 0 | 0 |

Outer predictions use the measured inner/outer covariance conditional on fitted inner residuals, as in the previous audit. These are plug-in diagnostics, not fully marginalized posterior predictive scores. A higher outer residual-square sum is worse under this fixed calculation; it is not a global exclusion probability. The prior empirical source is a descriptive reference with a different functional dependence, not the same morphology rule.

Interception yields a tiny lens RMS change relative to transparent capture while worsening both inner fits and outer stellar predictions. Its successful spiral-rotation comparison therefore does not establish joint success on these early-type systems. No capture amplitude or lensing multiplier was adjusted to hide the shortfalls. All systems and successful optimizer counts are retained in lensing-results.json.

## Consequence

The current capture/interception candidate has now faced three kinds of gravitational comparison: spiral rotation, Milky Way rotation/conditional vertical force, and early-type stellar motion/lensing. It is useful as a specified transport model, but has not matched them together. The morphology extension, spherical source, ordinary population mass assumptions, stable retention and common external exposure remain material limitations. This result does not single out which missing physics is responsible, and it does not reject every possible companion mechanism.

Before claiming an improvement, any revised model must predict a physical spatial distribution and be tested against motion and lensing together, with geometry and source assumptions stated in advance. Redshift/timing/brightness, actual stellar supply, Milky Way bulge structure, collisions and an untouched distinctive prediction remain open. All six goals remain open.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/lensing.py`. The script reuses the established deprojection, annular Jeans and covariance calculations, while replacing the extra source consistently in dynamics and bending. lensing-results.json includes all stellar predictions, fitted nuisance parameters, lens predictions, boundaries and the capture-input hash.
