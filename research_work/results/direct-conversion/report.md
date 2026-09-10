# Direct conversion: integrated branch and redshift benchmark

9 September 2026. [Active branch specification](../../../research_plan/direct-conversion-bound-companions.md), [fixed protocol](protocol.md), [code](run.py), [all predictions](all_predictions.csv), [results](results.json).

**Post-run user clarification:** comparable performance to expansion is sufficient; superiority is not required. The fixed protocol and numerical outputs are preserved. Their nearly equal redshift errors are encouraging for this limited comparison but do not prove statistical equivalence or parity across the full set of observations. See [prior art and the revised success criterion](prior-art-and-parity.md).

## Outcome

Adopt direct continuous photon-to-companion conversion as the working branch, without an added clock factor or environmental suppression of conversion. This removes the earlier clock-sign conflict from the starting propagation law. Retain light-speed companion transport, no void storage, and the proposed capture into extended bound configurations in wells. The capture mechanism and stable configurations are not established by choosing this branch.

The branch predicts a definite conversion redshift at 100 million light-years, but the recovered galaxy comparison does not show it outperforming the selected expansion model. Its stationary propagation closure also lacks whole-event time stretching. These are substantive findings, not reasons to conceal the results or keep adjusting individual objects until they match.

## Prediction at 100 million light-years

**Provenance: established fractional-loss mathematics under the proposed constant transfer law, with an empirically calibrated rate from an already exposed dataset.**

\[
z_{\rm transfer}=\exp(\alpha_0D)-1,\qquad
\alpha_0=0.0002488993286382367\ {\rm Mpc}^{-1}.
\]

At a stipulated geometric path of 100 million light-years, D=30.66013938 Mpc:

| Quantity | Conditional prediction |
|---|---:|
| Conversion redshift | 0.0076604806 |
| Wavelength multiplier | 1.0076604806 |
| Received wavelength for emitted 500 nm | 503.8302403 nm |
| Photon energy remaining | 99.2397756% |
| Energy transferred into traveling companions before capture | 0.7602244% |
| cz, used only as a reporting unit | 2296.5543 km/s |
| Full-event duration factor in stationary fixed-speed closure | 1 |

The cz value does not mean the source must recede at that speed in our model. The redshift is a propagation contribution, not an exact prediction for an arbitrary star or galaxy irrespective of its motion, host potential and measurement frame. The transmitted energy fraction is a per-packet identity, not a calculation of whether the universe supplies enough deposited energy.

## Actual recovered measurements near that distance

Select the five old test-partition group representatives closest in adopted distance to 100 million light-years, without using residuals to choose them. All historical partitions have already been exposed. PGC numbers identify the catalog representatives; no unverified common names are attached.

| Representative PGC | Group PGC | Distance, million light-years | Observed CMB-frame z | Direct-conversion z | Expansion-comoving benchmark z |
|---|---|---:|---:|---:|---:|
| 19476 | 19476 | 106.225 | 0.00915300 | 0.00813926 | 0.00815432 |
| 46115 | 46247 | 106.715 | 0.01033382 | 0.00817699 | 0.00819203 |
| 6983 | 6983 | 89.997 | 0.00462653 | 0.00689153 | 0.00690662 |
| 51787 | 51787 | 87.948 | 0.00637775 | 0.00673414 | 0.00674916 |
| 20047 | 20047 | 85.551 | 0.00785877 | 0.00655001 | 0.00656496 |

The two predictions are very close to one another, while the object residuals can be much larger. Direct-conversion residuals in cz units are -303.9, -646.6, +679.0, +106.8 and -392.4 km/s respectively. We have not assigned per-object peculiar velocities to remove these residuals. Independent motion, source-calibration and measurement information is necessary to explain scatter without fitting away failures. These examples do not constitute a fresh test or a claim to explain every individual galaxy exactly.

## Explicit expansion comparison

The comparison is spatially flat FLRW with fixed Omega_m=0.3, Omega_Lambda=0.7, no radiation term, and H0 as its one fitted scale. These parameters define a reference model, not an imported active cosmology or an exhaustive search over expansion theories. Conversion also has one scale parameter. Scales are fit on the old training partition for the old test result; each coarse-sky held-out region gets predictions from fits using the other seven regions.

**Provenance: established FLRW distance formulas, used solely as competing predictions.**

\[
E(z)=\sqrt{0.3(1+z)^3+0.7},\quad
D_C=\frac{c}{H_0}\int_0^z\frac{dz'}{E(z')},\quad
D_L=(1+z)D_C.
\]

[Hogg's distance definitions](https://arxiv.org/abs/astro-ph/9905116) distinguish these distances. In the geometric convention we compare the supplied number as static path length in conversion and present-day comoving distance in FLRW; those are alternative model interpretations, not the same photon travel length in expansion. In the luminosity convention both models interpret the supplied number as a luminosity distance. No catalog number is altered. Surface-brightness-fluctuation calibration is more specific than a bolometric standard candle and must be forward-modeled before treating the luminosity sensitivity as a faithful reinterpretation of the data.

For the static no-event-stretch model, standard bolometric accounting gives D_L=r*sqrt(1+z). This implies z=expm1(2*LambertW(alpha*D_L/2)). **Provenance: conditional algebra using the standard Lambert W function, not a new physical law.** It is only an interpretation sensitivity, not a replacement of the active geometric-path law.

| Distance interpretation | Model | Old test RMS, 25 groups | Coarse-sky held-out RMS, all 164 groups |
|---|---|---:|---:|
| Geometric | Direct conversion | 415.414 km/s | 450.716 km/s |
| Geometric | Expansion, comoving distance | 414.991 km/s | 450.502 km/s |
| Luminosity sensitivity | Conversion with its flux relation | 414.660 km/s | 450.358 km/s |
| Luminosity sensitivity | Expansion, luminosity distance | 413.684 km/s | 450.074 km/s |

The independent linear small-distance control gives old test RMS 414.659 km/s and coarse-sky RMS 450.358 km/s. It is not labeled a complete expansion model. Historical training-fit scale values in km/s/Mpc are 74.61814 for geometric conversion, 74.92205 for expansion-comoving, 75.17573 for luminosity-conversion, and 76.03761 for expansion-luminosity.

Paired bootstrap percentile intervals for conversion RMS minus expansion RMS, resampling the eight held-out sky regions 2000 times with fixed predictions:

- Geometric convention: -0.952 to +1.389 km/s.
- Luminosity convention: -2.983 to +3.583 km/s.

Both span zero. The numerical advantage in this run is slightly on the expansion side, but these data do not meaningfully distinguish the forms at these distances. The intervals are descriptive, conditional on exposed-data predictions and the chosen regions; they are not full theory uncertainties or proof of statistical equivalence. Repeated use of the dataset prevents a fresh-validation claim. We did not compute a Bayes factor or claim to rule out either framework.

## Can the whole system hold together?

| Link | What presently works or is specified | What remains to be established |
|---|---|---|
| Photon energy decreases with distance | A one-parameter empirical trend and definite prediction at 100 million light-years | Microscopic interaction, achromatic frequency shift, phase, linewidth, momentum and symmetry/frame behavior |
| Lost energy enters companions | Scalar transfer terms cancel in a common ledger | Full field/recoil/gravity conservation and production process |
| Companions cross voids without storage | Explicit free-transport/no-capture assumption | A real excitation with the stated speed, lifetime and interactions |
| Capture in deeper wells | Possible effective rate and target stored state | Capture kinematics, binding and response of the well |
| Extended bound patterns | A representation of energy spread among spatial modes | Derived mode spectrum, occupations, spatial support and longevity |
| Extra gravity | A conditional source equation can be written | One response that predicts both rotation and lensing from the generated distribution |
| Source-energy sufficiency | Deferred as requested | Causal formation history and sufficient supply under the selected gravity law |
| Whole-event timing | The minimal closure makes a clear prediction: duration factor 1 | Reproduce observed transient time stretching with the same interaction or record failure |
| Brightness and distances | Bolometric loss predicts a flux factor | Actual bandpass, SBF calibration, extinction and a common distance interpretation |
| Comparable to expansion | Very similar redshift errors in the stated recovered-data comparisons | Formal comparability on independent observations and across the claimed scope remains to be established; superiority is not required |

The timing issue is especially immediate. In this stationary fixed-speed closure, emissions separated by ten days arrive separated by ten days even though each photon has lost energy. At 100 million light-years the difference from a hypothetical factor 1+z is less than one percent, but at larger redshift it grows. Measurements of supernova light-curve stretching provide an empirical target, not merely an expansion assumption to discard. The [DES study of 1504 supernovae](https://arxiv.org/abs/2406.05050) reports such a relation; this run does not reanalyze those observations or their selection/systematics. No assertion that expansion is proven uniquely follows, but describing it as having no substantial observational support would be inaccurate.

## What to do next

Keep the constant direct-conversion rule as the simple reference branch. Do not add adjustable time fields or per-galaxy rates solely to improve its fit. The next physical task is a common photon-companion interaction that predicts both carrier energy and signal timing while preserving narrow lines and images. In parallel with that modeling, pursue the already identified independent distance/velocity provenance and fresh-data uncertainty requirements. Bound-mode capture and joint lensing/rotation remain necessary downstream goals; the source budget stays deferred until revisited explicitly.

## Verification

The run reproduced the prior conversion scale and old test RMS, preserved all 164 distances and spectral shifts, and verified every expansion training-fit prediction with independent adaptive integration plus root bracketing (328 checks, maximum absolute z difference 2.41e-15). The luminosity-conversion inverse reconstructs its input distance to 4.45e-16 relative precision. All five models have predictions for every coarse held-out region. No new astronomical observations, calibrated physical interaction, or complete model validation were produced.
