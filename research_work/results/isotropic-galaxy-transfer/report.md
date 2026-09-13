# All-direction interception improves capture transfer, but is insufficient

The preceding force-feedback run produced a poor, nearly saturated deposition fit. This experiment adds the user's proposed interception effect: incoming companions can be captured before reaching the center. The spatial deposition is computed from attenuated rays and a specified capture coefficient rather than assigned directly as an extra acceleration.

## Actual rotation predictions

Fit on the historical 89 training galaxies, then freeze every common parameter for 29 validation and 31 test galaxies. These partitions are already exposed, not blind confirmation. Each galaxy receives equal weight in squared log10 speed residual. The table uses finer-grid predictions at frozen parameters.

| Group | Empirical baseline RMS km/s | Transparent capture RMS km/s | Intercepted capture RMS km/s |
|---|---:|---:|---:|
| Training, 89 galaxies | 20.005 | 37.406 | 30.709 |
| Validation, 29 galaxies | 27.377 | 35.137 | 33.184 |
| Test, 31 galaxies | 17.200 | 42.987 | 23.994 |

The respective transparent/intercepted log10 RMS values are 0.151937/0.141004 (training), 0.123156/0.119980 (validation), and 0.107575/0.098545 (test). Both evaluation groups improve in both metrics when interception is included. Interception adds one fitted global parameter, and these are repeated exploratory comparisons, so no statistical significance or unique physical cause is established. The empirical acceleration relation still fits better.

## Formula and fitted quantities

The **postulated** capture coefficient is kappa=k0/[1+(r/a)^2]^2, with a=s R_disk using measured stellar scale length. **Known absorption mathematics** gives the intensity surviving along each ray as exp(-tau), with tau the integral of kappa from the incident boundary. The angular mean J produces rho_d=C J/[1+(r/a)^2]^2. **Known spherical gravity** integrates this density to predict rotation. The full formula and its energy bookkeeping are in [derivation.md](derivation.md). These are not claimed as new mathematical forms or a first-principles graviton interaction.

| Parameter | Transparent approximation | Interception model |
|---|---:|---:|
| C, Msun/kpc^3 | 1.35474e7 | 2.71737e7 |
| a/R_disk | 4.02214 | 4.03298 |
| k0, kpc^-1 | Not independently fitted | 0.122551 |

All starts converge away from bounds. The two branches are separately fitted; the observed improvement is the effect of allowing interception plus refitting shared parameters, not the effect of attenuation at fixed C. C=(k0/c) integral u_external dt remains fitted exposure, not a completed stellar energy-supply estimate or assumed cosmic age. The transparent control neglects depletion at fixed accumulated exposure; it does not claim nonzero capture at exactly zero opacity.

## What this changes

Accounting for companions intercepted on the way inward is useful within this specified candidate. It helps prevent every inner region from receiving the full undepleted external supply. This is stronger than a synthetic demonstration: the common fitted rule transfers better to both historical evaluation groups.

However, the opacity shape is still prescribed using stellar size. It is not derived from binding, and the model assumes deposits remain in place. Support/orbital redistribution, opacity feedback, nonspherical geometry, actual external stellar supply and collision behavior remain unresolved. Source amplitude is free; the redshift conversion rate does not independently set it. This therefore cannot establish photon origin or satisfy motion and lensing jointly.

Retain this as the current tested transport candidate for further physical comparison, while retaining the empirical formula as the better-fitting target. The next useful comparison is whether the same deposited spatial source predicts an independent gravitational observable, or whether a physically specified redistribution law improves its rotation prediction. Do not select an arbitrary per-galaxy profile from residuals. All six goals remain open.

## Reproduction and limits

`python research_work/results/isotropic-galaxy-transfer/run.py` writes parameters, all optimizer starts and group scores to results.json, and all 3,150 observed/predicted speeds per branch to predictions.json. Doubling radial/angular resolution changes group RMS by less than 0.00094 km/s at the fitted parameters, meeting the declared 0.1 km/s criterion. This verifies the prediction calculation, not observational uncertainties. Input hashes match the earlier sample; ordinary-matter templates, fixed mass/light, stipulated distances, spherical approximation and exposed-data caveats remain unchanged.
