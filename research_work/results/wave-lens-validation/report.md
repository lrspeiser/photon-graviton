# Frozen wave parameters improve the reused validation sample

The training-selected wave model improves both observables on all seven eligible validation galaxies, without refitting its particle mass, source fraction or stellar masses. Dispersion RMS falls from **92.12 to 38.22 km/s** and Einstein-angle RMS from **0.6437 to 0.2699 arcseconds** relative to the paired ordinary-matter benchmark.

This is a useful conditional validation result, but **not a strong case for the full photon-companion theory**. The sample has prior exposure under other models, substantial residuals remain, and the photometric mass/geometry assumptions and capture physics are unresolved. The final test-role predictions remain unopened.

## What stayed frozen

The selected common parameters are m=10^-24 eV/c^2 and deposited-source/stellar mass ratio f=0.6, chosen on the 32 training systems. Their calibration record and training hashes were verified unchanged after scoring. The validation protocol was written before these wave predictions.

The same pipeline supplies the supported Schrodinger-Poisson source, spherical isotropic Hernquist stellar tracer, 3-arcsecond aperture, 1.5-arcsecond Gaussian seeing and conditional static geometry. The wave density determines both the potential used for stellar motions and the light-bending integral. No separate lensing multiplier or per-galaxy source parameter is introduced. These field and dynamical equations are **known literature applied to a hypothetical companion identity**, not newly derived photon-conversion laws.

Stellar masses are fixed to the same conditional Salpeter population normalization used during training. The known luminosity-scaling conversion is

\[
M_{\rm conditional}=M_{\rm published}
\frac{D^2(1+z)^2}{D_{L,\rm publication}^2}.
\]

Here D is the model's archived static distance, and the two factors of 1+z represent photon energy loss and the stipulated event stretching. This is a **conditional normalization conversion**, not a full stellar-population fit or independent proof of those propagation effects. The publication's expanding-model luminosity distance is reconstructed only to remove its original normalization; it is not adopted as the fictional geometry. An independent integral checks that reconstruction to 1.1e-14 relative precision.

The published mass-to-light ratio, population age, dust, metallicity, IMF and stellar-evolution assumptions are retained. Their uncertainties are not erased by rescaling. They remain a major limitation on absolute comparisons.

## Sample and exposure

Eligibility uses the original validation role, E morphology, available measured dispersion, finite published Salpeter mass and positive updated I-band radius. Seven systems qualify. J0157-0056 is excluded for missing dispersion and J1103+5322 for non-E morphology; neither exclusion uses the new residuals.

These are the same validation-role galaxies already examined in the empirical-profile pilot. They were excluded from the new 32-system wave calibration, but are **reused validation**, not a pristine blind sample. The earlier empirical-profile failure remains archived and is not superseded as a result for that different model.

## Measured and predicted values

The following are central-value predictions from the frozen wave model. Measured dispersion errors are preserved in `predictions.json`; this table is not a claim that all measurement and model uncertainties have been propagated.

| Galaxy | Measured dispersion, km/s | Wave prediction, km/s | SIE angle, arcsec | Wave angle, arcsec |
|---|---:|---:|---:|---:|
| J0029-0055 | 229 | 190.1 | 0.960 | 0.947 |
| J0216-0813 | 333 | 342.7 | 1.160 | 1.286 |
| J0912+0029 | 326 | 336.1 | 1.630 | 1.864 |
| J0946+1006 | 263 | 189.8 | 1.380 | 0.734 |
| J1420+6019 | 205 | 182.5 | 1.040 | 1.035 |
| J1531-0105 | 279 | 258.0 | 1.710 | 1.845 |
| J1627-0053 | 290 | 242.9 | 1.230 | 1.286 |

| Model, same fixed stellar masses | Dispersion RMS | Median predicted/measured dispersion | Angle RMS | Median predicted/SIE angle |
|---|---:|---:|---:|---:|
| Ordinary matter | 92.12 km/s | 0.668 | 0.6437 arcsec | 0.620 |
| Ordinary matter plus supported wave | 38.22 km/s | 0.890 | 0.2699 arcsec | 1.045 |

Every system's absolute error is smaller with the wave model for both measurements. This is a descriptive paired count, not a statistical significance test. Mean wave residuals are -26.11 km/s and -0.0162 arcseconds. The near-zero mean lens residual does not mean every galaxy matches: positive and negative discrepancies partly cancel.

J0946+1006 remains a conspicuous mismatch, with an angle shortfall of about 0.646 arcseconds and dispersion shortfall of about 73 km/s. It is retained in all metrics. We do not remove it, assign it a special source fraction or reinterpret its input data to improve the summary.

The previously declared equal-weight log-error score is 0.04495 for waves versus 0.32263 for the ordinary benchmark. For comparison, the selected wave training score was 0.01467, with dispersion RMS 31.99 km/s and angle RMS 0.1440 arcseconds. Thus validation performance is worse than training, although it still improves substantially over the paired benchmark. The score is a diagnostic objective, not a likelihood with full uncertainties.

The published SIE angle is an image-model summary. The calculation still lacks a full imaging likelihood, non-spherical stellar/orbit modeling, per-exposure seeing and a jointly refitted population mass posterior. No chi-square acceptance, evidence ratio or cosmological model preference is claimed.

## Numerical verification

The seven source equilibria have mass-normalization errors below 6.7e-10 and normalized virial residuals below 1.4e-9. The largest central potential depth divided by c^2 is 1.21e-5, consistent with the weak-field regime used here. These checks do not prove collective stability.

Doubling aperture/radial quadrature for the first system and the extremes in wave parameter changes dispersion by at most 0.000006 km/s. Independent projected-mass and ray-deflection calculations agree within 4.9e-12 relative at the predicted angles. Ordinary-matter dispersion normalization also agrees with the archived fixed-shape mass-scaling identity within 1.5e-12 km/s. The observed residuals are far larger than these numerical differences.

Input hashes, exact eligibility, preserved shared parameters and sample roles are checked in `verify.py`. The runner accepts training or validation only; a validation call with a parameter-grid index is rejected. No test-role prediction is calculated.

## What this means for the full objective

We have progressed from an imposed gravity correction to a supported field source with common parameters that improves two kinds of measurements outside its calibration sample. That warrants continued testing. It does not explain how photon energy creates the massive wave state, how it settles into these profiles, why f is shared, or how redshift and supernova timing follow from the same physics.

Before calling this a strong case, the priorities are to assess the retained mass/geometry/orbital assumptions, test stability and the production/capture connection, and freeze a final test protocol. A further change motivated by these seven residuals is model development and must be labeled as such; it must not reuse these galaxies as fresh validation. The total photon-supply budget remains deferred and unpassed.

Reproduce with `python research_work/results/wave-lens-validation/prepare.py`, `python research_work/results/wave-lens-training/run.py --role validation`, the same run with `--verify-resolution`, and `python research_work/results/wave-lens-validation/verify.py`. Dependencies are NumPy, SciPy and Astropy. Full paired predictions, conditional mass inputs, exclusions and verification are archived here.
