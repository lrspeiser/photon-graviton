# First test result: improvement over fixed masses, no robust advantage over a mass correction

The frozen wave model improves the reserved test galaxies substantially relative to the original photometric-mass ordinary model. Against the stronger common-stellar-mass benchmark, however, the result is mixed: **waves have lower dispersion RMS but higher lensing RMS**, and only a small advantage in the predeclared combined score.

That small score advantage depends strongly on one system, J0935-0003. This test therefore does **not** establish robust predictive superiority over the mass-corrected ordinary model, much less a strong case for photon-fed companions.

## Protocol and first test exposure

The eligibility, models, common parameters, assumptions and scoring rule were committed and pushed before preparing test inputs or calculating their predictions: [protocol freeze at commit 97c7afd](https://github.com/lrspeiser/photon-graviton/blob/97c7afd/research_work/results/wave-lens-test/protocol.json). The protocol hashes the calibration records, source tables, geometry, source solver and prediction/preparation code. Those hashes pass verification after execution.

This is the first model-prediction scoring of the original test-role systems. Public catalog rows and some aggregate sample information had previously been accessed, so it is not pristine data blinding. No parameter, eligibility rule or stellar mass was changed in response to these test residuals.

Twelve of the 13 test-role systems meet the frozen input requirements. J0405-0455 lacks a Salpeter population mass and is excluded for that reason alone. All eligible systems completed the numerical calculation and remain in every primary metric. The existing training and validation roles are unchanged.

**The test residuals are now exposed.** These galaxies cannot serve as a fresh independent test for subsequent changes motivated by these results. A new independent sample or explicitly redesigned evaluation will be needed for further development.

## Frozen models and assumptions

The three models are:

1. Ordinary matter at the conditional photometric stellar masses.
2. The same ordinary profile multiplied by the training-fitted common factor lambda=1.6453621586.
3. Ordinary matter at the original conditional masses plus the supported stationary wave with m=10^-24 eV/c^2 and source fraction f=0.6, selected on training.

All retain isotropic luminous stellar orbits, Hernquist light/ordinary-mass shapes, the same rescaled Salpeter population estimates, static redshift-derived geometry, three-arcsecond aperture and 1.5-arcsecond Gaussian seeing. The source density predicts both stellar motion and lensing; there is no independent lensing factor.

The Schrodinger-Poisson, Jeans and weak-lensing equations are **known physics used for a hypothetical companion identity**. The common f and lambda are calibrated assumptions, not first-principles energy-transfer laws. The published mass rescaling is a conditional luminosity-normalization adjustment, not a full population refit that removes age/IMF/dust/stellar-evolution priors. No expanding-universe distance is adopted as the fictional geometry; the publication's luminosity distance is reconstructed only to undo its original mass normalization.

## Predeclared primary result

The same score used in calibration averages squared logarithmic fractional errors in dispersion and angle with equal weight. Lower is better. It is a diagnostic score, not an uncertainty-weighted likelihood, evidence ratio or significance statistic.

| Test model | Dispersion RMS, km/s | Lens-angle RMS, arcsec | Joint score |
|---|---:|---:|---:|
| Original ordinary masses | 78.31 | 0.4473 | 0.19085 |
| Common stellar-mass correction | 50.27 | **0.2381** | 0.04226 |
| Frozen wave | **36.14** | 0.2839 | **0.03828** |

The primary wave-minus-mass-correction score difference is -0.003989. Although waves have the lower combined score, they do not improve both observables over this benchmark. Relative to the mass-corrected model, waves have smaller absolute dispersion error for 6 of 12 galaxies and smaller absolute lens error for 5 of 12. This contrasts with the more uniform improvement seen in the reused validation sample.

For waves, median predicted/measured ratios are 0.915 for dispersion and 0.912 for lens angle. Mean residuals are -20.44 km/s and -0.0334 arcseconds. The mass-corrected benchmark has median ratios 1.004 and 1.055, respectively. Medians alone hide the substantial scatter and influential systems.

## Individual wave predictions

| Galaxy | Measured dispersion, km/s | Wave dispersion, km/s | SIE angle, arcsec | Wave angle, arcsec |
|---|---:|---:|---:|---:|
| J0252+0039 | 164 | 168.7 | 1.040 | 0.706 |
| J0737+3216 | 338 | 327.8 | 1.000 | 1.413 |
| J0935-0003 | 396 | 307.9 | 0.870 | 0.825 |
| J0959+0410 | 197 | 160.7 | 0.990 | 0.713 |
| J1016+3859 | 247 | 208.2 | 1.090 | 0.953 |
| J1023+4230 | 242 | 211.7 | 1.410 | 1.169 |
| J1134+6027 | 239 | 201.6 | 1.100 | 0.973 |
| J1250+0523 | 252 | 293.0 | 1.130 | 1.820 |
| J1443+0304 | 209 | 192.9 | 0.810 | 0.796 |
| J2238-0754 | 198 | 200.4 | 1.270 | 1.151 |
| J2303+1422 | 255 | 241.9 | 1.620 | 1.487 |
| J2321-0939 | 249 | 225.8 | 1.600 | 1.524 |

The source tables' dispersion errors and all three models' predictions are retained in `predictions.json`. SIE angles are image-model summaries, not raw image fits. No claim is made that these differences are within the full observational uncertainty; that likelihood has not been constructed.

## Influence diagnostic, with no sample removal

After seeing the primary result, we calculated the combined-score difference while omitting each system in turn solely to measure influence. This was a **post-test diagnostic**, not a predeclared model-selection rule. All 12 systems remain in the official score above.

Omitting J0935-0003 in that diagnostic changes the wave-minus-mass score from -0.003989 to +0.01590, reversing the ranking. Its individual score is 0.03307 for waves versus 0.25585 for the common mass correction. Conversely, waves perform poorly for other systems, including the substantial lens-angle overprediction at J1250+0523. The small full-sample advantage should therefore not be characterized as robust across galaxies.

This diagnostic is not a reason to discard J0935-0003, alter its observations or declare its measured dispersion unreliable. Any future investigation of influential objects must preserve the current test result and use independently justified data checks.

## Verification

The frozen input hashes remain unchanged. All three models use exactly the same 12 eligible galaxies. The common-mass dispersion obeys the independent sqrt(lambda) scaling identity to 2.3e-16 relative. Projected-density and ray-deflection calculations agree within 6.8e-12 relative at the predicted lens angles.

Doubling aperture/radial quadrature for the first system and the extremes in dimensionless wave parameter changes dispersion by at most 0.000043 km/s. Source normalization errors are below 4.2e-10 and normalized virial residuals below 3.0e-7. The maximum central potential depth divided by c^2 is 1.19e-5. These checks support the numerical weak-field calculation; they do not prove collective stability or physical formation of the wave source.

The differences between models and observations are much larger than these numerical discrepancies. Adjusting integration precision cannot turn the mixed test result into universal agreement.

## Big-picture consequence

This closes the first reserved-sample evaluation of the current conditional galaxy/lens pilot. It does not complete the full research objective. The supported wave profile is worth understanding, but the new test weakens a claim that it explains the data distinctly better than a coherent ordinary-mass correction. The wave model also had two calibrated parameters versus one in that benchmark, and neither comparison includes a complete population/orbit/geometry likelihood.

The next work should address that ambiguity and the missing physical connection, rather than continually retuning the now-exposed test galaxies. A common treatment of stellar-population and orbital uncertainties, independent systems and a specified photon-conversion/capture action are needed. Supernova event timing, broadband propagation, lossless companion transport, source accumulation, stability and the deferred total photon-energy supply remain unresolved. A successful descriptive gravity fit would not by itself establish any of them.

Reproduce using `prepare.py --role test` in the wave-lens-validation directory, `run.py --role test` and its `--verify-resolution` variant in wave-lens-training, then `python research_work/results/wave-lens-test/verify.py`. The test protocol rejects a parameter grid or orbital audit on this role. Inputs, all predictions, exclusions and the post-test influence diagnostic are archived here.
