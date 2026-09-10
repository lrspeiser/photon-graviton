# Much of the wave improvement is degenerate with stellar mass normalization

A single common multiplier of **1.64536** on the conditional photometric stellar masses explains a large part of the improvement previously attributed to the added wave source. The frozen wave setting still has smaller validation errors in both observables, but its advantage is much less pronounced against this benchmark than against unadjusted stellar masses.

This comparison does **not** establish that increasing the stellar masses by 64.5% is physically allowed. It isolates an important ambiguity: motion and lensing respond to gravitational mass, while the ordinary mass supplied by stellar-population modeling remains uncertain. The added-source interpretation needs to be assessed alongside that uncertainty.

## Declared benchmark and fit

Before this benchmark run, `protocol.json` specified one shared positive multiplier lambda in [0.5,3], fitted only on the 32 training systems using the same equal-weight squared-log residual objective as the wave grid. It rescales the entire ordinary Hernquist mass profile; light profile, distance, aperture, seeing and isotropic stellar-orbit assumptions stay fixed. There is no companion source in this benchmark.

The relation

\[
M_{b,i}=\lambda M_{\mathrm{population},i}
\]

is an **empirical nuisance-parameter model**, not a new gravity law, a derived population correction or an accepted IMF change. Every galaxy shares the same lambda. No individual mass is fitted to its own data.

For a fixed ordinary-matter shape, the **known Jeans scaling** gives

\[
\sigma_{\mathrm{ap},i}(\lambda)=\sqrt\lambda\,
\sigma_{\mathrm{ap},i}(1).
\]

The bending at a fixed impact parameter scales linearly with lambda, but the Einstein angle generally does not:

\[
\theta_E=\lambda\,\alpha_1(\theta_E).
\]

Here alpha_1 includes the same source/lens distance ratio as in the existing pipeline. The root is recomputed at every trial lambda; multiplying the old Einstein angle by lambda would be incorrect for an extended galaxy.

For efficiency, the known projected Hernquist mass is used, with direct quadrature near its removable singularity at r/a=1. It is independently checked against the line-of-sight integral over 100 impact parameters. These are established gravitational relations, not a new photon-companion formula.

The fitted lambda is 1.6453621586, an increase of approximately 0.216 dex. The fit is inside the declared interval. A 121-point scan checks the bounded optimizer's result. The calibration is written to disk before loading the validation prediction table; lambda is not changed afterward. The validation sample has prior exposure, and this new benchmark is a follow-up prompted by those earlier analyses, so it is not pristine blinding.

## Paired results

Each row below uses the same training or validation objects and measured quantities. The fixed wave model retains its previously selected m=10^-24 eV/c^2 and f=0.6; it is not refitted for this comparison.

| Sample/model | Dispersion RMS, km/s | Lens-angle RMS, arcsec | Joint log-error score |
|---|---:|---:|---:|
| Training: original ordinary masses | 62.77 | 0.4628 | 0.14085 |
| Training: shared stellar-mass multiplier | **31.04** | 0.2278 | 0.02430 |
| Training: frozen wave | 31.99 | **0.1440** | **0.01467** |
| Validation: original ordinary masses | 92.12 | 0.6437 | 0.32263 |
| Validation: frozen stellar-mass multiplier | 45.85 | 0.3541 | 0.06287 |
| Validation: frozen wave | **38.22** | **0.2699** | **0.04495** |

On training, the common mass correction is slightly better for dispersion, while the wave profile is better for lensing and the combined score. On validation the wave model is better on both RMS measures. Consequently the wave advantage is not entirely removed by this simple mass normalization, but much of the gain over the original low-mass benchmark is not unique to a wave source.

The wave model was selected from a two-parameter grid, while this mass benchmark has one fitted parameter. Neither comparison is a complete uncertainty-weighted likelihood. The score differences are not a statistical model preference, detection significance or evidence ratio. They cannot establish that a 0.216-dex mass shift is ruled out or allowed by population data.

The model-predicted median ratios in validation are 0.856 for dispersion and 0.959 for lens angle under the common mass correction, compared with 0.890 and 1.045 for waves. A median close to one does not resolve the per-object scatter. No object is removed from either summary.

## Independent checks

The projected-mass formula agrees with direct quadrature to better than 3.8e-14 relative over the checked impact range. Setting lambda=1 reproduces the original ordinary-model angles within 4.7e-13 relative. The optimizer's score is below the best point in the 121-point scan, and its actual evaluations are retained. The forward predictions use fixed-shape mass scaling and a new lens root, rather than approximating both observables by one common multiplicative change.

The script also records, for diagnosis only, the individual lambda that would match each galaxy's dispersion and the one that would match its lens angle. Those values are **not adopted as per-galaxy fitted masses**. A mismatch between them indicates that a scalar mass correction cannot fit both observables for that object within the assumed shape and orbits.

`frozen-calibration.json` records the common training factor and input hashes. `predictions.json` contains 39 paired-object predictions across the 32 training and seven reused-validation systems. The test-role predictions remain unopened.

## Consequence for the full theory

The evidence currently supports a narrower statement than “the data identify photon-fed companions”: a supported extra-source profile with shared parameters can fit these conditional mass-and-geometry inputs better than either the unadjusted ordinary profile or this one-parameter mass-rescaling benchmark. Independent stellar population and orbit information is needed to determine whether that advantage survives plausible ordinary-matter models.

The useful next step is a common treatment of mass and orbital uncertainties across both hypotheses, followed by a frozen final-test comparison. Allowing flexibility only in the wave model would overstate its advantage; assuming every mass shift is allowed would also be unjustified. The reported population errors do not automatically encompass a coherent IMF, age or luminosity-calibration systematic, so those cannot be treated as interchangeable without modeling them.

The photon conversion, supernova timing, companion capture/transport, source normalization, stability and total supply requirements remain unresolved. The shared wave fraction f and this stellar multiplier are alternative gravitational descriptions; neither derives an energy-conversion mechanism. The deferred photon budget remains unpassed.

Reproduce with `python research_work/results/stellar-mass-degeneracy/run.py`. Requires NumPy and SciPy. The full scan, calibration, per-object diagnostics and both samples' scores are archived alongside this report.
