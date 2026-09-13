# Hsiao source-template identity and assumption audit

The cached Hsiao07.dat is byte-identical to snflux_1a.dat in the author-hosted release. This resolves the source-version ambiguity for a conditional joint spectral/photometric pilot. It does not make the template an independent physical clock or a first-principles luminosity prediction.

## Verified data properties

- 254,506 rows: 106 daily phases from -20 to +85 days relative to B maximum, each with 2,401 wavelengths from 1000 to 25000 Angstrom at 10-Angstrom spacing.
- All values are finite and the tabulated spectral flux values are nonnegative.
- The first phase has almost wavelength-independent values near 1.3e-38. This numerical near-zero boundary must not be interpreted as a measured early spectrum or extrapolated explosion law.
- The accompanying UBVRIYJHK light-curve table has 111 phases from -20 to +90 days. Its final five days have no corresponding spectral grid. An initial grid-equality check identified this difference; the audit now preserves the distinct ranges.
- Cached spectrum and author member have SHA-256 6bd032004eccc7b5b52f7c021f0c4b541d801a70195028bf5af0a2578f554d78. The original local dependency manifest is checked as well.

## Scientific role

The source paper describes a composite assembled from observed spectra, with weights for coverage, source multiplicity, type and stretch. Its final spectral series is adjusted to a light-curve template. It also discusses ambiguities relating spectral evolution to intrinsic light-curve stretch. Consequently its temporal and color behavior are empirical inputs, not consequences of companion propagation. See [Hsiao et al. 2007, sections 5.1-5.4](https://arxiv.org/pdf/astro-ph/0703529).

The [author data page](https://astrophysics.physics.fsu.edu/~hsiao/data/) identifies the phase and wavelength coverage and the October 2007 update. It does not turn the file's normalization into an independently known luminosity for SN 2006mk.

This source can provide a declared reference history L_lambda(phase, wavelength) for an exploratory forward calculation. A successful fit would mean compatibility conditional on that history. It would not establish a microscopic interaction, prove that the history is universal, or constitute a blind test of an independently calibrated clock. Exact per-training-object timing conventions and their uncertainty are not resolved by comparing the files.

## Rules for the joint pilot

1. Interpolate only within the spectral phase and wavelength support. Keep pre-event/post-event photometric measurements for a separately specified baseline; do not silently drop them or extend the spectral file to +90 days.
2. Keep the propagation arrival-time mapping separate from any intrinsic source-width parameter. If both merely rescale the same phase axis, only their product is identifiable without additional information.
3. Predict all epochs and bands from the same source history and shared source parameters. Do not fit a separate luminosity to each data point.
4. Treat source normalization and dust/color variation explicitly. Do not call a fitted normalization an absolute luminosity prediction or an independent distance measurement.
5. Use the correct ESSENCE photon throughput and calibration, not the nearby-template UBVRI columns or DES filters as substitutes for observed R4m/I4m.
6. Preserve the source-model uncertainty when comparing to the provisional spectral aging estimate. The earlier SNID-library estimate and this empirical template are not guaranteed statistically independent.

These are analysis requirements. Interpolation and template-based forward modeling are known methods; no new physical formula is introduced here.

## Remaining input work

The full Narayan 2016 Table B2 response remains unacquired. The accessible public trees of gnarayan/source_synphot (3bc3d48217ad7ea5630131e68fd3c544d14d10f6) and gnarayan/sne-internal (cf5ecacded7510e0a588af11c733837ebdd55974) supplied no directly identified ESSENCE response files in this search. This does not establish that no public copy exists. The previous legacy-filter mismatch remains in force.

The next physical comparison needs the release-matched response and a declared source-variation treatment. All six research goals remain open; no observed fit is added by this audit.

## Reproducibility

run.py downloads the author archive only when missing, rejects a changed cached archive on subsequent runs, checks the old dependency hash, compares the spectral arrays exactly and checks both grids. The tar members are read as data without executing or extracting archive code. All hashes and ranges are in results.json.

Source archive: https://astrophysics.physics.fsu.edu/~hsiao/data/hsiao_template.tar.gz
