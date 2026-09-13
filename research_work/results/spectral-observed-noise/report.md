# Phase-difference recovery with actual distant-spectrum errors

This is a conditional synthetic calibration, not an age measurement of SN 2006mk. It uses the two real native wavelength grids and ERR arrays acquired from ESO, but substitutes labeled nearby spectral features for the supernova signal.

## Construction

- Select the earliest/latest eligible template epochs for each of 14 nearby objects with at least 20 labeled days of separation. Every match excludes all templates belonging to that object.
- Interpolate the processed library onto each observed rest-wavelength grid over 4000-5500 Angstrom, using catalog z=0.4754. Apply additional Gaussian broadening of resolving power 440; original template resolution is not deconvolved, so this is not an exact instrumental-resolution match.
- Set feature RMS to 0.1, 0.3 or 1.0 times the observed median SN flux. These are sensitivity settings, not measured feature amplitudes. Add independent Gaussian draws using actual ERR values on the native grid.
- Remove a quadratic continuum by weighted least-squares projection, correlate against all other objects and take the median phase of the best five distinct object matches. The projection is verified to remove an arbitrary quadratic to numerical precision.
- Use independent SeedSequence([seed, epoch]) streams. Draws are shared across amplitude settings for paired comparisons. Three seeds are an exploratory check, not confidence-interval calibration.

The observational wavelength grids contain [844, 835] retained samples. Median error divided by median flux is [0.131, 0.332]. This is not a feature-specific signal-to-noise measurement.

## Results

| Seed | Feature RMS / median flux | Median recovered / labeled age span | RMS span error (days) |
| --- | ---: | ---: | ---: |
| No added noise | 1.0 | 0.9822 | 2.723 |
| 3101 | 0.1 | 0.9747 | 5.541 |
| 3101 | 0.3 | 1.0038 | 4.687 |
| 3101 | 1.0 | 0.9737 | 4.034 |
| 3102 | 0.1 | 1.0193 | 7.244 |
| 3102 | 0.3 | 1.0086 | 2.947 |
| 3102 | 1.0 | 0.9780 | 2.333 |
| 3103 | 0.1 | 1.0101 | 7.238 |
| 3103 | 0.3 | 0.9874 | 5.119 |
| 3103 | 1.0 | 0.9754 | 2.811 |

Even with no added noise, object-to-object template differences produce a 2.72-day RMS span error. Across the final independent-epoch noise draws, weak-feature cases have RMS errors of 5.54-7.24 days, intermediate cases 2.95-5.12 days and strongest-feature cases 2.33-4.03 days. No recovered span is nonpositive in these cases. These are errors over labeled synthetic pairs, not confidence bounds for the real supernova.

For SN 2006mk, unchanged event intervals and a duration factor of 1+z differ by approximately 10.335 source days for the acquired pair. The present calibration does not establish whether its real spectra distinguish those hypotheses. Real feature contrast, continuum treatment, host residuals, correlated extraction noise and template uncertainty still matter. Nearby template phase labels also retain their historical calibration assumptions.

## Formula provenance

Weighted polynomial projection, Gaussian noise, Gaussian spectral smoothing and normalized correlation are known statistical/signal-processing operations. The amplitude sweep and matching choices are analysis assumptions. No new propagation or companion equation is derived here. No distant phase or timing exponent is fitted and no scientific holdout is opened.

## Next required step

Estimate the observed feature contrast with an explicit continuum/host model, test sensitivity to continuum complexity and error correlations, and calibrate at that range before fitting the real pair. Additional spectra and matched photometry are still needed for the joint redshift/timing/brightness objective. All six research goals remain open.

Input provenance is preserved in results.json and the [acquisition report](../essence-spectrum-acquisition/report.md).
