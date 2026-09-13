# SN 2006mk: photometry linked to the actual spectral pair

The published ESSENCE q106 light curve is acquired and preserved with its original header and SHA-256 hash. The updated catalog identifies q106 as SN 2006mk. All 72 observations are retained: 44 R4m and 28 I4m, including 22 negative-flux measurements. Negative difference flux is not negative emitted energy and is not discarded or converted into magnitudes.

This joins actual brightness measurements to the same object used in the exploratory spectral-aging analysis. It is not a source-model fit or validation of a redshift-distance law.

## Closely timed measurements

Each selected photometric observation lies within 0.14 day of its spectral midpoint. They are close in time, not simultaneous; no interpolation or assumption of zero evolution across that offset is applied. All other measurements remain available for a future full temporal likelihood.

| Band | Early flux | Late flux | Late/early | Magnitude change |
| --- | ---: | ---: | ---: | ---: |
| I4m | 12.5212 +/- 1.2862 | 8.0325 +/- 1.7215 | 0.642 +/- 0.152 | 0.482 +/- 0.258 |
| R4m | 11.8771 +/- 0.5890 | 4.1999 +/- 0.6124 | 0.354 +/- 0.054 | 1.129 +/- 0.167 |

Fluxes use the published zeropoint 25 convention; they are not erg/s/cm2. The familiar magnitude relation m=25-2.5 log10(F) is known photometric mathematics, not a new model equation. Quoted errors are first-order propagation using the listed marginal errors as independent. They omit shared template/zeropoint covariance and the small epoch mismatch, and are not complete uncertainty intervals. No significance claim is made.

The difference in R-I between these epochs is +0.647 mag, with a conditional diagonal uncertainty 0.307 mag. This describes a changing observed band ratio, not an independently established new propagation effect; intrinsic supernova evolution is part of the forward prediction. Constant within-band multiplicative calibration cancels from each temporal ratio, but shared additive difference-flux baselines do not.

## Combined target now available

- Measured host redshift: updated catalog z=0.4754, uncertainty 0.0001; the older photometry header says 0.4750. Retain that version difference instead of silently overwriting either file.
- Actual spectral midpoint interval: 32.073696 observer days.
- Exploratory matched spectral age span: 20.1 days, with substantial calibration limitations documented separately.
- Observed R4m and I4m temporal flux ratios: 0.354 and 0.642, with the close-but-not-identical observing epochs preserved.

## Why this does not determine the propagation normalization

For a fixed path and constant propagation parameters over the event, write the predicted band flux as F_b(t_o)=K times B_b(t_e; z), where B_b is the source spectrum integrated through the appropriately shifted detector response. This is a conditional forward-model decomposition, not a new physical law. The time map determines t_e. The common time-independent factor K cancels from F_b(t_2)/F_b(t_1).

Consequently these ratios test source evolution and its time mapping, but cannot by themselves establish the absolute energy loss, distance or capture efficiency. Using z to infer distance and then calling that distance an independent redshift test would be circular. Absolute brightness needs a justified source-luminosity model, independent distance information where available, the actual ESSENCE passbands and their calibration uncertainties. The existing DES passbands are not interchangeable with R4m/I4m.

The next joint calculation must predict the two-band evolution and spectra from a shared source history under each proposed arrival-time mapping. It must not separately tune a luminosity at each epoch to match the observations. No new companion parameter or age/size constraint is imposed by this acquisition.

## Provenance and checks

Source: https://cdsarc.cds.unistra.fr/ftp/J/ApJS/224/3/lcs/q106.W6yr.clean.nn2.Wstd.dat

Catalog documentation: https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJS/224/3?format=html&tex=true

Narayan et al. 2016, ApJS 224, 3, ESSENCE light curves. Original processing header reports pruning, error adjustment and baseline settings; the released light curve is a processed measurement product. No full measurement covariance is present in this six-column file.

The parser verifies observation IDs are unique, all supplied errors are positive and the selected measurements are within the stated spectral-epoch tolerance. Cached input changes are rejected by hash on subsequent runs. Re-running reproduces the output exactly. All six research objectives remain open.

[Spectral pair](../spectral-real-pair/report.md) | [Spectral acquisition](../essence-spectrum-acquisition/report.md)
