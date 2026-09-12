# Observed source-clock audit: SN 2011fe

## Finding

A single radioactive exponential is not an adequate description of this nearby event's late R-band light curve. Treating its measured decline as a direct nuclear clock would manufacture an apparent stretch of about 2.09. That number is **not a measurement of propagation stretch**. It demonstrates why a source model and spectral information are needed before using this band as an independent clock.

## Observations and provenance

We extracted all 14 LBT R-band measurements from Table 1 of [Dimitriadis et al. (2017), version 2](https://arxiv.org/html/1701.07267v2), retaining published Vega magnitudes and uncertainties. The underlying observations are attributed there to Shappee et al. (2016). The source HTML is cached with SHA-256 and the extracted rows are recorded in observations.json. These are published processed photometry, not newly reduced detector images or a bolometric dataset.

We use the actual MJD dates, 56337.48 through 57425.49. Subtracting 55814.30 only defines a plotting/window origin: no time rescaling, distance, expansion model, redshift, or photon-conversion parameter enters this fit. All selected rows remain, including the early brightening. This source and its behavior have been inspected, so this is exploratory evidence and cannot become the untouched distinguishing test.

## Formula status and calculation

**Known mathematics:** a constant exponential luminosity in an unchanged band corresponds to a straight line in magnitude:

m(t_o) = a + s (t_o - t_center).

We fit a and s by weighted least squares using the tabulated magnitude errors. No intrinsic scatter or extra uncertainty is fitted. The quoted chi-square is a descriptive diagonal-error discrepancy; absent systematic covariance prevents treating it as a complete significance calculation. Vega's constant magnitude zero point changes a, not s.

**Conditional consequence, not a new law:** if the band luminosity followed only exp[-t_o/(A tau)], the slope would imply A = 2.5 / [ln(10) tau s]. We use the same illustrative tau=111.3 day benchmark as the preceding feasibility calculation. It predicts 0.975 mag per 100 days when A=1. Neither single-isotope dominance nor a constant band fraction has been established here.

| Observer-date subset | Points | Fitted mag / 100 days | Diagonal chi-square / degrees of freedom | Misleading apparent A |
|---|---:|---:|---:|---:|
| Entire selection | 14 | 0.4663 | 1142.13 / 12 | 2.092 |
| Before origin + 900 days | 5 | 0.2330 | 738.07 / 3 | 4.187 |
| After origin + 900 days | 9 | 0.4086 | 13.34 / 7 | 2.387 |

The two windows are exploratory diagnostics, not selected validation samples. Formal slope errors and every prediction/residual are in results.json. In particular, the large whole-sample residuals invalidate reading its precise formal slope error as a precise clock calibration.

## Consequence for our research

We must not fit A freely to a single optical tail and report the result as evidence for our time mechanism. A constant stretch of a single exponential remains a single exponential, so that maneuver cannot explain the whole curve's structure anyway. This failure concerns an oversimplified emitting-source model; it does not distinguish expansion from companion transport, and it does not rule out a properly constrained radioactive clock.

The next observational requirement is multiband flux and spectra with a source energy-deposition model, including additional power and redistribution where needed. Keep nuclear lifetimes as physical hypotheses and let source constraints be checked nearby before applying the model to distant events. The DES artificial timing batch remains unchanged; this audit supplies real-data evidence for the source-model requirement, not completion of any of the six demonstrations.

Reproduce with `python research_work/results/sn2011fe-clock-audit/run.py`. Cached HTML is reused; deleting that single cache file deliberately reacquires the pinned-version URL and produces a new recorded hash if its representation changes.
