# Conditional joint source-history pilot for SN 2006mk

The spectral-pair source ages (-3.4, +16.7 days) and an unchanged empirical Hsiao source history predict R/I temporal flux ratios 0.367/0.618. Observed ratios are 0.354/0.642, with incomplete diagonal errors 0.054/0.152. No source amplitude, width, dust, color or age was fitted to these photometric ratios. This is an encouraging conditional cross-observable comparison, not a release-calibrated validation, blind prediction or new physical mechanism.

## Why a preliminary calculation is useful

Full 2016 ESSENCE throughput remains unavailable in this workspace. An additional arXiv-source acquisition confirms that its TeX table is abbreviated too. The legacy curves are therefore used only as declared comparison inputs, with illustrative wavelength shifts and negative-lobe handling retained. This changes no previous conclusion about their calibration mismatch. It does not claim that those sensitivities bound the missing response uncertainty.

Three source trajectories share the same starting phase inferred from the early spectrum:

1. Use both independently matched spectral phases (-3.4, +16.7).
2. Use late phase = early phase + observer interval/(1+z), with z=0.4754.
3. Use late phase = early phase + observer interval.

The photon wavelength transformation is the same in every case; only the event time mapping differs. The first case is a spectral-to-photometric comparison; the other two are conditional timing controls. The data were already exposed, so none is a preregistered prediction on untouched observations.

## Numerical predictions

| Source trajectory | Band | Unshifted legacy prediction | Range across declared filter cases | Measured ratio |
| --- | --- | ---: | --- | ---: |
| spectral_pair | R | 0.3669 | 0.3600-0.3738 | 0.3536 +/- 0.0545 |
| spectral_pair | I | 0.6185 | 0.6159-0.6188 | 0.6415 +/- 0.1525 |
| event_stretch_1_plus_z | R | 0.3120 | 0.3054-0.3188 | 0.3536 +/- 0.0545 |
| event_stretch_1_plus_z | I | 0.5681 | 0.5667-0.5681 | 0.6415 +/- 0.1525 |
| unchanged_arrival_intervals | R | 0.1319 | 0.1275-0.1366 | 0.3536 +/- 0.0545 |
| unchanged_arrival_intervals | I | 0.3461 | 0.3427-0.3480 | 0.6415 +/- 0.1525 |

The sensitivity cases shift the photon throughput by -50, 0 and +50 Angstrom and either retain or clip the small negative interpolation lobes. They are illustrative choices, not measured passband uncertainties. No case is selected because it fits better.

## Formula provenance and computation

Use the known synthetic-photometry integral C_b(t) proportional to integral L_lambda(t, lambda_o/(1+z)) lambda_o T_b(lambda_o) d lambda_o. Common distance, amplitude and constant redshift prefactors cancel in the temporal ratio. The legacy response includes a wavelength factor, so divide it by its original wavelength to recover a proportional photon throughput before shifting it; multiply by the shifted observer wavelength once during integration. This is known detector bookkeeping, not a novel propagation law.

The source is linearly interpolated within the verified author grid, with no extrapolation. The evaluator reproduces an exact tabulated spectral phase. Source identity is checked against the author audit and every input file is hashed. A repeated complete run reproduces results.json exactly.

## Interpretation and limits

The two-band decline is compatible with the spectral-age-based empirical history in this pilot. The unchanged-interval trajectory predicts more fading, especially in R. Do not convert this discrepancy into a sigma rejection: source diversity, phase uncertainty, template correlations, dust, exact passbands, shared calibration and non-simultaneity are omitted. The hypotheses are not exhaustive source models.

An event stretch of 1+z also produces reasonably close ratios in this restricted comparison, but has been imposed as a time map. It is not derived from photon energy transfer. No distance-redshift fit, absolute luminosity prediction, companion supply calculation, deposit model or novel mechanism-specific test is supplied by these ratios. All six research goals remain open.

Next obtain release-matched throughput and propagate source/phase variations jointly across the full light curve, rather than tuning each epoch. This pilot justifies that work by showing that the spectral and photometric measurements can be connected numerically.

Sources: [verified source template](../hsiao-source-audit/report.md), [spectral ages](../spectral-real-pair/report.md), [observed brightness](../sn2006mk-joint-photometry/report.md), [legacy-filter mismatch](../essence-passband-audit/report.md).
