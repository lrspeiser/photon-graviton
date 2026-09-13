# Timing interval sensitivity to event-grid resolution

For both exposed artificial samples, we profiled mean duration and intrinsic scatter at fixed timing exponent, using the same bounds and finite-support constraint as the revised estimator. Three optimizer starts were retained per evaluation. We located lower and upper crossings of the nominal chi-square(1) 95% likelihood-ratio threshold on both 321- and 641-node grids.

The threshold and profile likelihood are known statistical methods, not new physics. These endpoints are nominal local profile crossings, not empirically calibrated confidence limits. The search does not establish absence of distant disconnected likelihood regions or prove global nuisance optimization.

| Sample | Grid nodes | Lower crossing | Upper crossing |
|---|---:|---:|---:|
| split_gaussian-snr5-b1-seed902 | 321 | 0.68736417 | 1.36674287 |
| split_gaussian-snr5-b1-seed902 | 641 | 0.68675737 | 1.36759087 |
| split_gaussian-snr5-b0-seed903 | 321 | -0.39908655 | 0.29675415 |
| split_gaussian-snr5-b0-seed903 | 641 | -0.39908664 | 0.29596363 |

Maximum absolute crossing change: 0.00084799974.

All input arrays and source hashes are recorded. Every optimizer trial is retained, including unsuccessful or inferior fits. These two numerical checks do not establish population-wide coverage, nor do they validate a physical interpretation of redshift. The full 160-case development recalibration remains separate and unchanged. Wider grid checks and fresh-seed calibration remain necessary. All six goals remain open.
