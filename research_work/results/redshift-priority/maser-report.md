# Fixed-parameter maser external diagnostic

The existing exponential rate was applied without refitting to all six galaxies in [Pesce et al. 2020, Table 1](https://arxiv.org/abs/2001.09213). Source table and frame note were visually checked on PDF page 4. The paper's distances are disk-model angular-size distances; they are adopted as fixed path distances under the user's stipulated-distance contract. We do not apply the paper's expansion-based distance-redshift equation or its flow corrections. Disk-model assumptions and matter/clock consistency remain relevant to a completed alternative theory.

## Formula provenance

The exponential fractional-change solution and linear control are established mathematics, not novel formulas. Their fixed coefficients are calibrated in this project on the earlier SBF sample. This diagnostic adds a comparison, not a physical derivation of the time-field rate or an originality claim.

## Results

| Galaxy | Distance (Mpc) | Observed z | Exponential prediction | c times residual (km/s) |
|---|---:|---:|---:|---:|
| UGC 3789 | 51.5 | 0.011074 | 0.012901 | +547.7 |
| NGC 6264 | 132.1 | 0.033999 | 0.033426 | -171.7 |
| NGC 6323 | 109.4 | 0.026023 | 0.027604 | +473.9 |
| NGC 5765b | 112.2 | 0.028439 | 0.028320 | -35.5 |
| CGCG 074-064 | 87.6 | 0.023924 | 0.022043 | -563.9 |
| NGC 4258 | 7.58 | 0.002266 | 0.001888 | -113.2 |

The exponential RMS is 384.3 km/s, versus 390.8 for the linear control. Exponential MAE is 317.6 km/s and signed bias +22.9 km/s. These are residual reporting units, not inferred velocities. The 6.5 km/s RMS difference is not a demonstrated significant improvement with six objects, shared calibration risks and unresolved nuisance effects. Two objects lie inside the prior 10.2-93.2 Mpc range; four lie outside, including one closer object. Within-range exponential RMS is 555.8 km/s; outside-range RMS is 258.9 km/s. Small, heterogeneous subsets cannot establish a general distance trend or all-distance accuracy.

Published distance errors are retained in the CSV but do not move the adopted values to improve residuals. Velocity measurement errors are far smaller than these residuals, so instrumental spectral precision alone is not an adequate prediction-error model. No validated interval coverage, acceptance threshold or fundamental-theory fit is claimed.

## Exposure and protocol accounting

The source's aggregate abstract result was viewed before the protocol; target values were subsequently viewed together during table extraction. Predictions and labels are separated in saved artifacts, but the prediction-file seal was created after that viewing. It is an integrity hash, not proof of blinding. This deviation from the ideal feature-before-label workflow is recorded in maser-seal.json. The fixed formula, coefficients and all-six selection were not changed after extraction.

An exact normalized-name repository scan found a prior NGC 4258 reference and no direct name hits for the other five outside the current audit. That is not an alias or group crossmatch. PGC identifiers, earlier uses under other names and common calibration links remain unresolved; no target is certified fresh. All six remain reported, with this limitation rather than a favorable subset selected. NGC 4258 can be a distance-calibration anchor for other methods, so geometric distance measurement does not imply independence of every cross-catalog calibration.

The source explicitly states CMB-frame optical velocities; Table 4 confirms its measured column is c*z. Only that measured column is used, not its inferred expansion velocities or peculiar-motion corrections. The latter were visible during source extraction but did not enter calculations.

## Reproducibility and next work

Run `python research_work/results/redshift-priority/maser_check.py`. The script writes source/protocol/prediction hashes, feature values, all-row CSV and metrics. The source PDF SHA-256 is recorded, and table values are transcribed explicitly with no optimized parameters. The script requires NumPy; no cosmological model library is used.

This advances external distance testing but does not satisfy genuinely unexposed validation, independently constrained motion errors or demonstrated predictive improvement. Finish alias/group matching and identify a separate sealed sample. In parallel, the proposed time-field mechanism must supply a rate and physical clock response before an empirical exponential fit can be interpreted as that mechanism's prediction.
