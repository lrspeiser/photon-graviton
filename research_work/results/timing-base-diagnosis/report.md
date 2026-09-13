# Frozen timing calibration: completed base-case diagnosis

All 160 predeclared base simulations are complete. This report concerns those base cases only; the eight selected numerical refinements are still being completed by the original live run. No running code, protocol, threshold or saved base outcome has been changed. The base-case content is separately hashed so subsequent refinement writes cannot silently alter this diagnosis.

## Outcomes

| Shape | SNR | Injected b | Mean bias | Invalid fits | Accepted 95% truth membership | LR-only truth membership | Preset cell screen |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| split_gaussian | 5 | 0 | 0.0489 | 7 | 13/20 | 20/20 | False |
| split_gaussian | 5 | 1 | 0.0159 | 4 | 16/20 | 20/20 | False |
| split_gaussian | 20 | 0 | 0.0079 | 0 | 20/20 | 20/20 | True |
| split_gaussian | 20 | 1 | 0.0056 | 0 | 20/20 | 20/20 | True |
| split_gaussian_shoulder | 5 | 0 | 0.0121 | 4 | 16/20 | 20/20 | False |
| split_gaussian_shoulder | 5 | 1 | 0.0066 | 1 | 18/20 | 19/20 | False |
| split_gaussian_shoulder | 20 | 0 | 0.0019 | 0 | 20/20 | 20/20 | True |
| split_gaussian_shoulder | 20 | 1 | 0.0108 | 0 | 19/20 | 19/20 | True |

The four high-signal cells pass the declared base screening criteria; the four low-signal cells fail. Passing a 20-seed cell does not establish precise 95% coverage. Seeds are paired across cells, so the counts cannot be pooled as independent trials. These are artificial fluxes on real cadence patterns, not measured cosmological timing.

## What fails

All 16 invalid fits have successful optimizers but inferred intrinsic log-width scatter at the imposed floor of 0.03, violating the preset interior-solution requirement. No result is reclassified as valid here. The fitted timing slopes can be close to their injected values even when the scatter estimate lies at this boundary.

The accepted-truth count requires both a valid fit and a profile-likelihood ratio below the nominal chi-square threshold. LR-only membership removes the first condition solely as a diagnostic. It is not replacement coverage and is not grounds to remove the validity rule. Neither count by itself supplies a calibrated confidence region for a boundary model.

## Required repair

The existing continuous-scatter integration development already permits exactly zero scatter and avoids forcing narrow populations onto width-grid nodes. That component must be incorporated with proper full-domain support handling and reliable optimization, then recalibrated as a new version. Its earlier two-sample test and narrower numerical search bounds do not certify all 160 cases.

A repaired estimator must retain comparisons with the original outcomes, include zero scatter without undefined log(0), handle finite width support explicitly, and validate its profile statistic at the boundary. Simply reducing the floor, discarding failed cases, or presenting LR-only counts as success would not establish reliability. Final refinement results are still needed to separate integration sensitivity from the population-scatter issue.

This is a failure of the current pipeline screening for faint simulated events, not evidence against the proposed photon-companion physics. No new observed timing estimate is reported. All six scientific goals remain open.

## Verification

All frozen source/protocol/input hashes are checked against the live checkpoint, and exactly 160 base records are required. The asserted scatter-floor diagnosis is checked on every invalid fit. A repeated run requires the same canonical base-case hash. Standard profile-likelihood and covariance concepts used here are known statistical methods, not new physical equations.

[Frozen experiment](../timing-coverage-calibration/protocol.json) | [Continuous-scatter development](../timing-continuous-scatter/report.md)
