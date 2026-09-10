# Independent foreground reddening versus exposed residuals

All 164 fixed-coordinate IRSA requests succeeded. The primary reference-pixel reddening input spans E(B-V) = 0.0053 to 0.1965 mag. No galaxy was dropped or assigned a new distance, velocity or conversion rate. The comparison uses the already saved coarse-sky out-of-fold constant-baseline predictions, with signed residual defined as prediction minus observation in c-times-redshift units.

| Comparison with primary E(B-V) | Pearson correlation | Spearman correlation |
|---|---:|---:|
| Signed residual | -0.2295 | -0.2382 |
| Absolute residual | +0.1667 | +0.1301 |
| Region-demeaned signed residual | -0.1377 | -0.1456 |
| Region-demeaned absolute residual | +0.0895 | +0.0499 |

For the last two rows, both reddening and the relevant residual statistic have their own predefined-region mean subtracted. This describes within-region variation; it is not a physical correction or a full confounder adjustment. The eight individual region summaries are saved. Region 2, with the largest negative mean residual, does not have the largest mean reddening. No cut or region was selected to strengthen the correlation.

The modest association weakens after regional mean removal. This does not establish that dust causes the errors, nor that dust is irrelevant. Catalog calibration, distance distribution and other sky-dependent inputs remain potential confounders. Ordinary correlation significance tests treating all rows as independent would be inappropriate given the preceding regional sensitivity; no such p-value or causal significance is claimed.

## Source and physical meaning

The [IRSA DUST interface](https://irsa.ipac.caltech.edu/applications/DUST/docs/dustProgramInterface.html) supplies foreground Milky Way reddening at requested coordinates. Its SandF fields use Schlafly-Finkbeiner calibration and its SFD fields use Schlegel-Finkbeiner-Davis calibration. We retain both reference-pixel estimates and local means/scatters; they are not independent replicate measurements. The primary statistic was fixed before the batch. The local spatial scatter is not a per-galaxy uncertainty on the reddening correction. Only XML statistics were retrieved; generated image and extinction-table links were not followed.

This input is independent of our redshift residual construction: it is obtained from an external sky map at catalog coordinates. It is still a calibrated map estimate, not an assumption-free direct measurement. It concerns Galactic foreground dust, not intergalactic void density or gravity. It cannot be inserted as the required line-of-sight environmental term in the time mechanism merely because it is a sky-dependent quantity.

Dust changes the brightness and color used in distance work. But the published SBF distances may already include appropriate corrections, and the user's adopted values remain unchanged. A second correction would require tracing the original passbands, calibration and extinction treatment. The present correlation does not provide that correction, an independent galaxy motion, or a new photon-energy transfer law.

**Formula provenance:** Pearson correlation, Spearman rank correlation and subtraction of group means are established statistical operations. No new physical formula is proposed in this diagnostic.

## Reproduction and next action

Run `python research_work/results/redshift-priority/foreground_dust.py` for offline analysis of the saved metadata. Add `--fetch` only to refresh the external queries; current remote responses may differ. The script verifies the original 164-row source hash, complete unique group coverage, exact query coordinates, returned coordinates to one millionth of a degree, and finite numeric inputs. Source/protocol and response hashes are preserved, with parsed map values for every object. Raw HTTP response bytes are not archived; their hashes identify the retrieved responses but cannot be recomputed from the parsed projection alone.

Next trace the original SBF extinction/calibration treatment and other independent sky-dependent inputs, while retaining the same exposed sample and failed smooth-rate candidate. Do not fit dust into alpha or reinterpret foreground reddening as void exposure. No fresh outcome was opened, the ELVES screen remains five excluded/24 pending, and predictive improvement remains unproven.
