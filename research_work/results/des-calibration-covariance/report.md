# Published DES five-year zero-point covariance

The published Fragilistic zero-point covariance has been acquired, verified and connected to the historical calibration transform. This supplies measured cross-filter uncertainty instead of the artificial covariance used in the earlier software check. It does not complete a brightness likelihood or its full systematic-error budget.

## Provenance and selection

The [DES model documentation](https://des-sn-dr.readthedocs.io/en/latest/2_LCFIT_MODEL.html) explicitly links the [Fragilistic covariance release](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz). The associated [Brout et al. calibration paper](https://arxiv.org/abs/2112.03864) describes simultaneous cross-calibration of survey magnitude offsets. We use that calibration information without importing the paper's cosmological inference.

The file is pinned to commit c447f0fea703fcd0fff57de5000947b5ca81286b and verified against both Git blob SHA-1 and SHA-256. It loads with pickle disabled. Its actual dimensions are 102 by 102, with unique named entries. The full matrix is symmetric to 5.1e-21 and positive definite, with smallest eigenvalue 1.11e-6 in magnitude-squared units.

Both DES3YR and DES5YR entries exist. We select the exact labels DES5YR g/r/i/z, at zero-based indices 30-33, rather than relying on filter letters alone. The full matrix remains cached because comparisons with other surveys need cross-survey correlations, not just this four-by-four subblock.

## Extracted uncertainty

**Known covariance algebra, not a new physical law:** sigma_delta,j=sqrt(V_jj). The first-order fractional flux uncertainty from these offsets is 0.4 ln(10) sigma_delta,j.

| DES five-year band | Zero-point standard deviation (mag) | First-order fractional flux uncertainty |
| --- | ---: | ---: |
| g | 0.006051 | 0.5573% |
| r | 0.005695 | 0.5245% |
| i | 0.005726 | 0.5274% |
| z | 0.005686 | 0.5237% |

These describe this released zero-point component, not the total error of each supernova measurement. The matrix has nonzero cross-band terms; the report data retain all of them. Measurements in the same survey/filter share the corresponding offset uncertainty, so it is not independent noise that shrinks as the square root of the number of supernovae.

The previous Jacobian implementation was exercised with this measured subblock and illustrative predicted fluxes. It preserves cross-observation covariance and positive semidefiniteness. The illustrative fluxes are not observed supernova values. No data were corrected or fitted.

## What remains

The DES documentation also describes changing passband wavelengths and retraining source models when propagating calibration uncertainty. A zero-point-only matrix does not replace those effects. The compatibility of the pinned later photometry with the historical calibration still needs checking, as do the exact prescription for additional calibration components and correlations with a source model. Do not arbitrarily inflate, diagonalize or add an overlapping uncertainty component.

This acquisition enables a component of the joint redshift/timing/brightness calculation. It neither removes the source-evolution ambiguity nor supplies a physical cause for photon-companion transfer. The timing calibration continues independently, and all six scientific demonstrations remain incomplete.

```powershell
python research_work/results/des-calibration-covariance/acquire_and_check.py
```

The command reconstructs the immutable source cache, checks its hashes, writes the selected numerical subblock, and reproduces the propagation audit. It does not execute downloaded code.
