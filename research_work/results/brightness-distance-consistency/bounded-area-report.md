# Bounded beam-area response improves farther brightness transfer

The previous time-only test showed that assigning the brightness correction to extra event stretching worsens timing. This run instead retains b=1 and photon number while testing one bounded optical-area postulate. It improves both the existing brightness training and farther prediction groups. The physical mechanism remains unclosed.

## Formula and status

Known photon energy bookkeeping gives the transferred fraction

    f=1-E_obs/E_emit=1-1/(1+z).

The additional **phenomenological postulate** is

    D_G^2=D^2(1+eta f),
    F=L/[4pi D^2(1+z)^2(1+eta f)],
    D=ln(1+z)/alpha.

The two factors of 1+z represent energy loss and stipulated event stretching. The extra term describes beam spreading, not removal of photons. Its area increase saturates at 1+eta. The flux bookkeeping is known mathematics; the proportional area response is a proposed ansatz, not a derived metric, a graviton interaction or a claimed unique formula in the literature. Neither photon-number conservation nor this scalar flux equation establishes a conserved local stress-energy law for the field causing the spreading.

The fitted eta=0.615115, with alpha held at 0.0002488993286/Mpc and b=1. At z=1, area is 1.30756 times the Euclidean reference; the asymptotic ratio is 1.61512. These are conditional predictions, not measured beam areas. Source calibration includes the same proposed correction at calibrator distances, so the comparison does not insert a free brightness normalization per supernova.

## Actual standardized brightness comparison

The same 77 calibrator rows set common absolute magnitude; 466 noncalibrator rows at 0.1<=zHD<0.3 set eta. It is frozen for the 494 farther rows. Existing full measurement/calibration covariance and frame conventions are preserved. These are previously exposed standardized data, not a newly blind sample or a raw bolometric-flux fit.

| Group | Previous fixed model score | Bounded-area score |
|---|---:|---:|
| Training, 466 rows | 468.030 | 445.421 |
| Farther prediction, 494 rows | 424.482 | 402.335 |

Lower is better on this unchanged covariance. These are point-prediction GLS scores, not parameter-marginalized confidence levels or evidence ratios. The old eta=0 results are reproduced before scoring. The optimum is away from the declared bounds.

Correlated descriptive mean magnitude residuals:

| Redshift interval | Rows | Previous residual | Revised residual |
|---|---:|---:|---:|
| 0.1-0.3 | 466 | 0.12487 | 0.02448 |
| 0.3-0.6 | 365 | 0.18695 | 0.02326 |
| 0.6-1 | 104 | 0.20482 | -0.03700 |
| 1-3 | 25 | 0.34974 | 0.02811 |

These bins are not independent measurements. The bounded response avoids the large farther-distance overcorrection of the previous unbounded opacity/time reinterpretation, but the function was investigated after other failures on this dataset. Its improved transfer is exploratory and requires untouched evaluation before claims of predictive discovery.

## What this achieves and what it does not

There is now a simple conditional brightness relation with better farther predictions while leaving the redshift rate and adopted event-stretch relation unchanged. That is useful progress toward a joint phenomenological description. It does not derive the event timing, the beam spreading, or their connection from one local interaction. It also does not show that galaxy images remain sharp, that all spectral bands behave correctly, or that angular-size and lensing relations are consistent.

A global optical/metric model must generate this area relation with source and detector boundary conditions and energy accounting. An arbitrary distance-dependent area cannot be assumed to arise from a regular spacetime. The previously derived plane-lapse beam example is not automatically this new relation. Source standardization, color and population assumptions also need a forward-model treatment; the released magnitudes may retain fiducial modeling choices even though no expansion distance is used here.

Retain this as an unconfirmed candidate rather than replacing the paper's physical claims with a successful curve fit. Deposit formation/support, joint lensing, Milky Way structure, clusters and untouched mechanism-specific prediction remain unresolved. All six goals remain open.

Reproduce with `python research_work/results/brightness-distance-consistency/bounded-area.py`. bounded-area-results.json records scores and calibration; bounded-area-predictions.json records all 960 standardized magnitudes and both predictions. Input hashes are verified against the existing release.
