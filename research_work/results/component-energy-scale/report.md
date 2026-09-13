# Is the Cepheid mismatch only an overall deposit-strength problem?

**A shared amplitude helps substantially but does not reproduce every bin.** Scaling the fixed added source by 1.7194 reduces RMS from 28.41 to 13.16 km/s on the nine exposed Cepheid proxy bins. Individual bins require scales from 1.2834 to 2.5777. The minimum possible largest absolute residual for one nonnegative scale is 20.454 km/s. These are conditional calibration results, not energy-supply predictions or significance statements.

## Formula provenance and fitting

Known linearity of Poisson gravity and the unchanged bin-averaging prescription give v_i(k)^2=B_i+k D_i, where B_i is the ordinary-matter squared balance speed and D_i is the previous added squared-speed contribution. The source shape, ordinary matter and all seven component ratios are fixed; k is a shared multiplier. The identification of k with exposure or stored companion energy requires the assumed linear gravity response and fixed distribution. No microscopic energy conversion efficiency is measured here.

For positive B,D and observed proxy y, mean squared residual has derivative mean[D(1-y/sqrt(B+kD))] and positive second derivative mean[y D^2/(2(B+kD)^(3/2))]. These known calculus identities give the unique nonnegative least-squares optimum, verified independently with scalar minimization. Each bin's exact-match multiplier is (y_i^2-B_i)/D_i. The minimax solution balances the largest positive and negative residual; monotonic predictions establish the global minimum worst residual for this one-parameter family.

| Radius mean (kpc) | Required individual scale | Inferred proxy (km/s) | Shared RMS-fit prediction | Residual (km/s) |
|---|---:|---:|---:|---:|
| 6.574 | 1.5890 | 241.70 | 246.01 | +4.31 |
| 7.516 | 1.4272 | 236.02 | 246.13 | +10.11 |
| 8.529 | 1.3800 | 231.81 | 243.79 | +11.98 |
| 9.466 | 1.2834 | 224.36 | 239.87 | +15.51 |
| 10.514 | 1.4658 | 225.41 | 234.20 | +8.78 |
| 11.517 | 1.8356 | 231.79 | 227.95 | -3.84 |
| 12.552 | 2.1556 | 234.85 | 221.14 | -13.71 |
| 13.389 | 2.1242 | 228.00 | 215.56 | -12.44 |
| 14.313 | 2.5777 | 234.32 | 209.45 | -24.87 |

The distinct minimax optimum is k=1.86467; its worst residual is 20.454 km/s. It minimizes a different loss than RMS and is not substituted into the RMS comparison. Repeating both calculations using the coarse gravity predictions changes the RMS optimum by about 0.000245 and the optimized RMS by 0.000695 km/s.

## What can and cannot be concluded

The fixed-source failure partly reflects normalization. Its shape can approximate these exposed velocity proxies after one more fit, but no single amplitude matches all bins exactly. Whether residuals are acceptable requires the missing shared distance, matter-model, selection and dynamical uncertainties. The earlier empirical-field RMS of 15.36 km/s used a different calibration history; the newly fitted 13.16 is not a fair blind improvement claim.

All nine bins are retained and already exposed. Final-test stars are unopened. A 72% increase in effective source strength is not evidence that a sufficient photon reservoir exists, nor proof that capture can form this spatial distribution. This curve also does not predict vertical motions or lensing. Future evaluation must carry the same calibrated source into those observables and independently determine its supply and dynamics. All six scientific objectives remain open.
