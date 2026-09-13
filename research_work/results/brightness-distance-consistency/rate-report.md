# Can one recalibrated conversion rate repair brightness and redshift together?

**Result:** a lower alpha fitted to the lower-redshift supernova group improves frozen central brightness predictions for the farther group, but shifts away from the nearby-galaxy rate and leaves the brightness trend between redshift groups unchanged. A single constant-rate adjustment is not a demonstrated joint solution.

## Formula and calibration

Retain the existing phenomenological exponential shift, photon-number conservation and equal event/wavelength stretching:

    1+z=exp(alpha D),  S_t=S_E,  F=L/[4 pi D^2 S_E^2].

These are the previously specified model assumptions with known integration and flux bookkeeping, not new fundamental laws. No opacity parameter is included. The fit still does not derive the interaction causing alpha.

Use the same existing standardized Pantheon+ magnitudes, full covariance and conditional frame convention as the [opacity comparison](opacity-report.md). Calibrate luminosity on 77 Cepheid-host rows, treating their distance moduli as stipulated geometric distances and including the proposed propagation dimming there. Fit alpha on 466 noncalibrators at 0.1<=zHD<0.3. Freeze it for 494 rows at higher redshift, then apply it to all 164 already exposed nearby galaxy-distance measurements. All samples are retrospective; no final untouched holdout is opened.

Let Dbar_c be the GLS-weighted calibrator distance. Relative to alpha_0, every high-redshift predicted magnitude changes by exactly

    Delta m(alpha) = -5(alpha-alpha_0)Dbar_c/ln(10)
                    -5 log10(alpha/alpha_0).

This follows from M(alpha)=M_base-5 alpha Dbar_c/ln(10) and D=ln(1+zHD)/alpha. The shift is independent of the target's redshift. Fit one GLS offset in the training group and invert this monotone expression; the rate is constrained to c alpha between 20 and 120 km/s/Mpc. That declared bound is not reached.

## Results on actual retained measurements

Brightness training prefers alpha=0.0002350796784/Mpc, or **c alpha=70.4751 km/s/Mpc**. The formal delta-chi-square-one interval is 69.4712 to 71.4934. The original galaxy fit gives 74.6181, with its previously computed training-tile bootstrap 95% interval 72.4172 to 76.8874. These are different uncertainty constructions, with shared calibration and model dependencies; nonoverlap must not be converted into an independent Gaussian tension claim. The fitted rate is about 5.6% smaller.

| Brightness sample | Rows | Original rate chi-square | Frozen brightness-rate chi-square |
|---|---:|---:|---:|
| 0.1<=z<0.3, fitted | 466 | 468.030 | 452.170 |
| z>=0.3, transfer comparison | 494 | 424.482 | 391.842 |

Both columns use the same covariance for central predictions; the improvement is not produced by widening the prediction covariance. These scores are not a full predictive likelihood, and the transfer column does not include the fitted rate's uncertainty. Training uses one fitted coefficient. A lower chi-square than the row count is not proof of physical correctness.

After the fit, GLS mean magnitude residuals are 0.0000, +0.0621, +0.0799 and +0.2249 in the respective redshift ranges 0.1-0.3, 0.3-0.6, 0.6-1 and 1-3. Positive residuals mean observed standardized sources are fainter than predicted. The differences from the first group have formal covariance-derived errors 0.0127, 0.0237 and 0.0646 mag. They are correlated and inherit the published reduction; they are not three independent new detections.

**Those inter-group contrasts cannot change under any constant alpha adjustment in this model.** Their persistence is an exact consequence of the common magnitude shift, not a failed optimizer. Other physical effects or reduction assumptions could affect them, but have not been fitted here.

## Transfer back to nearby galaxy redshifts

Keep every published galaxy distance fixed and predict z=exp(alpha D)-1 with the brightness-fitted rate. No galaxy data refit is performed.

| Historical exposed galaxy partition | Rows | Original rate RMS, km/s | Brightness rate RMS, km/s |
|---|---:|---:|---:|
| Training | 104 | 457.577 | 489.419 |
| Validation | 35 | 437.065 | 425.982 |
| Test | 25 | 415.414 | 406.897 |

Units are c times redshift residual, not inferred speeds. The galaxy training fit worsens while the two smaller reused partitions improve slightly. We must report all three rather than selectively call this a galaxy success or universal failure. The brightness-derived rate produces negative average redshift residuals in all three, approximately -178, -60 and -91 km/s.

## Interpretation and limits

This tests the same coefficient across brightness and distance-redshift observations. It improves the farther brightness comparison without invoking photon removal, but does not establish one rate that resolves both observables or their residual structure. The previous stationary-conversion timing failure also remains: the flux law presupposes a working event-stretching mechanism.

As in the previous test, standardized B-band magnitudes are not direct bolometric measurements. SALT2, host/color, selection, velocity corrections and covariance retain published assumptions. The Cepheid calibration is stipulated here rather than rederived. The distant path lengths are inferred from our formula, not independent distance measurements; no expansion-based distance law is used. These qualifications limit interpretation and cannot be hidden by fitting a new constant.

**Decision:** retain this as a measured calibration tradeoff, not a replacement confirmed rate. Any further change must predict the shape across redshift or independently justify a calibration change; a constant alpha alone cannot repair the inter-group trend. rate.py records both sets of per-object predictions, covariance-based scores, and input hashes. All six objectives remain open.
