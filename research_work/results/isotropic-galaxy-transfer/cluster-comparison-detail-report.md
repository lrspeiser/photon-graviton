# Coma comparison: why the NFW result looked worse

The earlier outer-only score is not a full-profile goodness-of-fit score. The authors fitted all six bins; we initially fitted only three. This diagnostic retains every measurement, including the negative outer shear, and uses the same two-parameter shapes. Point mass is a diagnostic control, not a physical cluster candidate.

## published

| Shape | Fit all six | Predict outer three from inner three | Leave-one-bin-out sum |
|---|---:|---:|---:|
| transparent | 3.6839 | 2.7018 | 6.8573 |
| strong_interception | 3.7245 | 2.5024 | 10.7919 |
| nfw | 3.8549 | 4.7106 | 6.8209 |
| mond_point_baryons | 4.4575 | 4.8299 | 19.6181 |
| point_mass_control | 7.3607 | 2.3261 | 32.1738 |
## blank_subtracted_sensitivity

| Shape | Fit all six | Predict outer three from inner three | Leave-one-bin-out sum |
|---|---:|---:|---:|
| transparent | 2.2554 | 2.0964 | 3.0469 |
| strong_interception | 2.2605 | 1.9081 | 3.6862 |
| nfw | 2.5883 | 3.8452 | 8.5905 |
| mond_point_baryons | 3.1834 | 3.3011 | 25.8690 |
| point_mass_control | 6.3649 | 1.9777 | 36.1069 |

## Interpretation

NFW is competitive: its full-profile score is 3.8549 versus 3.6839/3.7245 for the companion shapes. Its leave-one-out score 6.8209 is essentially tied with transparent capture 6.8573 and better than strong interception 10.7919. These differences do not establish superiority. The point-mass control deteriorates to 32.1738 on leave-one-out prediction.

For the original inner-only NFW fit, outer bins 4 and 6 lie 1.28 and 1.74 plotted standard errors below its predictions; bin 5 differs by only 0.20. NFW fitted to all six prefers scale 0.431 first-radius units, versus 1.006 from the inner three. This is why extrapolating the inner fit looked worse: the three inner bins do not fix the scale tightly. It is not a failed NFW projection.

All sums use diagonal plotted errors. Lower is closer under these assumptions; leave-one-out predicts each omitted bin using five others and sums those six errors. It is an interpolation/stability diagnostic, not an independent next-cluster prediction. No probability or model-selection significance is claimed. The blank-subtracted rows assume an additive control offset and independent errors, and are sensitivity results rather than a mandatory correction.

The published NFW full-profile chi-square is 3.87 for four degrees of freedom. Our reconstructed-bin, bin-center fit is only an approximate reproduction; exact source weighting/bin averaging is unavailable. We independently checked the numerical NFW projection against the known analytic projected profile, with maximum relative error 7.77e-12. The prior inner fits are reproduced to the checks recorded in the script.

A point mass places all lens mass inside the measured radii and has shear proportional to 1/R^2. It can be an exterior approximation, but these six points do not demonstrate that Coma is a point mass or that its extended gas and galaxy motions can be explained that way. Keep it out of the physical candidate ranking.

The companion alternatives and compact-baryon MOND remain freely normalized shapes. Neither is a baryon-calibrated cluster prediction; the one-third retention factor is still absorbed by the fitted amplitude. See the previous protocol.

Sources: [Kubo et al., section 3.2](https://arxiv.org/pdf/0709.0506); [known analytic NFW lensing, Wright and Brainerd](https://arxiv.org/abs/astro-ph/9908213). All six research goals remain open.
