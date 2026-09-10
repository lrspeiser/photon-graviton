# Why foreground reddening is not a ready-made distance correction

## What the source audit establishes

[Tonry et al. 2001](https://arxiv.org/abs/astro-ph/0011223), page 4, describes replacing earlier extinction estimates with SFD maps. Both fluctuation magnitudes and colors were corrected. The distance zero point was also recalibrated using revised extinction estimates for calibrators, limiting the overall distance change. Their already corrected photometry must not be treated as uncorrected observations.

[Cosmicflows-4](https://arxiv.org/abs/2209.11238), sections 5.1-5.2, combines ground-based optical Tonry data, HST Virgo/Fornax data, the Cantiello extension and infrared HST data. It compares overlapping galaxies and aligns zero points; the full integration uses its global calibration procedure. Therefore the adopted CF4 SBF column is not one raw I-band survey.

The source-level chain is clearer, but we have not yet mapped every one of our 164 representatives to its contributing surveys, passbands, original extinction estimates and calibration weights. No source is assigned to a row merely from its distance or sky direction. Local primary PDFs are identified by SHA256 4870ea7a4902e14b78642a46047d1c011ed6b540414a061122954b8af40c9bca (Tonry) and 5562be287b6cadd3034a1e4c63a4488de16eb0591de794989506e068ec4477aa (CF4).

## Established calibration algebra, not a new formula

Write observed fluctuation magnitude as m_bar_obs and observed color as C_obs. Let E = E(B-V), and let R_V and R_I convert E into extinction in the two relevant bands. A linear SBF absolute-magnitude calibration is M_bar = a + b(C_0 - C_pivot), where C_0 = C_obs - (R_V-R_I)E. The familiar I-band color-calibration form is established in the [SBF survey calibration literature](https://arxiv.org/abs/astro-ph/9609113). Filter definitions and coefficients must match the actual source.

**Provenance: established SBF and extinction relations, with an algebraically derived sensitivity; not unique to this research.** At fixed observed photometry:

mu = m_bar_obs - R_I E - a - b[C_obs - (R_V-R_I)E - C_pivot],

d mu/dE = K = b(R_V-R_I) - R_I.

Brightness correction alone contributes -R_I, while the color-based luminosity calibration contributes b(R_V-R_I). Their competition means the sign cannot be inferred from dimming alone. For fixed b and extinction coefficients, changing the adopted reddening by delta_E and zero point by delta_a gives the exact finite difference:

delta_mu = K delta_E - delta_a.

**Provenance: conditional algebra using the same established calibration.** If calibrators share sensitivity K, their distance scale is held fixed, and their effective reddening revision is delta_E_cal, then the required zero-point revision is delta_a = K delta_E_cal. Hence:

delta_mu = K(delta_E - delta_E_cal).

A common reddening revision cancels under those conditions. A different passband, changing slope, heterogeneous calibrator weights or a different absolute distance scale changes the calculation. We have not assumed those conditions hold for the whole CF4 compilation.

**Provenance: established distance-modulus definition.** The resulting conditional distance ratio is D_new/D_old = 10^(delta_mu/5). It is not applied to the user-adopted distances here.

## Consequence for the residual audit

The dust-map audit measured total foreground E, not the error or revision delta_E in the original survey's correction. It also did not measure the calibrator revision or covariance. Even a stronger residual correlation with total E would not supply those missing quantities. Treating the correlation as a distance correction would combine an unknown calibration change with a residual-based adjustment, obscuring rather than independently testing the proposed propagation law.

This provides a specific next data requirement: recover per-row survey provenance and the original photometric/extinction calibration before evaluating any dust-revision sensitivity. Preserve all catalog distances and the frozen unsuccessful candidate comparisons. Foreground dust still does not provide intergalactic void exposure or an independently measured galaxy motion.

Run `python research_work/results/redshift-priority/sbf_extinction_response.py`. Four symbolic identities check the derivative, finite change, zero-point substitution and common-shift cancellation. No target outcome is read and no coefficient is fitted. These checks validate the algebra under its assumptions, not the physical completeness of our time model or the accuracy of a catalog correction. The full fresh-prediction and credible-uncertainty goals remain open.
