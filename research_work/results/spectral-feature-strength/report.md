# Observed spectral feature-strength audit

This audit measures residual structure in the two acquired SN 2006mk spectra. It does not measure supernova age, establish a timing law, or isolate supernova features from every systematic effect.

## Method and formula provenance

All operations below are known weighted least squares and noise bookkeeping, not new physical laws. In the 4000-5500 Angstrom approximate rest window, fit polynomial continua of degree 2, 3 and 5, with and without an additional fixed basis vector from the corresponding host extraction. The host coefficient is a nuisance parameter, not a physical contamination fraction.

With C=diag(sigma_i^2), Q an orthonormal basis for the whitened design, and P=I-QQ^T, the residual is r=C^(1/2) P C^(-1/2) f. Its expected noise power is tr(C P). Subtract this from sum(r_i^2) to calculate a signed unweighted excess variance. Preserve negative estimates; the separately reported nonnegative RMS clips them at zero only for display.

Also report a differently weighted diagnostic: [sum(r_i^2/sigma_i^2)-(n-k)] / sum(1/sigma_i^2), divided by median flux squared. This is a weighted residual-energy measure, not the same estimator or an intrinsic feature variance under arbitrary wavelength dependence. Both assume diagonal noise and a fixed design. The host extraction has noise and may share extraction errors, so that assumption is conditional.

## Results

RMS fractions below are relative to median observed SN flux. They describe structure remaining after continuum projection.

| Epoch MJD | Polynomial degree | Host basis | Unweighted excess RMS | Weighted excess RMS | chi-square / dof |
| --- | ---: | --- | ---: | ---: | ---: |
| 54031.20605 | 2 | False | 0.145 | 0.200 | 3.298 |
| 54031.20605 | 2 | True | 0.141 | 0.197 | 3.219 |
| 54031.20605 | 3 | False | 0.145 | 0.200 | 3.293 |
| 54031.20605 | 3 | True | 0.141 | 0.196 | 3.217 |
| 54031.20605 | 5 | False | 0.114 | 0.168 | 2.631 |
| 54031.20605 | 5 | True | 0.111 | 0.167 | 2.614 |
| 54063.27974 | 2 | False | 0.000 | 0.346 | 2.446 |
| 54063.27974 | 2 | True | 0.000 | 0.345 | 2.442 |
| 54063.27974 | 3 | False | 0.000 | 0.345 | 2.445 |
| 54063.27974 | 3 | True | 0.000 | 0.344 | 2.441 |
| 54063.27974 | 5 | False | 0.000 | 0.324 | 2.277 |
| 54063.27974 | 5 | True | 0.000 | 0.324 | 2.279 |

The early spectrum has unweighted excess RMS about 0.11-0.15 of its median flux. The late spectrum has a negative unweighted excess variance in every case, while its weighted statistic remains positive. This is not evidence that the late spectrum lacks features: high-error pixels dominate the unweighted expected noise power, whereas lower-error pixels dominate chi-square. Neither diagnostic alone calibrates its age precision.

Under the quadratic/no-host case, weighted residual RMS fractions are 0.200 and 0.346. These suggest that the previous strongest-feature injection (RMS equal to median flux) cannot simply be assumed representative. However, the injected amplitude was defined before polynomial removal; these measured residual amplitudes are after removal. A direct numerical equality between them would be unjustified. Template-dependent transfer through the same projection is needed.

Host-basis inclusion changes the quadratic weighted diagnostic from 0.200 to 0.197 in the early spectrum and 0.346 to 0.345 in the late spectrum. This small change does not exclude host contamination: the fixed host basis cannot represent every extraction residual, and its uncertainty is not propagated. Changing polynomial degree removes some real or systematic broad structure.

## Verification

Direct least-squares residuals agree with the orthogonal projection for all 12 cases. For the final design, 4,000 independent Gaussian draws reproduce the analytic residual-noise power to relative error 0.000148. These checks verify linear algebra under the stipulated noise model; they do not validate the archive errors or observational likelihood.

## Required next step

Map the measured weighted residual diagnostic to each candidate template through the same continuum projection, then calibrate recovery using epoch-specific amplitudes and alternative covariance/host treatments. A real phase estimate must retain template and preprocessing uncertainty. No measured aging rate is claimed here. All six research goals remain open.

Sources and FITS provenance: [ESO acquisition audit](../essence-spectrum-acquisition/report.md). Input hashes and all signed statistics are in results.json.
