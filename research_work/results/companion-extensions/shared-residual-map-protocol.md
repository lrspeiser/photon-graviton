# Quick shared-law residual diagnostic

Fixed before calculating this diagnostic, 13 September 2026. Existing samples
have already been inspected; this is not preregistration for untouched data.
No new fits or downloads.

Use saved matched rotation predictions for ordinary matter, exact-third
companions, fitted simple MOND and restricted shared NFW. Keep all original
galaxies, radii and splits. Bin by R/Rdisk: below 1, 1 to below 3, at least 3;
central disk surface brightness: below 100, 100 to below 500, at least 500
solar luminosities/pc squared; atomic-gas proxy fraction: below 0.2, 0.2 to
below 0.5, at least 0.5. The proxy is 1.33 MHI/(1.33 MHI+0.5 L3.6), not a full
gas fraction: molecular gas and distinct bulge mass-to-light are omitted.

Average residuals within each galaxy/bin first, then weight galaxies equally.
Report predicted-minus-observed mean speeds, fractional errors, speed RMSE and
log-speed RMS, numbers of galaxies and points. Report descriptive 95% bootstrap
intervals for mean signed speed using 1000 galaxy draws and fixed seed 13092026.
These intervals describe scatter across this sample, not full measurement or
calibration uncertainty. Preserve separate train/validation/test summaries plus
a clearly labeled all-exposed-sample diagnostic. Compute paired outer-minus-inner
fractional residuals for galaxies represented in both radial bins.

Rank each galaxy's contribution to the equal-galaxy velocity and log-speed losses
without excluding any objects. For six lenses, compare saved full-covariance
all-motion chi-squared and signed velocity residuals for transferred companions
and the matched stellar-only control. Keep Chabrier and Salpeter proxies separate;
do not count them as independent galaxies or pool lens and rotation objectives.

Deliver one diagnostic report, source hashes and executable extraction. Determine
whether a uniformly stronger extra component is supported, or whether sign/shape
conflicts make such a change inappropriate. Do not introduce a new parameter in
this quick check. A promising residual association is not a derived interaction.
