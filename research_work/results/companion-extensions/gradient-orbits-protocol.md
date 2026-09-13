# Combined stellar gradient and radial orbital fit

13 September 2026. Protocol written before executing the combined fit.

Use the same six exposed lens systems, conditional nonexpanding geometry, fixed light profiles, PSF, covariance, adjusted-reference companion amplitude and exact-third response as the preceding separate experiments. At each stellar gradient, use the catalogue lens angle to fix total stellar mass. This consumes lensing calibration; it is not a predicted lens angle.

Combine Upsilon(r)=Upsilon_out[1+h/(1+(r/Re)^2)] with beta(r)=beta0+(beta_infinity-beta0)r^2/(r^2+Re^2). Keep Re fixed to the measured projected half-light scale. Gradient h lies in [-0.8,9], central beta in [-2,0.45], and outer beta in [-2,0.95]. The expanded outer bound was selected in prior exposed-data work and is not an independently specified physical limit. No new companion parameters are fitted.

Fit all three nuisance parameters to inner motion bins only. Evaluate the covariance-conditioned outer-bin residual without fitting it. Use direct Jeans integration and multiple starting points, including both previous separate-fit solutions. Require reproduction of the h=0 radial-orbit and equal-endpoint gradient limits. Require positive predicted second moments and positive stellar mass/density. Verify the imposed angle algebra. Report all cases, boundaries, convergence and numerical limitations. A lower inner residual is expected from added freedom; the distinguishing outcome is whether outer motions improve without unphysical or unbounded parameters. Do not label a necessary orbital criterion as sufficient, and do not carry the earlier orbital positivity check over to newly fitted parameters.

Compare with both separate extensions and the free-mass motion control. No branch replaces the reference from this exposed-sample diagnostic. Distribution-function positivity, stability, stellar-population evidence and nonspherical modeling remain separate requirements.
