# Brightness-distance consistency on exposed data

2026-09-10. The preceding goal turn made progress by separating chosen rate/law assumptions from the messenger-delay result and recording the photon-generated-wave source branch. Complete the already planned link between stage 3 photometry and stage 2 independently estimated distances.

Use every row of redshift_paper/all_164_groups.csv and the existing alpha, unchanged. No data acquisition, holdout opening, rate fitting, galaxy-specific motion fitting or alteration of the primary stipulated catalog distances. All original train/validation/test rows are previously exposed; retain their labels only for historical diagnostics.

Conditional bolometric indicator sensitivity: conserved photons, spectral stretch S=exp(alpha D), duration stretch S^b, static geometry, fixed populations, ordinary endpoint standards. Then D_app=D*exp(p*alpha*(D-D_anchor)), with p=(1+b)/2. Use b=1 (shared spectral/event transport) and b=0 (the timing-deficient energy-only control). Evaluate declared calibration-anchor distances 0,1,10 Mpc in each; zero is the ideal unaffected-calibration limit, not an actual calibrator. Retain all six scenarios and baseline p=0. Do not choose an anchor by outcome agreement.

Use known Lambert W inversion and independently verify with bracketed roots. Propagate the published modulus errors only as conditional local derivatives, verifying against finite differences; do not turn them into an incomplete chi-square or a claimed full covariance. Compute c*(predicted z-observed z) RMS, bias and MAE for all rows and original partitions, retaining negative residuals. Summarize the distance and predicted-redshift changes.

These are idealized sensitivity scenarios, not the actual SBF correction pipeline: passbands, source colors, K corrections, Doppler/endpoint effects, standardization and shared calibration are not modeled. Observed redshift is never used to invert a distance. Alpha was nevertheless calibrated from this same exposed sample, so this is not independent validation or a refitted joint photometric model. Known algebra and empirical calibration are not novel physical laws.
