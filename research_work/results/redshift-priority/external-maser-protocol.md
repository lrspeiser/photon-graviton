# Maser external comparison protocol, before table extraction

Source fixed to Pesce et al. 2020, arXiv:2001.09213v2, six-galaxy compilation. Its abstract and aggregate result have already been seen. This is an external diagnostic, not a pristine blind sample. No target selection based on residuals; retain all six and separately flag overlap. No target-level measurements from this source have yet been extracted in this pass.

Use the paper's geometric distances, not CF4 maser moduli registered onto a common MCMC distance scale. Read velocity conventions and frames before interpreting numbers. Use uncorrected CMB-frame optical c*z if supplied, or document the exact supported conversion; do not import a cosmology-derived flow correction. If conventions cannot be established, do not score that row until resolved.

Freeze existing rates: constant exponential alpha=0.0002488993286382367/Mpc, linear c*alpha=75.17651383428588 km/s/Mpc, c=299792.458 km/s. No refit, no new candidate and no per-object motion. Report all signed residuals, RMS, MAE, median absolute error, bias and actual sampled range. Report results separately within 10.2-93.2 Mpc and outside, where available; do not characterize extrapolation as prior validation.

Distances remain stipulated fixed values for the primary comparison. Published distance uncertainties are provenance/sensitivity information, not permission to move each target to the curve. The source's 250 km/s motion uncertainty is an assumption; an optional interval using it must be labeled illustrative, never certified predictive coverage or a success threshold. No statistically justified absolute-accuracy acceptance criterion is available yet. This diagnostic cannot complete the fresh-data or uncertainty requirements.

Before scoring: save source hash, extracted features and frozen predictions separately from observed labels; record all six identities and overlap with the repository exposure registry. Quarantine here means preventing further tuning after this protocol, not retroactively erasing the abstract exposure. A genuinely unexposed validation sample remains a separate requirement.
