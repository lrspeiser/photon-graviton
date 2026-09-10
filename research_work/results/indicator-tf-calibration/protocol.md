# Indicator-anchored i-band calibration protocol

Declared before this calculation: use the staged 109 exact indicator/measurement overlaps, without querying or fitting recession velocities. This is an exploratory calibration, not a fresh holdout or a physical void map. The previous goal response was a clarification, not new research evidence; this calculation resumes stage 2.

Use finite icmag, Ai, Wmx, e_Wmx, Inc and e_Inc; require Wmx>0, e_Wmx>0, 45<=Inc<=90, e_Inc>0 and e_Wmx/Wmx<=0.2. Record all excluded matches. No luminosity, redshift, residual or inferred-distance rejection. These initial quality cuts are analyst choices, not a claim of survey completeness.

Known empirical Tully-Fisher form: x=log10(Wmx/sin Inc)-2.5; M=a*x+b. Corrected apparent i magnitude m=icmag-Ai follows the source ReadMe. Compare two declared photometry scenarios: p=0 as a conventional inverse-square diagnostic and p=1 for the proposed shared frequency/event stretch. With fixed prior alpha=0.0002488993286382367/Mpc, fit y=m-[5log10(D)+25]-5*p*alpha*D/ln(10). The latter term assumes static dilution, conserved photons and S=exp(alpha D); it is conditional transport photometry, not a first-principles derivation. Stipulate the published indicator D for this initial calibration; its own brightness/zero-point dependencies remain unresolved.

Fit ordinary least squares with no clipping, and leave out entire CF4 identifier groups for exploratory cross-validation. Read only PGC/group columns from the cached CF4 file. Group membership may itself use velocity information; here it only prevents some shared-group leakage, never supplies distance or predictor values. Shared indicator zero points still cross folds. Also report each indicator method's residuals. Do not select the better scenario or retune alpha.

Invert each held-out predicted apparent modulus by solving D*exp(p*alpha*D)=10^[(m-a*x-b-25)/5]. Compare against the stipulated indicator distance. Report predictive scatter and multiplicative distance error, not a full errors-in-variables likelihood. Propagate linewidth/inclination uncertainties for a diagnostic only; no complete dust/selection/calibration covariance is available in the staged subset. OLS slope attenuation, common zero points, selection and correction dependencies prevent precision or unbiased-distance claims.

Independent verification: bracketed root versus Lambert W inversion, full rank of each fit, finite values and one and only one excluded-group prediction per eligible galaxy. Publish all rows and both scenarios. Do not deploy the estimator as a void map until its calibration quality and selection are adequate.

Source conventions: staged J/ApJ/902/145 ReadMe notes 2, 4, 6 and 7; https://arxiv.org/abs/2009.00733 . K corrections use measured spectral shifts and are not automatically expansion-derived distances, but their inherited implementations and dust estimates still need joint modeling.
