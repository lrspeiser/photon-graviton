# Option 3: frozen exploratory protocol

Written before fitting or scoring this run. This is a computational galaxy holdout, not a historically blind discovery: SPARC and the RAR are established, public data, and this user's previous work has used SPARC. Earlier redshift fit is supplied by the user and is not revalidated here.

Use the original SPARC Rotmod_LTG.zip and SPARC_Lelli2016c.mrt downloaded from https://astroweb.case.edu/SPARC/. Treat published radii, distances, inclinations and baryonic mass models as inputs. No expansion or dark matter term enters the orbital calculation, but some adopted catalog distances are Hubble-flow based. This is conditional on catalog geometry, not a cosmology-independent test. Rotation tracers are principally gas, not individual stars.

Keep galaxies with quality Q<=2, inclination>=30 degrees, positive disk scale length, and at least five finite rows with positive radius, observed speed, speed error and baryonic radial acceleration. Preserve signed component v*abs(v), including negative gas contributions. Fixed disk mass-to-light ratio 0.5 and bulge 0.7 solar units; no per-galaxy adjustment. Do not use Vflat as input. Minimum point count and positivity are data-validity filters, not model-residual cuts.

Sort eligible galaxy names by SHA256('option3-v1:'+name), split floor(0.6*N) train, floor(0.2*N) validation, remainder test. All radii of each galaxy stay together. Fit shared coefficients on training only; choose candidate by validation loss; no train+validation refit. Save coefficients, chosen model and test predictions before scoring test outcomes. Model code may load the archive to prepare data, but test speeds are not passed to fitting/selection and are not printed or plotted before the freeze. Not an externally sealed blind experiment.

Set k=0.000077315 per million light-years; aK=c^2*k in SI; x=gbar/aK. Postulate gtotal=gbar+aK*f, dln(T)/dr=aK*f/c^2, vpred=sqrt(r*gtotal).

Candidates fixed in advance:
1. constant: f=A.
2. radial: f=A/(1+r/Rdisk).
3. square_root: f=A*sqrt(x).
4. power: f=A*x^p, p in [0,1].
Fit log10(A) in [-5,1]. Baseline ordinary: f=0. A reference benchmark, excluded from candidate selection, is the established RAR: gtotal=gbar/[1-exp(-sqrt(gbar/a0))], with a0=alpha*aK and log10(alpha) in [-5,1], trained identically. These are phenomenological functions, not derived field equations; the square-root candidate has the familiar deep-MOND scaling. Positive added acceleration cannot fix galaxies whose ordinary-matter prediction already exceeds measured speeds.

Training objective and primary score: sqrt(mean over galaxies(mean over radii(log10(vpred/vobs)^2))). This gives each galaxy equal weight; no use of an independent-point chi-squared. Also report galaxy-balanced RMSE in km/s and mean absolute fractional error. The objective intentionally measures predictive scatter, not a full observational-error likelihood. A low score is not evidence of fit within observational uncertainties.

After freeze, score all fixed models on test once. Use 2000 paired galaxy-bootstrap samples, fixed seed 20260908, for conditional 95% intervals on baseline-minus-selected and RAR-minus-selected primary score differences. No parameter refits in bootstrap. Report outer-third-of-radius-range scores and distance-method 2/3/5 test subset. Show all held-out curves and all galaxy errors to avoid cherry-picking. No post-test model changes or optional parameter searches.

Identifiability: with free A or alpha, fixed k does not test a numerical link to redshift. Rescaling k can be absorbed by the coupling (for power, A transforms as k^(p-1)). No photon transport, lensing, vertical motions, conservation law or solar-system behavior is established. A static local clock field produces endpoint clock shifts, not an automatic path-accumulating cosmological redshift. The two mechanisms still require a common dynamical derivation.
