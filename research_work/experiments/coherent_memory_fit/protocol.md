# CMF-1: matched comparisons, outward profile memory and bounded lens geometry

Declared before calculation, baseline main 62d81d3. This is the next finite campaign toward the user's shared stellar-motion and cluster-lensing solution. Success is not assumed from the instruction to pursue it. Existing archives remain unchanged.

## Requirements and data
Use IH-1's 149 SPARC objects and all 3,152 positive-radius rows, 89/29/31 object split, fixed catalog distances, disk/bulge mass-to-light 0.5/0.7, and all signed gas contributions as already reduced by data.py. Retain zero ordinary-acceleration rows. Use its 12 X-COP gas/stellar sources and 246 pressure points, 6/3/3 split and one nonnegative outer pressure nuisance per cluster. No inferred halo/hydrostatic-mass column is evaluated. Use the same six exposed Coma figure bins and two fixed source brackets.

No dark matter, expansion, distance fitting, changed G/c, object-specific response coefficients or target-derived input features. All data are historically exposed: these partitions are not a new blind test. R500 and other archived reduction conventions remain conditional input assumptions.

## Candidate equations
Let x=tanh[-ln(max(g_b,1e-30 SI)/1e-10)/4], y=tanh[ln(r/10 kpc)/4], z=tanh[ln(M/1e11 Msun)/4].
For logarithmic radius t=ln(r), outward memory satisfies ell dH/dt=x-H, starting H=x at the innermost supplied ordinary-source radius. Integrate exactly for piecewise-linear x(t). h=H-x. ell is dimensionless radial persistence, NOT time, a propagation speed, photon binding probability or a vector swirl solution. Boundary dependence must be reported.

With P=sum(theta_j F_j), propose
    D=1/[1+(g_b/1e-7 SI)^2]
    g=g_b exp[8 D tanh(P/8)].
This gives positive total force for g_b>0, permits enhancement and suppression, and recovers g/g_b -> 1 at high acceleration. At g_b=0 it predicts zero; retain those residuals. Fixed exponent cap 8 is a numerical/physical hypothesis, not a measured constant. No claim that this force is funded by an emission/binding action.

Ordered bases:
L: 1,x,x^2,x^3
R: L plus y,x*y,y^2
M: L plus z,x*z,z^2
RM: R plus z,x*z,z^2,y*z
RMH: RM plus h,x*h,z*h, with ell=0.3,1,3.
These seven structural variants have 4/7/7/11/14 coefficients. Coefficients lie in [-8,8], all signs permitted. Ridge strengths lambda=0,0.001,0.01 multiply mean(theta^2). Joint training weights w=0.1,1,10 multiply cluster chi2/point/10 in the IH-1 galaxy equal-object MSE/20^2 objective. Also run galaxy-only weight 0: 7*3*4=84 fits with three deterministic-seeded starts each (252 attempts). Do not call these 252 new theories. Preserve all attempts and solver termination flags. Maximum 400 function evaluations per attempt.

Primary joint selection minimizes the fixed validation score G_RMSE^2/20^2+C_chi2_per_point/10 among optimizer-success candidates with w>0. A secondary gated selection, if available, additionally requires G_validation <= matched MOND and C_validation <=10; it remains separately labelled. Galaxy-only selection minimizes validation RMSE. Save selection before test-labelled scores and Coma. Never select using Coma.

## Matched controls and interpretation
Refit the standard simple-MOND reference g=(g_b+sqrt(g_b^2+4*a0*g_b))/2 to the very same training rows/object weights; one universal a0 in [1e-12,1e-8] SI. Attribute it to MOND, not this project. Evaluate IH-1's frozen released law and ordinary matter on identical arrays. MOND is a comparison, not an adopted mechanism. No active dark-matter benchmark is run.

Report validation/test RMSE and per-object paired differences. Use 10,000 paired object-bootstrap resamples of the already-exposed test partition, seed 200926; confidence intervals are conditional uncertainty diagnostics, not a correction for model search or prior exposure. Report parameter counts, training/validation gap, cluster residuals and all failures. A lower point score alone is not superiority.

## Light and geometry
Keep IH-1's stipulated shared weak-field light coupling with no separate photon strength:
    alpha(b)=4/c^2 integral_0^infinity g(sqrt(b^2+u^2))*b/sqrt(b^2+u^2) du.
Apply the response to 9 Mpc then continue the extra force as r^-2; 3/30 Mpc are declared reach sensitivities. Compute outward memory on the ordinary-source radial profile, not on the photon's observed deflection. Source radial grid starts 0.1 kpc; vary to 0.01/1 kpc for the primary chosen law. Preserve IH-1's r=h^-1 radius/0.7 convention.

In static Euclidean thin-lens geometry, gamma_t=(D_l*beta/2)*(alpha/b-dalpha/db), beta=1-D_l/D_s in [0,1] for background sources. Adopt the published Coma distance D_l=100 Mpc (Alabi et al. 2020, section 1) as a fixed input under the universe contract. Do not import that paper's expansion model or convert source redshifts into distances. The archived radial conversion and adopted distance are an approximate cross-catalog geometry, not a newly calibrated angular data reduction. Report this limitation.

Calculate fixed beta=0.5,0.9,1 and a physically bounded beta fit in [0,1], plus the formerly unbounded shape normalization for diagnosis. beta=1 is an upper amplitude envelope when shear is positive; it is not a measured source geometry or a universal pointwise envelope for signed shear. A best unbounded beta>1 cannot be repaired by moving sources in this geometry. Missing source weights, bin edges, covariance and angular calibration still preclude an absolute observational validation claim. No separate fitted lens amplification is permitted.

## Verification and completion
Before the campaign run, test analytic Jacobians against finite differences; pressure matrix against direct trapezoids; exact memory evolution against an independent ODE solver; unchanged distances/row counts; point-mass optics and high-acceleration limits. After selection, refine cluster radial grids 600->1200 and light quadrature/derivative, and repeat memory-grid/start-radius diagnostics. Save results even on failure and do not hide nonconverged fits.

Operational target: matched MOND galaxy validation/test point scores improved and validation cluster pressure <=10, plus physically admissible Coma amplitude consistent with its uncertainties, with no separate fitted optical strength. A bootstrap interval spanning zero is inconclusive. Full solution additionally requires new independent observations, source/field energy accounting, source generation, reciprocal matter response and a derived photon attachment/light law; numerical fitting cannot mark those complete.

Commit this declaration, commit implementation before fitting, checkpoint evidence, then commit the report/audit and push main. New findings may motivate a separately declared follow-up; do not keep changing this campaign until a preferred verdict appears.

## Attribution and research inspected
- SPARC: https://arxiv.org/abs/1606.09251
- MOND and its historical acceleration-law limit: https://adsabs.harvard.edu/pdf/1983ApJ...270..365M
- X-COP pressure profiles: https://arxiv.org/abs/1805.00042
- Coma shear: https://arxiv.org/abs/0709.0506 ; section 3.1 obtains its critical density using a cosmological redshift-distance mapping, which is not adopted here.
- Fixed Coma distance convention: https://academic.oup.com/mnras/article/496/3/3182/5859958
- Established line-of-sight optics: https://arxiv.org/abs/astro-ph/9912508
- Earlier symbolic acceleration-law searches: https://doi.org/10.1093/mnras/stad597

Polynomial response fitting, exponential relaxation, ridge regularization, bootstrap and lens integration are established tools. No claim of historical originality for them, or that this scalar profile memory constitutes a newly derived graviton swirl, is made.
