# OP-1: finite outer-support persistence
Declared 19 September 2026, before implementation and calculation.

## Question and scope
Can an ordinary-source-dependent radial retention of the existing CMF extra force improve outer galaxy rotation without losing cluster pressure or optical transfer? This executes the first priority of the companion-stream research brief. It is an effective radial closure, not a derived vector swirl, microscopic graviton interaction, or complete energy-funded theory.

No dark matter, expansion, changed G/c, distance fitting, or object-specific force strength. Reuse the frozen 149-galaxy / 3152-row, 12-cluster / 246-pressure-row and six reconstructed Coma-bin reduction. All observations were previously exposed; no blind-validation claim. Keep CMF's 11 shared coefficients fixed. Source luminosity/mass assumptions, cluster boundary pressure fitting and Coma geometry limitations are inherited and must be reported.

## Declared equations
Let g0 be the frozen CMF shared acceleration and gb ordinary acceleration.
q(r) = max[r (g0(r)-gb(r)), 0], in squared-speed units.
L(M) = L0 (M / 1e11 solar masses)^eta.
On increasing radii,
H[0] = q[0]
H[i] = max(q[i], H[i-1] (r[i-1]/r[i])^delta exp[-(r[i]-r[i-1])/L]).
W = 1 / (1 + (gb_SI / 1e-7)^2).
g_new = g0 + A W (H-q)/r.

The recurrence is a decaying running maximum. It stores a support envelope, not a conserved energy. The continuum discrete-source counterpart is max over earlier nodes j of q[j]*(r[j]/r[i])^delta*exp[-(r[i]-r[j])/L]. Both implementations will be checked independently. Flat squared-speed retention yields approximately 1/r extra force until finite decay; negative delta tests outward growth. This is not a proof of inward self-attraction or a lifetime inferred from data.

Evaluate on the union of supplied radii and a 512-point log-radius grid over their range; interpolate gb linearly in log radius. No extrapolation of missing inner galaxy sources. For A=0 reproduce the original baseline exactly. Source interpolation and unknown interior support are limitations, with declared sensitivity tests.

## Candidate list, fixed before results
delta in {-0.5, 0, 0.25, 0.5, 1}
L0 in {30, 100, 300} kpc
eta in {0, 1/3, 1/2}
A in {-0.5, 0.25, 0.5, 1, 2}
225 combinations plus the unchanged baseline. They are variations of one closure, not 225 independently invented theories. Negative A is a declared sign-reversal control. Reject nonfinite or negative total force, retaining reason and candidate. No continuous optimizer; the grid selects four universal settings in addition to the frozen 11 coefficients. No per-galaxy retuning.

## Ordered execution and selection
1. Save source/input/Git hashes, versions, inventory and controls.
2. Map baseline signed residuals across all historically exposed data, recording split, radius, fractional radius, acceleration and source mass. This is diagnostic; it must not amend the declared grid.
3. Evaluate every candidate on train/validation only. Use existing equal-object galaxy RMSE and pressure chi2 per point with one nonnegative boundary pressure per cluster, diagonal errors only.
4. Select 'protected': smallest validation galaxy RMSE among cases whose train AND validation pressure scores are no worse than 1.05 times their baseline scores; baseline is eligible. Break ties by training galaxy RMSE then candidate id.
5. Select 'joint': smallest validation [(galaxy_RMSE/20)^2+pressure_chi2_per_point/10], then candidate id, baseline eligible. Save selection before current-run test/Coma evaluation.
6. Evaluate frozen selections on the exposed test-labelled sample, baseline, matched simple MOND, and paired object bootstrap against each baseline. Count shared coefficients and scanned settings.
7. Coma: both frozen source brackets, unchanged same-force optical coupling, beta restricted to [0,1], fixed 100 Mpc distance, no redshift-to-distance expansion law. Primary reach 9000 kpc; 3000 and 30000 kpc are sensitivities, not reselected fits. Evaluate baseline and both selections.
8. Archive failures as well as successes and report all selections. No new hypothesis after scores without a separate dated amendment.

## Verification and decision rules
Controls: independent quadratic-cost envelope versus recurrence (relative 1e-11); A=0 baseline; zero source; zero correction on a constant q with positive decay; analytic retained-peak decay; large-acceleration release; same source distances/inventory; independent pressure integral versus matrix.
Selected cases: 512 versus 1024 internal grid at every galaxy (relative force change <=0.005 using per-object maximum force normalization); cluster source grid 600 versus 1200 (pressure <=0.005); combined refined Coma source/force/optics (shear <=0.005). Record failures, do not erase or relabel them. Inner-domain removal of the first 10 percent of galaxy rows is a sensitivity, not an accuracy check or revised fit. Include repeat-step no-double-counting for identical sources and dimensional/rescaling checks where applicable.
Descriptive success requires lower validation AND test galaxy RMSE than frozen CMF, train/validation/test pressure within 5% of CMF, and numerical checks passing. Beating MOND requires the analogous matched galaxy score plus paired bootstrap interval entirely below zero on reused data; this is still not independent confirmation. Do not claim lensing success from profile pressure, free beta>1, or an arbitrary optical multiplier. Report bounded Coma chi2 and preferred beta rather than invent an absolute acceptance criterion for incomplete data.

## Energy and attribution
The support envelope and its strength are hypotheses. H is not energy density; conserved test orbits in an imposed radial potential do not close source fuel, field production, carrier recoil or photon energy. If useful, its force requirement becomes a target for the separately proposed energy-counted vector model. If it fails, report which assumption fails without excluding all persistent swirls.
Running maxima, exponential relaxation, numerical integration, regression/selection and bootstrap are established mathematics. Existing CMF/IH reports credit SPARC, X-COP, Kubo et al., Milgrom and standard weak-field optics. The brief credits collective alignment to Vicsek et al.; this campaign does not implement or claim that mechanism. No historical novelty claim. Numerical verification is separate from empirical adequacy.
