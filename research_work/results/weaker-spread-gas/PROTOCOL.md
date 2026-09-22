# JR-4: weaker gravity, redistributed companion structure, and gas

21 September 2026, America/Los_Angeles. Baseline main 7386427ad7228b8d57319b1be78bbea3c47f31e8. Exploratory construction responding to the user's proposal; not independent confirmation. Declaration precedes the new fits. Original observations, R10 parameters, source masses, distances, inclinations, stellar populations and orbital parameters stay fixed.

## Theory first

The companion hypothesis permits environment-dependent strength and spatial organization. Test whether faint/gas-rich galaxies need a weaker total galactic response, a more extended companion state, or both. A changed effective galaxy-scale coupling is not automatically a change to laboratory G or to a cloud's self-gravity. The causal claim that weaker gravity leaves more gas is separate from reproducing rotation and lensing.

## Fixed mathematical alternatives

Write the R10 baryonic and companion inward radial accelerations as gb and gchi. Use

    Sigma = Mb/(2*pi*Re^2), in solar masses/pc^2
    w = Sigma0/(Sigma+Sigma0)
    s = exp(-a*w), lambda = exp(b*w)
    gchi_new(r) = s/lambda^2 * gchi_old(r/lambda).

Sigma uses the original R10-normalized stellar plus gas source mass and the fixed photometric half-light radius. It is a global compactness proxy, NOT a measured gas surface density. Gas fraction and observed velocity residuals are not inputs to the correction. Shared bounds: a,b in [0,ln(10)], log10 Sigma0 in [-1,3]. The exponent of the environmental gate is fixed to one. R10 source coefficients are not refitted.

The dilation is an exact mass-preserving redistribution when s=1: A_new=A/lambda, rc_new=lambda*rc, rt_new=lambda*rt; hence Mchi_total=A*rt/G is unchanged. s<1 reduces the effective gravitational normalization by s. That reduction is not an explanation of where stored energy goes; coupling/storage dynamics remain open.

Compare baseline, one constant overall suppression, environmentally weaker companion only, companion spreading only, weaker-and-spread companion, weaker total galactic response, and weaker total response plus companion spreading. In the total-response variants g_new=s*gb+s/lambda^2*gchi_old(r/lambda). Ordinary-matter distribution is not dilated. These are static effective constructions, not a completed three-dimensional modified-gravity field equation.

## Data and scoring

Use all 149 archived SPARC galaxies and 3152 positive-radius rows, preserving the 89 training / 29 validation / 31 comparison split. Fit universal constants on the 89 only, with equal-galaxy mean fractional squared velocity residuals as the primary development objective. Retain original velocity-error chi-square and physical-unit errors. A separate standardized-residual fit of the same declared forms is a sensitivity, not a license to select an objective from comparison results. All samples were examined earlier: none is newly blind.

Select the primary candidate by the 29-galaxy fractional-squared validation score; retain every form. Save fitted parameters and selection before scoring the 31 and the six lenses. Do not revise the equations or thresholds after their exposure. The lens predictions use the same new force in both the inherited annular stellar Jeans operator and photon deflection, with no new lens nuisance fit. R10 conditional distances, equal-potential assumption and spherical lens geometry remain limitations.

Run individual strength/spread fits separately as inverse diagnostics, allowing no suppression stronger than tenfold or dilation larger than tenfold. They are not universal predictions. Summarize the original 54 >20-percent outliers and the pre-existing overpredicted/underpredicted split; retain all 149. Correlation with gas fraction is descriptive and uses a mass model also entering Sigma, so it is not independent evidence of causality.

## Gas implications and what is not in these data

The archived inputs do not include a joint star-formation history, molecular-cloud density/velocity dispersion model, or gas outflow calculation. Do not claim to predict gas fractions or star-formation rates from a better rotation fit.

Calculate the distinct conditional limits explicitly. For a fixed-density isolated cloud, changing its self-gravity to Gcloud=sG gives t_ff proportional to s^(-1/2), Jeans length proportional to s^(-1/2), and Jeans mass proportional to s^(-3/2), at fixed sound speed. Confirm collapse times with a radial equation integration. In a fixed-shape galactic potential scaled by s, escape speed scales as sqrt(s), so weaker gravity also makes gas easier to lose. It does not create gas.

For a razor-thin gas disk with unchanged Sigma_g and sound speed, Q=cs*kappa/(pi*Gcloud*Sigma_g). If only the background rotation force is weakened by s, kappa scales as sqrt(s) and Q decreases. If both background and perturbing gas self-gravity scale by s, Q increases as s^(-1/2). Thus rotation-derived weakness alone does not establish slower cloud formation. The missing physical quantity is the perturbation response, not another fitted gas fraction. For redistributed profiles, calculate the companion potential-depth and radial-force ratios directly rather than apply a uniform-scaling identity incorrectly.

## Reproducibility and interpretation

Verify baseline predictions against packaged R10 arrays, positive monotone companion enclosed mass, self-similar mass normalization, independent lens quadrature, and refinement of selected lens predictions. Preserve input hashes, all starts, failures, and all per-object outputs. No cluster result, microscopic energy-budget closure, nested Solar-system solution, or causal star-formation validation is claimed. A failed declared form rejects that form under its assumptions, not the wider companion or environmental-gravity hypothesis.
