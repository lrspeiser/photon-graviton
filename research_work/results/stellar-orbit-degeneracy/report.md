# The wave advantage survives three shared orbital assumptions

The frozen wave model has a lower combined error than the common-stellar-mass benchmark under all three tested orbital assumptions. This narrows one uncertainty: the comparison is not reversed by changing the stars from modestly tangential to modestly radial orbital preferences.

This is a **limited sensitivity result**, not a determination of actual stellar orbits, a new validation claim or full model selection. The ordinary-mass multiplier and wave parameters stay fixed. The tested alternatives do not span every allowed orbit distribution, and neither model is jointly refitted with population-mass uncertainty.

## What changes physically

We use the known constant-anisotropy parameter

\[
\beta=1-\frac{\sigma_\theta^2+\sigma_\phi^2}{2\sigma_r^2}.
\]

Beta=-0.3 represents a modest preference for sideways/tangential motion, beta=0 isotropic velocity dispersions, and beta=+0.3 a modest radial preference. Each value is shared across galaxies within a comparison. The values are declared sensitivity choices, not observed priors or a new orbital law.

The same mass can produce different measured line-of-sight speeds when stellar orbits favor different directions. In this calculation, changing the luminous tracer's orbital distribution leaves the source density and lensing unchanged. It is distinct from changing the hypothetical companion reservoir's own support mechanism.

The **known spherical Jeans solution** is

\[
j(r)\sigma_r^2(r)=r^{-2\beta}\int_r^\infty j(s)g(s)s^{2\beta}ds.
\]

The projected velocity variance includes the known factor 1-beta R^2/r^2. The aperture integration implements this through

\[
W_\beta(r)=\int P(r\sin\theta)(1-\beta\sin^2\theta)\sin\theta\,d\theta,
\]

\[
K_\beta(s)=s^{2\beta}\int_0^s r^{2-2\beta}W_\beta(r)dr.
\]

P is the Gaussian-seeing aperture-inclusion probability. The projected variance numerator is the integral of j(s)g(s)K_beta(s); the denominator uses the ordinary luminosity weight W_0. These are standard moment equations and a reordered integration, not a new photon-companion formula. The same expressions were derived and checked in the earlier [constant-beta lens analysis](../lens-training-pilot/anisotropy-report.md); here they are applied to the newly supported wave profiles.

The second verification route integrates the radial Jeans pressure first and then projects it. That checks the reordered-kernel implementation against an independent integration order for the actual source profiles.

Across all 234 source/orbit calculations, the two variance calculations agree within 7.4e-6 relative. The beta=0 dispersions reproduce the archived isotropic predictions exactly, and lens angles are identical across beta. These numerical checks do not establish a physical orbital population.

## Shared assumptions and preserved calibration

The wave model retains m=10^-24 eV/c^2 and source fraction f=0.6. The ordinary benchmark retains its previously training-fitted common stellar-mass factor lambda=1.6453621586. Both use the same underlying rescaled Salpeter populations, sizes, static geometry, aperture and seeing assumptions.

For the mass-corrected ordinary benchmark, its constant-beta dispersion scales as sqrt(lambda) at fixed density shape. Its lens prediction is taken from the correctly recomputed shared-mass benchmark; no linear scaling of Einstein angle is assumed.

No mass, source fraction or beta is selected anew on the reused validation sample. All three alternatives are reported on training and the already-exposed validation galaxies as sensitivity results. No final-test predictions are opened. In particular, the original isotropic calibration remains the frozen prediction; the somewhat better radial result is not silently substituted for it.

## Results

| Sample | Beta | Ordinary mass-correction dispersion RMS, km/s | Wave dispersion RMS, km/s | Ordinary joint score | Wave joint score |
|---|---:|---:|---:|---:|---:|
| Training | -0.3 | 31.56 | 33.83 | 0.02443 | 0.01578 |
| Training | 0 | 31.04 | 31.99 | 0.02430 | 0.01467 |
| Training | +0.3 | 31.36 | 30.68 | 0.02462 | 0.01387 |
| Reused validation | -0.3 | 50.52 | 42.02 | 0.06633 | 0.04838 |
| Reused validation | 0 | 45.85 | 38.22 | 0.06287 | 0.04495 |
| Reused validation | +0.3 | 40.04 | 34.81 | 0.05921 | 0.04154 |

The joint score is the previously declared equal-weight mean squared logarithmic error for dispersion and lens angle. It is not a likelihood with population, orbit and imaging uncertainties.

Lensing RMS remains 0.2278 arcseconds for the ordinary mass-correction benchmark and 0.1440 for waves on training; on reused validation it remains 0.3541 and 0.2699 respectively. Those values are exactly unchanged across beta in the stored predictions, as required when the gravitating density is unchanged.

The ordinary benchmark is slightly better on training dispersion at beta=-0.3 and 0. The wave model has a better joint score at every tested beta because of its lower lensing error, and lower dispersion RMS at each corresponding beta in reused validation. Thus the result is more qualified than saying waves win every measurement under every assumption.

## What this does not settle

The stellar-mass factor was originally fitted under isotropy and is not reoptimized at each beta here. Likewise the wave mass and source fraction were selected under isotropy and stay fixed. This is a comparison of sensitivity around frozen choices, not the best possible ordinary model versus the best possible wave model or a fair-evidence calculation with fully matched nuisance priors.

Positive Jeans moments alone do not prove a globally nonnegative, dynamically stable stellar distribution function for every new potential and anisotropy. The prior analytic tracer proof applied to the empirical force family; it should not be claimed automatically for these wave profiles. A distribution-function or orbit-population construction remains needed before promoting a nonzero-beta alternative to a completed predictive model. This is separate from stability of the wave source itself.

Spatially varying anisotropy, flattened/triaxial structure, distinct stellar populations and population-mass gradients can change the result further. They need observational constraints and comparable treatment in both hypotheses. Testing three constants does not marginalize over those possibilities.

The next decisive evidence remains a better constrained joint mass/orbit analysis and a frozen final-test comparison, together with the missing photon-production and capture physics. Neither the orbital sensitivity nor improved galaxy scores derives a photon-to-wave conversion rate, explains supernova event stretching or passes the deferred energy-supply budget.

## Reproduction

Run `python research_work/results/wave-lens-training/run.py --grid-index 4 --orbit-audit` and the same runner with `--role validation --orbit-audit`, then `python research_work/results/stellar-orbit-degeneracy/summarize.py`. The new orbital outputs are saved separately from all original isotropic predictions. `results.json` and `predictions.json` here contain the paired hypotheses; input hashes and the independent-pressure check are preserved. Dependencies are NumPy and SciPy.
