# Projected deposition field and the point-source limitation

## Outcome

A spatially resolved evaluator now sums the unchanged 814 THEMIS catalog sources and projects deposited power along the Earth-to-center direction through the stipulated 1 Mpc receiver. The map contains 1,264 point samples. It retains source directions rather than substituting a spherical density profile.

The calculation exposes a real limitation of treating galaxies as point emitters: the projected density directly over an internal emitter diverges logarithmically as the integration resolves smaller distances. A finite quoted peak from this model would be resolution-dependent. No lensing peak or observed lensing agreement is claimed.

## Equations: established transport, conditional interaction

For source luminosity L, source-to-point distance ell, pre-entry distance e and in-receiver distance s=ell-e, the incoming energy fractions are P0=exp(-alpha e), C0=1-P0. Internal emitters have e=0. Known coupled linear transfer gives

    C(s) = C0 exp(-kappa s)
           + P0 alpha [exp(-alpha s)-exp(-kappa s)]/(kappa-alpha).

Known geometric dilution then gives deposited power density

    q(x) = sum_sources kappa L C(s)/(4 pi ell^2).

The projected power is integral q(x_perp,z) dz through the sphere. These are conditional consequences of the optional constant-rate conversion/capture postulates, not novel fundamental laws. Luminosities are SED-model inferences and the archived alpha is empirical. The receiver, steady source histories, absence of intervening absorption, fixed geometry, point emission and permanent stationary storage retain prior limitations. No gravitational response has been applied.

For a point source inside the sphere, C(s)=alpha s+O(s^2) near the source. Consequently q approximately L kappa alpha/(4 pi s). Its three-dimensional volume integral is finite because the volume element contains s^2 ds. Its line-of-sight integral at exactly zero projected separation behaves as integral dz/abs(z) and diverges. This is a mathematical artifact of the point-emitter assumption, not infinite total deposited energy. Finite source extent or a specified finite observational aperture regularizes the projection; neither may be chosen to force a desired lensing peak.

## Resolution checks

The point grid uses 40 coordinates per axis within the circular receiver projection, avoiding exact center sampling. Increasing line-of-sight quadrature from 64 to 128 nodes changes values at points more than 0.05 Mpc from all internal projected source positions by at most 4.79e-5 relative (1,222 samples). Across all samples the maximum change is 1.252%. The 0.05 Mpc label is only a numerical diagnostic; it does not remove physical source regions or authorize their exclusion from a later fit.

A separate exact-center diagnostic deliberately probes the cusp:

| Line-of-sight nodes | Projected power at center (Lsun/Mpc^2) |
|---|---:|
| 32 | 2.25536e+08 |
| 64 | 2.6626e+08 |
| 128 | 3.07236e+08 |
| 256 | 3.48398e+08 |

The increasing values are consistent with the analytic logarithmic divergence; no convergence assertion is made there. Off-source agreement is numerical sensitivity, not catalog uncertainty or a proof of global integration accuracy. All evaluated values are finite/nonnegative at the selected quadrature points. This map is point-sampled, not pixel-averaged; summing samples is not an independent total-power audit. The preceding ray-integrated ledger remains the reference total.

## Consequence and next step

The spatial map advances the source calculation, but finite measured emitter sizes or resolved luminosity profiles are needed to predict central surface structure. A specified observing aperture can support a finite aperture-integrated test even before arbitrarily fine structure is resolved. Neither a deposit centroid nor a sampled surface-power maximum is automatically a lensing peak. Accumulation histories, supported reservoir motion, both metric potentials and forward image distortions remain to be calculated. No final astronomical holdouts were used; all six physical objectives remain open.

Data: map.json records coordinates, power and resolution diagnostics; results.json records the summary and source hash. Run `python research_work/results/cluster-projected-deposition/run.py`.
