# Direct conversion branch: fixed comparison protocol

9 September 2026, before running this new comparison. User direction: choose a coherent branch without a sign-changing time prescription or screening of conversion, integrate it with companion transport and bound storage, and test redshift at 100 million light-years. The existing data and partitions are already exposed; this is not fresh validation.

## Branch decision

Use direct, continuous fractional photon-energy transfer per traveled distance. Set alpha=alpha_0 constant for this first branch; no extra clock multiplier, no time-field source and no environmental screening of conversion. Short distances produce small accumulated effects, not an onset threshold. Retain ordinary local clocks/c and baseline gravity as the initial coupling, forward light-speed companions without further loss, no permanent storage in voids, and the proposed deep-well capture into extended bound modes. Capture's environmental dependence is distinct from screening photon conversion. Capture microphysics, stable modes and joint gravitational response remain unresolved. Energy supply is deferred, not solved.

## Data and frozen choices

- Source: `redshift_paper/all_164_groups.csv`, required SHA-256 8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0. Use only identifiers, split, sky coordinates, adopted distance and observed CMB-frame cz/z. Do not import its archived model-prediction columns.
- Preserve every distance and observed spectral shift; no clipping, per-object velocity adjustments, fitted distances, new environmental terms, or new cosmological prior estimates.
- The old alpha=0.0002488993286382367/Mpc and old linear rate are reused for the historical train/validation/test illustration. Each comparison model fits one scale parameter on the 104 old training rows; retain all old partitions. Verify the direct-conversion fit reproduces the saved rate before reusing it.
- At 100 million light-years use D=100e6*c*(31557600 s)/(1e6*648000/pi*149597870700 m), with c=299792458 m/s. Report redshift, wavelength stretch, energy fraction transferred, and unchanged-event-duration prediction. Illustrate an emitted 500-nm wavelength. This is the branch's geometric-path prediction, not an actual star observation.
- Select the five recovered test-partition rows nearest to that distance using distance only, with PGC integer as tie breaker, before calculating residuals. These are group representatives, not confirmed named stars or necessarily isolated galaxies. Show all five and their measured versus predicted shifts; no nearest-fit selection.

## Expansion comparison, for evaluation only

Expansion is an external comparison, not the active explanation for the fictional universe. Use spatially flat FLRW with Omega_m=0.3, Omega_Lambda=0.7, radiation neglected and H0 as the only fitted parameter. These fixed shape parameters define an illustrative comparison, not a comprehensive best-fit cosmology. Define E(z)=sqrt(0.3(1+z)^3+0.7) and I(z)=integral_0^z dz'/E(z'). Use established distance formulas from [Hogg](https://arxiv.org/abs/astro-ph/9905116).

Distance meaning matters. Evaluate two explicitly separate comparisons; do not select the favorable one afterward:

1. **Geometric-distance convention.** Active conversion uses z=expm1(kD/c). FLRW comparison uses D=(c/H0)I(z), its present-day line-of-sight comoving distance. In the static branch D is path length; in the expanding model it is not photon path length. These are conditional interpretations of the same supplied numbers, not proof that both are the original catalog's distance definition.
2. **Luminosity-distance convention sensitivity.** FLRW uses D_L=(c/H0)(1+z)I(z). For the static number-retaining, no-event-stretch branch, bolometric flux implies D_L=r*sqrt(1+z), with 1+z=exp(alpha*r). Therefore z=expm1(2*LambertW(alpha*D_L/2)). This is a separate interpretation test, not a change to any adopted distance or the branch's geometric prediction. Actual surface-brightness-fluctuation calibration must be forward-modeled before treating this generic bolometric sensitivity as a faithful reinterpretation of that catalog.

Fit scales k=c*alpha or H0 in [0,150] km/s/Mpc by squared cz residuals, with bounded scalar optimization and explicit boundary evaluation. Include a one-parameter linear small-distance control, but do not equate it with a full expansion model.

Also perform leave-one-coarse-sky-region-out predictions for all five models, using the existing fixed rule: RA quadrants plus four if declination>=0. Every fit uses only the other regions. Report pooled RMS, MAE, bias and each region's size. Use 2000 paired resamples of whole held-out regions, seed 2026090923, to describe RMS differences for each matched-distance comparison. These are conditional fixed-prediction bootstrap intervals, not new holdouts or complete fit/selection uncertainties.

## Numerical verification and reporting

Use 32-point Gauss-Legendre integration and iterative inversion for the expanding models, then independently verify every historical all-row prediction at its fitted scale with adaptive quadrature and bracketing roots (maximum absolute z difference below 1e-11). Verify the luminosity-conversion inverse by substituting its redshift into the original distance relation (relative error below 1e-11). Reproduce the saved conversion scale within 1e-5 km/s/Mpc and the saved test RMS within 1e-5 km/s. Preserve the data hash. No claimed superiority unless the declared comparison actually supports it, and even a numerical improvement remains exploratory.

## Whole-system requirement

Write the source -> photon loss -> traveling companion -> bound store -> gravity chain and an explicit status for each link. A scalar photon-loss equation does not specify a coherent wave, exact achromaticity, momentum exchange, capture or an extended halo. The stationary fixed-speed branch predicts event-duration stretch 1, not automatically 1+z; flag this against supernova timing observations rather than treating them as optional. Relevant primary evidence: [DES supernova time-dilation study](https://arxiv.org/abs/2406.05050). No new timing data analysis or energy-supply calculation is part of this run.
