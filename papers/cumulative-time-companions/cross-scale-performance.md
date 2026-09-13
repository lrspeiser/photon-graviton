# Cross-scale performance of radiation-conditioned companion deposition

Working research assessment, 13 September 2026. This is an evolving supplement to the historical manuscript, not a peer-reviewed result or a rebuilt PDF. Numerical claims below refer to the linked executable records. Results are conditional on a hypothetical nonexpanding universe and the explicitly stated geometry and matter inputs.

## Abstract

We investigate a phenomenological proposal in which photons transfer energy to traveling companions and some companion energy is retained in gravitational reservoirs. A bounded retention prescription with exponent one-third improves galaxy rotation predictions relative to an ordinary-matter-only baseline. On 60 comparison galaxies, it improves absolute speed errors for 56 objects but outperforms the tested simple MOND prescription for only 20. Its aggregate rotation error remains larger than MOND's. Two freely normalized companion density shapes describe the six reconstructed Coma shear measurements comparably to NFW; this does not test the galaxy-calibrated retention normalization. A new transfer of the exact-one-third galaxy parameters into six lens galaxies leaves approximately 14% lens-angle RMS discrepancies. These results support continued investigation of an empirical extra-gravity prescription, but do not establish a unified photon-to-gravity mechanism or superiority to MOND or dark matter.

## Model definition and provenance

### New extension assessment

The [companion-extension study](../../research_work/results/companion-extensions/report.md) retains the exact-one-third inventory while testing two distinct modifications. A 60-case axisymmetric gravitational-permittivity sweep gives no selected improvement in 18 outer Milky Way bins after choosing parameters on 20 inner bins. A 2,706-case ideal-Bose phase-selector sweep reduces inner RMS from 3.68/11.83 to 1.81/2.85 km/s for two ordinary-matter baselines; outer RMS changes from 8.36/9.78 to 8.36/9.29. Comparisons use the same spherical gas-monopole approximation, which differs from the earlier full-gas-disk transfer. Data were previously inspected; these are conditional spatial-transfer diagnostics, not blind or cross-galaxy validation.

The phase selector combines known ideal-gas statistics with a proposed companion kinetic scale and conservative contraction map. Its preferred effective masses (17.78/13.34 eV/c²) and compact fractions (13.9/32.3%) are model-dependent fit parameters. Transferring the two parameter sets between ordinary-matter baselines exposes substantial sensitivity. No thermalization, massive bound-state formation, feedback equilibrium, support or full energy balance is derived. The original light-speed traveling channel would need a distinct capture transition into such massive states.

An analytic Poisson-loss diagnostic shows that discrete fractional energy transfers broaden spectral lines unless losses are sufficiently small and numerous. For the archived redshift coefficient and a 100-million-light-year path, an illustrative added relative linewidth of 10^-5 requires per-event loss below 1.31e-8 and at least about 582,000 events. This is a design constraint, not a measured exclusion. Ordinary stationary whole-photon mixing changes intensity rather than redshifting surviving photons. Scalar/axion and string-inspired frameworks remain possible sources of explicit interactions, not derivations of the empirical law. No additional lensing success is claimed; the weak-field force pilot does not determine the relativistic lensing potential.

### Retained empirical foundation

The proposed physical interpretation is a lossless traveling companion channel with capture and retention. Conservation requires accounting for photon, traveling-companion, deposited, and released energy, together with boundary flux and any work done. A closed-volume reservoir equation without flux cannot be applied unaltered to companions traversing a cluster. Retention is not a source of new energy. Ordinary matter continues to contribute gravitational mass; these galaxy fits concern additional gravity, not replacement of the entire gravitational field by newly supplied photons.

The current galaxy retention prescription is

\[
X=\frac{L_{3.6}/(10^9 L_\odot)}{(R_d/\mathrm{kpc})^2},\qquad
\eta(X)=\frac{X^{1/3}}{1+X^{1/3}}.
\]

**Status:** the Hill/logistic functional form is established mathematics. The choice of this luminosity-density proxy and its interpretation as companion retention are project hypotheses. The exponent one-third is a simplifying candidate close to a fitted value 0.33999; it has not been derived as a fundamental constant.

The spherical effective capture prescription is

\[
a=s R_d,\qquad \kappa(r)=\frac{k_0}{[1+(r/a)^2]^2},\qquad
J(r)=\frac12\int_{-1}^{1}e^{-\tau(r,\mu)}d\mu,
\]
\[
\tau(r,\mu)=\int_{\mathrm{incoming\ ray}}\kappa\,d\ell,
\qquad
\rho_{\rm dep}(r)=\frac{2 C_{1/2}\eta(X)J(r)}{[1+(r/a)^2]^2}.
\]

**Status:** ray attenuation, angular averaging and mass integration use known mathematics. The opacity profile, storage prescription and gravitational interpretation are proposed. The exact-one-third galaxy fit gives C_(1/2)=4.72858924e7 solar masses/kpc^3, k0=0.2039029004/kpc and s=2.770766589. C is an effective accumulated normalization, not a measured external radiation inventory or a derived universe age. The nonretained channel is assumed to escape; detailed release and support dynamics remain unspecified.

The rotation observable is calculated using established Newtonian bookkeeping,

\[
M_{\rm dep}(<r)=4\pi\int_0^r\rho_{\rm dep}(s)s^2ds,
\qquad v^2(r)=v_b^2(r)+\frac{G M_{\rm dep}(<r)}{r}.
\]

This last relation assumes the deposited component has the stipulated ordinary gravitational response. It is not a derivation of how microscopic companion energy produces a stationary field.

## Galaxy rotation: improvement, but MOND remains better

The matched comparison uses 149 SPARC galaxies and 3150 radial measurements, with 89 training, 29 validation and 31 test galaxies. Distances, inclinations and stellar mass-to-light choices are held common. All partitions have been exposed during development; they are not blind discovery evidence.

| Model | Validation RMS speed error | Test RMS speed error |
|---|---:|---:|
| Ordinary matter | 58.22 km/s | 47.77 km/s |
| Exact-one-third companions | 32.49 km/s | 23.59 km/s |
| Simple MOND, training-fitted a0 | 26.88 km/s | 16.40 km/s |

Companions have three fitted shared constants; the tested MOND prescription has one. Across the 60 comparison galaxies, companions improve absolute speed RMS relative to ordinary matter in 56, and beat MOND in 20. For logarithmic errors, the counts are 60 and 25 respectively. These descriptive counts complement rather than replace the aggregate errors, and do not justify selecting only favorable galaxies. The shared halo scaling and per-galaxy NFW fits have different information budgets and are reported separately in the [matched comparison](../../research_work/results/isotropic-galaxy-transfer/model-comparison-report.md).

## Coma: comparable shapes, no normalization transfer

| Shape | All-six-bin residual sum | Predict-each-omitted-bin residual sum |
|---|---:|---:|
| Transparent companions | 3.68 | 6.86 |
| Strong interception | 3.72 | 10.79 |
| NFW | 3.85 | 6.82 |

The NFW all-bin result approximately reproduces the published 3.87. The earlier inner-three/outer-three ranking did not persist in the omitted-bin comparison. The data are six vector-reconstructed plotted shear values with diagonal errors, not a raw shape likelihood. The point-mass curve is an inadequate diagnostic control, not a candidate explanation of the extended cluster.

Crucially, both companion Coma shapes fit their amplitude and radial scale to Coma. A constant eta is absorbed by the fitted amplitude, so these fits provide no evidence for the one-third exponent, available energy supply, or transfer of galaxy calibration into clusters. Their T=0 and T=100 optical depths are fixed exploratory choices, not values predicted by the galaxy law for Coma. Nor is the compact-baryon MOND shape a full gas-and-star MOND cluster calculation. See [Coma diagnostics](../../research_work/results/isotropic-galaxy-transfer/cluster-comparison-detail-report.md).

## New exact-one-third lens test

The same exact-one-third galaxy capture constants were transferred to six lens galaxies using the existing regular optical branch. Only ordinary stellar mass and orbital anisotropy were fitted to the inner stellar measurements; outer motions and lens angles were predicted. Independent population-mass proxies were required for missing rest-frame 3.6-micron luminosities.

| Population proxy | Outer motion residual-square sum | Lens-angle fractional RMS |
|---|---:|---:|
| Chabrier | 44.3222 | 13.9231% |
| Salpeter | 42.7040 | 14.0937% |

The exact exponent does not repair the mismatch. The previous fitted-exponent version gave 13.9377%/14.1097%, and the earlier intercepted-source model gave 13.1117%. Better outer motions relative to the older model do not amount to simultaneous lensing and motion success. Population and morphology proxies remain limitations, not license for per-lens retuning. [Executed audit](../../research_work/results/isotropic-galaxy-transfer/cross-test-audit-report.md).

## Exact-one-third redistribution follow-up

Four conservative redistribution variants have now been tested, retaining the one-third law and all three capture constants. They modify only the location of the existing deposits through the known identity rho_new(r)=(1-f)rho_0(r)+f*rho_0(r/s)/s^3. The proposed physical interpretation is companion migration; the mathematical operation conserves the deposit inventory but does not derive transport, support or gravitational work.

Shared dilation and retention-conditioned dilation provide negligible improvement. Partial redistribution first preferred its compact-radius boundary; a separately declared wider-range follow-up converged across six starts to an interior solution, f=0.00369849 and s=0.10217474. Thus approximately 0.37% of the inventory is concentrated inward. All redistribution parameters were fitted on training galaxies only, then frozen for comparison galaxies and lenses.

| Quantity | Original exact-third | Compact redistribution |
|---|---:|---:|
| Validation speed RMS | 32.495 km/s | 31.666 km/s |
| Test speed RMS | 23.591 km/s | 22.591 km/s |
| Chabrier lens RMS | 13.9231% | 13.9149% |
| Salpeter lens RMS | 14.0937% | 14.0827% |

This is a modest rotation improvement with two additional shared parameters. The tiny lens change does not resolve the discrepancy or establish significant improvement. MOND remains better on the matched rotation benchmark. The initial partial variant worsens lensing, and every tested variant is retained in the [redistribution report](../../research_work/results/isotropic-galaxy-transfer/redistribution-report.md). This post-boundary exploration does not constitute a new blind test, and the original model remains separately archived. No new cluster normalization-transfer result follows from these tests.

## Other domains must remain in the assessment

The [halo deposition map](../../research_work/results/isotropic-galaxy-transfer/halo-deposition-map-report.md) records the full fitted NFW parameters and compares their finite-aperture masses with the existing one-third deposit inventory. For the standard-geometry targets inside a diagnostic 5 effective radii, four of six have enough current inventory in principle; J1112 and J1630 require more than that fixed inventory. This is not a verified photon-supply test, and 5 effective radii partly extrapolates beyond the stellar observations. The standard target is compared with the original physical deposit profile, not a new standard-geometry companion calculation. Positive density matching would reproduce gravity only with the appropriate stress/support and consistent geometry; no capture-and-migration mechanism has yet been shown to yield those profiles. Untruncated NFW has no finite total mass, so its normalization is not a total halo inventory.

The [matched NFW/geometry diagnostic](../../research_work/results/isotropic-galaxy-transfer/nfw-geometry-report.md) now fits the same six systems with stars plus NFW and exact catalog lens angles. Standard flat-FLRW benchmark distances improve all six NFW stellar scores relative to our regular optical geometry. For J0037/J1204/J1402, scores fall from 88.59/54.38/195.85 to 20.20/4.05/43.07. The latter still has substantial residuals, and several solutions have orbital or halo boundaries. This implicates the current optical branch alongside mass-profile assumptions; it neither establishes complete NFW success nor rejects every nonexpanding geometry. FLRW is used only as a comparison model. No companion standard-geometry transfer has yet been computed, so this is not a complete two-by-two test. Halo amounts are target-fitted and cannot be equated with a fixed companion energy budget.

The subsequent [lens/profile inverse diagnostic](../../research_work/results/isotropic-galaxy-transfer/lens-profile-compatibility-report.md) keeps this exact-one-third inventory but permits target-specific positive mixtures of seven conservative dilations. Forcing the catalog lens angle then leaves severe stellar-motion errors for J0037, J1204 and J1402 in both population cases. J1621 admits a low-residual solution at the orbital bound; J1112 and J1630 have much smaller lens-induced tension. This is an inverse fit, not a prediction. The limited profile basis, fixed geometry and constant orbital anisotropy do not justify ruling out every redistribution. The result motivates testing optical geometry and stellar/orbital modeling before adding an arbitrary light-bending response.

| Domain | Present assessment | Why it matters |
|---|---|---|
| Redshift | Fractional-loss law z=exp(integral alpha ds)-1 is calculable; alpha is empirical | The exponential is known mathematics; fitting redshift does not derive its physical cause |
| Supernova timing and brightness | A postulated event stretch and regular optical response improve conditional brightness predictions | Stretching each photon's wavelength does not alone derive the spacing between arriving events |
| Optical source energy | A restricted homogeneous temporal completion has a severe supply mismatch | A favorable optical fit does not guarantee its required radiation source exists |
| Solar-system source/retention tests | No tested common photon-feeding scaling reproduces the full Sun/Earth/Moon gravity across bodies | Full ordinary gravity and extra galactic gravity are different targets |
| Milky Way three-dimensional motions | Not completed with independent matched observations and the current exact-third model | Radial rotation alone does not establish the vertical field or bulge behavior |
| Cluster collisions and spatial lens maps | No completed prediction | A radial profile cannot explain offsets between gas, galaxies and lensing structure by itself |
| Microwave spectrum and angular correlations | Candidate mechanisms recorded; no joint predictive fit | A persistent background proposal must reproduce more than its existence |

The [regular optical report](../../research_work/results/brightness-distance-consistency/regular-area-report.md) and [restricted source-growth failure](../../research_work/results/brightness-distance-consistency/optical-growth-report.md) distinguish empirical fits from mechanism requirements. No fixed universe age or size is inserted to force success. All six project goals remain open.

## Near and distant companion streams

The [source-flow diagnostic](../../research_work/results/isotropic-galaxy-transfer/companion-streams-report.md) applies the proposed conversion to known inverse-square transport, F_c=L[1-exp(-alpha D)]/(4*pi*D^2), while retaining the one-third factor. In the small-conversion limit a single source gives F_c proportional to 1/D. Near/far illumination contrast is approximately (D+R)/(D-R); sources at 100 receiver radii produce only about 2% contrast. A uniform population in a shell contributes j[1-exp(-alpha D)]dD, allowing distant populations collectively to dominate without individually brighter sources. This is conditional geometry, not a measured cosmic energy supply.

For the cached M87 source field, 797 external galaxy entries give about 53.1% of catalog center intensity from 1-5 Mpc, with angular second-moment eigenvalues 0.127, 0.359 and 0.514. The catalog extends only to about 72 Mpc and cannot establish the distant population's contribution. Integrated galaxy luminosities are estimated inputs, not resolved stars or source histories. No cluster retention proxy is invented; a common one-third factor cancels from fractional shares.

Streams crossing do not automatically merge or bind. Lightlike bending scales are of order 10^-5 radians for the existing lenses; full ray tracing and companion-specific scattering remain unexecuted. The next mechanism test must map source paths through capture into a radial and angular deposit distribution before comparing it to halo targets. This calculation supplies that test's inputs and does not improve the reported rotation or lens scores.

## Straight-stream capture follow-up

The [executed capture test](../../research_work/results/isotropic-galaxy-transfer/stream-capture-report.md) holds the one-third law, k0 and original capture scales fixed and integrates external companion rays inside each diagnostic 5 Re region. Five finite source distances and a distant background are compared with all existing geometry/population halo targets. Boundary illumination is conditional; internal sources, ongoing photon conversion inside the boundary, exterior opacity, steering and migration are excluded. Constant history is required to equate deposited-power shape with accumulated density. Normalizing both enclosed profiles at the boundary isolates placement and does not verify energy supply.

With the standard comparison targets, separately fitted nonnegative source-distance mixtures reduce maximum cumulative fraction errors to 2.44, 5.52 and 5.83 percentage points for J1112, J1621 and J1630. They require predominantly or entirely the nearest tested shell, D/R=1.1. J0037, J1204 and J1402 remain at 17.97, 49.94 and 65.11 points, with inadequate central deposition. Under the retained geometry none of the six mixtures improves on distant illumination within this basis. The errors are normalized radial differences, not lens residuals or observational significance. Both population proxies and all cases are archived.

This test distinguishes a source-location adjustment from a new capture or migration mechanism. Nearby illumination can move deposited weight outward, but does not supply the missing central concentration in the difficult targets. Boundary-fitted source weights are not observed stellar populations and cannot be transferred as predictions. The original halo-fit limitations and incomplete source histories remain; no new rotation/lens fit or complete halo-formation mechanism is claimed.

## Post-capture inward migration

Four [migration variants](../../research_work/results/isotropic-galaxy-transfer/migration-variants-report.md) were tested against the finite halo targets: uniform contraction, a mobile fraction plus a stationary fraction, relaxation toward a support radius, and a retention-conditioned relaxation rate. These use known relaxation and cumulative-distribution transformations with a proposed companion interpretation. Capture constants and the exact one-third law remain fixed. Both distant illumination and a hypothetical nearby shell were tested under both geometries and population proxies. Parameters are shared within each branch; fits to five targets were frozen for evaluation on the sixth. The targets have already been inspected, so this is conditional cross-validation rather than blind observational validation.

For the retained geometry and distant input, partial migration has the strongest improvement: cumulative-fraction RMS is 26.74 percentage points before migration, 16.80 in the shared fit and 20.01 in omitted-target transfer. The shared fit moves 35.8% of deposits to about 15.0% of their previous radii. It improves J0037, J1204 and J1402 but worsens J1112, J1621 and J1630. This fraction is fitted to normalized halo profiles and does not replace the separate 0.37% rotation-trained result. Uniform and retention-conditioned variants improve less; the supported variant selects zero support radius and reduces to uniform contraction.

With standard-geometry targets and distant input, no tested variant improves aggregate omitted-target RMS; partial migration gives 24.32 versus the 23.55 no-migration baseline. A nearby-input partial variant improves its own worse baseline but remains worse than the distant baseline. Accordingly, inward mobility can address central deficits but has not yielded a shared solution across halo types. The transport maps conserve inventory, not the complete gravitational/kinetic energy budget, and do not derive a stopping stress or evolve self-gravity. No new velocity or lensing likelihood fit, absolute supply verification or completed research goal follows.

## Saturable bound-state capacity

The [capacity-limited settling test](../../research_work/results/isotropic-galaxy-transfer/saturation-report.md) adds three proposed receiving capacities: constant density, density proportional to a stellar Hernquist proxy, and density proportional to baryonic acceleration. The proxy uses original companion-fit stellar masses rather than target halo mass. These are phenomenological extensions built from known density and energy bookkeeping; no microscopic capacity law is derived. Existing capture and one-third retention remain fixed.

Packets fill available states between their proposed migration radius and original radius. Overflow remains in a tracked traveling reservoir rather than being deleted or normalized into the bound halo. Total supplied inventory is the original deposited mass inside diagnostic 5 Re, not the target halo mass. The calculation compares bound cumulative mass divided by supplied inventory with normalized targets, so it does not validate absolute mass or include the gravitational field of traveling overflow. Processing order is a declared assembly assumption.

With retained geometry and distant sources, omitted-target cumulative RMS is 24.64, 28.08 and 25.46 percentage points for constant, stellar and acceleration-linked capacity respectively, versus 26.69 without migration and 19.98 for archived partial migration evaluated on the same radial grid. The shared distant fits have no active capacity effect and reduce to unsaturated contraction. A near-source acceleration-linked case gives 23.62 points, better than its 31.66 baseline but worse than the prior near-source partial prediction of 20.47. Neither constitutes a common halo solution.

Packet accounting closes to numerical precision and bound cells respect capacity. Doubling the grid changes selected cumulative profiles by less than 0.8 percentage points, while reversing packet arrival order changes one by 43.8 points. Some stellar-capacity transfers leave roughly 28% of a target inventory unbound. These are material formation-history and traveling-field limitations, not evidence of disappearing energy. Capacity alone has not supplied a unique stable halo, a support stress, or a complete gravitational-work budget. All goals remain open.

## Wave and orbital support candidates

The [wave/particle comparison](../../research_work/results/isotropic-galaxy-transfer/wave-particle-report.md) replaces an arbitrary stopping radius with two restricted equilibrium constructions. A Gaussian massive scalar trial state minimizes gradient, stellar-potential and self-gravitational energy, with and without repulsive P=K rho^2 support. Randomly oriented circular particle orbits instead balance the stellar-plus-deposit attraction under a shared angular-momentum prescription. These are known physical/mathematical structures applied to a proposed photon-companion origin; no massless graviton-halo mechanism is established. Local scalar support concepts are discussed by [Hui et al.](https://arxiv.org/abs/1610.08297), whose cosmological formation assumptions are not adopted here.

With retained geometry, omitted-target profile RMS is 39.08 points for the simple wave, 33.25 for repulsive waves, and 25.78 for distant-input circular particles, versus 26.74 without migration and 20.01 for earlier partial migration. The nearby-input particle result is 26.85 versus the earlier partial 20.51. Standard-geometry transfer does not repair this result. Repulsive shared fits favor the smallest allowed gradient coefficient, so they do not identify a preferred wave scale. None beats the previous profile benchmark.

Selected wave states have positive curvature against Gaussian radius perturbations, and selected particle orbits have positive epicyclic frequency and numerical force balance. These are restricted support checks, not arbitrary-mode or collective stability proofs. Inconsistent ordered particle shell configurations were rejected before scoring; all selected transfers were valid. Wave radial quadrature refinement changes radii by at most 0.022%, and particle shell refinement changes cumulative fractions by under 0.05 percentage points. The supplied inventory remains fixed, including material outside the comparison region; its absolute adequacy and the formation energy/angular-momentum ledger remain unresolved. No new lensing or stellar-motion fit is claimed.

## Local compact/extended exchange

Three [local exchange rules](../../research_work/results/isotropic-galaxy-transfer/reservoir-exchange-report.md) use the known two-state stationary fraction f=k_on/(k_on+k_off), with proposed rate ratios set by stellar binding energy, binding-energy release over inward transport, or binding with a local density penalty. Original stellar/capture inputs set the rates; no target halo discrepancy drives them. The selected fraction moves to s*r while the remainder stays extended. Capture and one-third retention remain fixed. Both the rate scale and s are phenomenological; all rates use the initial environment, so this is not a self-consistent time-dependent equilibrium.

For retained geometry and distant illumination, omitted-target cumulative RMS is 20.34, 20.15 and 21.85 percentage points respectively, versus 20.01 for the earlier partial-migration model. The released-binding shared fit uses s=0.1391 and v_ex=720.84 km/s, yielding compact fractions of 30.9-37.6% from local inputs rather than individual fitted fractions. This gives a possible rationale for the earlier roughly 36% split at similar accuracy, but does not improve the primary benchmark. Several extended halos still worsen relative to no migration.

Under standard comparison geometry the occupancy penalty gives 23.47 points with distant input, compared with 23.55 without migration and 24.32 for earlier partial migration. This small difference is not a significance claim or new observational success. Inventory, positivity, inward ordering and resolution checks pass. Work released during settling, final-state support, updated occupations, self-gravity and the absolute photon supply remain unresolved. No new lensing or velocity fit or completed research goal follows.

## Feedback from the deposited gravitational field

The [coupled feedback test](../../research_work/results/isotropic-galaxy-transfer/feedback-report.md) updates the exchange-driving potential using stars plus both compact and extended deposits; the occupancy penalty uses current density. Capture, original finite inventory and baryonic gravity remain fixed. Iteration solves a coupled fixed point, with compact radii always referenced to the original source locations. It does not simulate repeated contraction, physical time, evolving baryonic orbits or photon arrival histories.

With old parameters fixed, retained-geometry omitted-target errors worsen to 21.69, 22.66 and 25.63 percentage points for binding, released binding and occupancy rules. Refitting shared parameters gives 20.04, 20.59 and 20.87, respectively. Binding is essentially tied with the earlier 20.01 partial-migration benchmark; the difference is below the refinement scale. Standard-geometry occupancy gives 22.93 versus 23.47 previously, without establishing a new joint observational fit. The updated shared binding rule predicts compact fractions of 33.3-42.5%.

All executed fixed-point and refit-grid cases converge; selected opposite-seed solutions agree within 2.5e-8 in cumulative fraction. This is consistency of the numerical closure, not global uniqueness or dynamical stability. Selected radial refinement shifts cumulative fractions by up to 0.162 percentage points. Turning feedback off reproduces the earlier prescribed-field profiles within 0.118 points. Binding-work diagnostics expose energy that requires a support/heat/radiation channel; they are not an added photon supply.

Frozen traveling/stored inventory ratios of 0, 0.01 and 0.1 are recorded as prescribed Newtonian energy-density sensitivities, with additive enclosed traveling mass separately retained. Their amount and source history are not measured. Directional radiation stresses, a relativistic field equation, evolving ordinary matter and the complete formation energy ledger remain absent. No all-energy gravitational model or new lens/motion likelihood success is claimed.

## Current Milky Way circular-speed transfer

A [new frozen Milky Way transfer](../../research_work/results/isotropic-galaxy-transfer/milky-way-current-report.md) uses the exact-third capture constants and recent migration/feedback parameters without fitting Milky Way speeds. It compares 38 cached [Eilers et al.](https://arxiv.org/abs/1810.09466) circular-speed summaries across approximately 5-25 kpc, preserving two previously declared ordinary-star/gas baselines. These are Jeans-inferred circular speeds rather than individual star velocities; the paper's separate halo inference is not imported into the prediction.

Exact-third RMS is 6.76/10.41 km/s for ordinary baselines I/II, compared with 52.57/62.34 for ordinary matter alone. The small SPARC rotation-trained partial redistribution gives 6.70/10.02, while the separately halo-trained large partial fraction gives 35.57/28.56. Latest binding feedback gives 15.08/8.32 and released-binding feedback 12.03/5.75. Thus the new feedback improves one ordinary baseline but worsens the other; no variant dominates. At 8.19 kpc under baseline I, observed circular speed is 228.86 km/s versus exact-third 226.60 and ordinary-only 188.43.

The calculation uses Rd=2.6 kpc as a star-count proxy and estimates 3.6-micron luminosity from archived stellar mass with nominal M/L=0.5. These are conditional inputs, not new independent infrared measurements. All predeclared Rd and luminosity sensitivity cases are preserved without selecting a favorable fit. Feedback uses the spherical angular average of the archived stellar potential and a thin-gas monopole; the final baryonic disk force remains the archived one. This is a spherical migration approximation rather than full barred-Galaxy dynamics. Resolution changes fiducial speeds by under 0.001 km/s, but that numerical precision is not an observational uncertainty claim.

The radial result is useful evidence of conditional predictive performance, not verification of photon supply, capture microphysics, individual orbits or the vertical field. These observations were seen previously, so frozen transfer is not blind discovery. No dark-halo source, expansion assumption, Milky Way velocity refit or completed broader goal is introduced.

## Comparison claims that the evidence permits

A common law that predicts both galaxy and cluster observables with fewer unsupported additions would be a meaningful advance. The present results do not yet show that: MOND is stronger in the matched galaxy test, while the cluster exercise is freely fitted and the lens transfer remains mismatched. A cluster shape fit cannot cancel a galaxy prediction deficit or establish a missing energy budget.

The identity of dark matter remains unresolved, with particle searches placing constraints rather than an established identification. That is distinct from absence of gravitational evidence. Our traveling companions and retained states have also not been detected. Moreover, a nonluminous deposited component gravitating through ordinary stress-energy would itself function as an effective dark component; the proposed distinction is its origin and transport. The fair comparison is predictive accuracy, common parameters, energy accounting and independent discriminating observations, not which unobserved entity has a preferred name. See the [PDG dark-matter review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-dark-matter.pdf).

The next empirical target is a gas-and-star-calibrated cluster calculation followed by a different cluster with shared parameters frozen. LoCuSS and CLASH provide relevant published observations, but no new ensemble profile fit was executed for this supplement. Existing reports of good NFW fits must be included in the comparison, not discarded. See [the next-sample specification](../../research_work/results/isotropic-galaxy-transfer/cluster-comparison-next-sample.md).
