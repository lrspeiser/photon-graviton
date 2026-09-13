# Cross-scale performance of radiation-conditioned companion deposition

Working research assessment, 13 September 2026. This is an evolving supplement to the historical manuscript, not a peer-reviewed result or a rebuilt PDF. Numerical claims below refer to the linked executable records. Results are conditional on a hypothetical nonexpanding universe and the explicitly stated geometry and matter inputs.

## Abstract

We investigate a phenomenological proposal in which photons transfer energy to traveling companions and some companion energy is retained in gravitational reservoirs. A bounded retention prescription with exponent one-third improves galaxy rotation predictions relative to an ordinary-matter-only baseline. On 60 comparison galaxies, it improves absolute speed errors for 56 objects but outperforms the tested simple MOND prescription for only 20. Its aggregate rotation error remains larger than MOND's. Two freely normalized companion density shapes describe the six reconstructed Coma shear measurements comparably to NFW; this does not test the galaxy-calibrated retention normalization. A new transfer of the exact-one-third galaxy parameters into six lens galaxies leaves approximately 14% lens-angle RMS discrepancies. These results support continued investigation of an empirical extra-gravity prescription, but do not establish a unified photon-to-gravity mechanism or superiority to MOND or dark matter.

## Model definition and provenance

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

## Comparison claims that the evidence permits

A common law that predicts both galaxy and cluster observables with fewer unsupported additions would be a meaningful advance. The present results do not yet show that: MOND is stronger in the matched galaxy test, while the cluster exercise is freely fitted and the lens transfer remains mismatched. A cluster shape fit cannot cancel a galaxy prediction deficit or establish a missing energy budget.

The identity of dark matter remains unresolved, with particle searches placing constraints rather than an established identification. That is distinct from absence of gravitational evidence. Our traveling companions and retained states have also not been detected. Moreover, a nonluminous deposited component gravitating through ordinary stress-energy would itself function as an effective dark component; the proposed distinction is its origin and transport. The fair comparison is predictive accuracy, common parameters, energy accounting and independent discriminating observations, not which unobserved entity has a preferred name. See the [PDG dark-matter review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-dark-matter.pdf).

The next empirical target is a gas-and-star-calibrated cluster calculation followed by a different cluster with shared parameters frozen. LoCuSS and CLASH provide relevant published observations, but no new ensemble profile fit was executed for this supplement. Existing reports of good NFW fits must be included in the comparison, not discarded. See [the next-sample specification](../../research_work/results/isotropic-galaxy-transfer/cluster-comparison-next-sample.md).
