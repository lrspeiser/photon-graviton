# Causal companion tests and closed-universe circulation

Research checkpoint — 9 September 2026

## Scope and conclusion

This is a hypothetical nonexpanding universe using real astronomical observations as comparison data. The calculations do not establish a photon–companion interaction, replace dark matter in our universe, or demonstrate a stable static cosmology.

The strongest new numerical lead is a baryonic-depth dependence in an empirical extra-acceleration law. A version defined consistently at each radius improves the galaxy validation score and approximately matches the eleven cluster mass estimates without a cluster normalization fit. Its cluster uncertainty is broad, and it does not yet connect gravity to independently supplied companion energy.

Measured bolometric luminosities make the photon energy-budget problem concrete. Catalogued external illumination, ordinary energy storage, and a simple deeper-well retention law do not close the causal connection. Closed geometry can recycle trajectories, but cannot multiply the energy carried by a photon. With the currently fitted propagation history extrapolated indefinitely, light has only a finite remaining path length even in a static space.

## Data and method

- SPARC: 175 published disk-galaxy mass models; 149 pass the existing quality/inclination/radius cuts, providing 3,150 velocity points. Disk and bulge mass-to-light ratios remain fixed at 0.5 and 0.7. The existing 89/29/31 training/validation/test assignment is reused. These are not blind new galaxies.
- DustPedia / Nersesian et al. (2019): newly retrieved CDS catalog J/A+A/624/A80. There are 814 usable SED fits among 875 catalog rows. Exact canonical-name matching gives 26 usable SPARC overlaps: 19 training, 3 validation, 4 test. Their current total bolometric luminosities, positions and distances are new inputs to these tests. Full radiation histories have not been reconstructed.
- The THEMIS SED fit supplies the primary luminosity; DL14 supplies a model-choice comparison. Target luminosities are rescaled by (SPARC distance / DustPedia distance)^2 to match the rotation-curve distance. Neighbor luminosities and positions retain the catalog distance system. This is a limitation of the environmental calculation, not a self-consistent distance reanalysis.
- X-COP: the eleven previously transcribed clusters with overlapping hydrostatic mass and gas-fraction inputs. Stellar fraction is assumed to be 0.02, with 0.01–0.03 sensitivity. Published hydrostatic/geometry assumptions remain. These are not independent lensing measurements.
- Fits minimize the mean across galaxies of each galaxy's mean squared log10 velocity residual. Reported km/s RMSE is also galaxy-weighted. It is not a reduced chi-square. Distance, inclination, stellar mass-to-light, and radial covariance are not marginalized.
- For the small DustPedia overlap, five deterministic name-hash folds provide out-of-fold predictions for every object. Every candidate is reported; no claim of pristine model-selection holdout is made. The old 4-object test scores are in results.json but too small to drive the conclusions.

## 1. Independently estimated energy supply

At the last measured radius R, the required Newtonian-equivalent extra mass is

M_X = R [v_obs^2 - v_bar^2] / G.

For the ordinary, cold deposited-energy benchmark, E_X = M_X c^2. This is deliberately an ordinary mass-energy test; a relativistic gravitational source needs its full stress tensor. An unspecified stronger gravitational response cannot be counted as extra supplied energy.

Use an explicit illumination duration T = 10 billion years and constant current bolometric luminosity L. Perfect capture and retention are optimistic assumptions. This duration is a scenario choice, not an imposed age of the fictional universe.

Full-conversion energy: E_full = L T.

Local conversion: E_local = L T [1 - exp(-k R)], where k = 7.7315e-5 per million light-years and R has matching units. This is a generous first-transfer budget; it neglects secondary losses before deposition.

For an external source at geometric source-to-target distance d, under the nearby exponential approximation S = exp(k d), shared photon/companion secondary redshift gives

F_companion = L log(S) / [4 pi d^2 S^2].

This follows from per-emitted-photon E_companion/E_emitted = log(S)/S and the arrival-rate factor 1/S. It assumes no intervening capture, no absorption and isotropic emission. It does not use the redshift from that galaxy to Earth as the target-to-source redshift. The constant-k approximation is appropriate only locally; the revised field history is not exactly an exponential of distance everywhere.

Sum the scalar incident flux over the 813 other usable catalog sources. Distances use Euclidean positions from catalog RA, declination and distance, with a 30-kpc minimum separation to regularize unresolved close pairs. Then E_external = pi R^2 F_sum T. The geometrical cross-section expression is also correct for a summed isotropic bath; one must not multiply the already intercepted power by a second area factor.

| Energy scenario | Median required / supplied energy, 26 galaxies | 10–90% range |
|---|---:|---:|
| Every photon converted completely | 5,428 | 1,876–11,453 |
| Local conversion at the fixed redshift coefficient | 8.89e8 | 4.21e8–3.32e9 |
| Captured companions from catalogued neighbors | 3.74e9 | 1.85e9–6.91e9 |
| Local conversion plus catalogued neighbors | 6.85e8 | 3.94e8–1.41e9 |

These ratios are calculated at each galaxy's last measured radius; that radius is not asserted to be a physical halo edge. The median SED-reported fractional luminosity uncertainty is 5.0%; median DL14/THEMIS luminosity ratio is 1.019. Neither addresses a deficit of this scale. Median external companion energy density is 3.17e-20 J/m^3. The catalog is incomplete and nearby-selected: this is a catalog-limited prediction, not an upper limit on all conceivable external radiation.

Illustrative independent 20% log-distance perturbations change a typical object's computed illumination by roughly 0.90–1.08 at the 95% interval endpoints (median endpoints across targets). This deliberately simple sensitivity is not a published-error likelihood and does not rescale neighbor luminosities with perturbed distances. The no-stretch approximation changes median illumination by 0.83%. The dominant uncertainty remains missing sources, long-term histories and unproved capture physics.

## 2. Can transport or retention predict the radial gravity?

For a restricted positive-density deposition model, let halo radius R_h = 10 R_d, where R_d is the observed stellar disk scale. R_h is prescribed independently of the measured velocity excess. Truncate all deposited mass there.

The local source uses a softened inverse-square production rate,

q_local(r) = k L / [4 pi (r^2 + R_d^2)],

with consistently converted units. This is an approximate spherical source prescription, not a derived stellar radiative-transfer solution. With no escape and duration T, its enclosed ordinary mass is

M_local(<r) = (k L T / c^2) [r' - R_d atan(r'/R_d)],   r' = min(r, R_h).

This explains why a cored density can arise from inverse-square illumination plus local conversion, but its physical normalization is far too small.

For external energy assume uniform deposited density inside R_h, total intercepted power pi R_h^2 c u_c, and perfect retention. Its enclosed mass fraction is (r'/R_h)^3. Predicted v^2 = v_bar^2 + G M_extra(<r)/r. Free universal gains are allowed to test profile/scaling information, but are not evidence that the energy exists.

A retention candidate uses f_ret = x^m/(1+x^m), x = v_depth/(200 km/s), m >= 0. The global depth proxy is v_depth = sqrt[2 G (0.5 L_3.6 + 1.33 M_HI)/R_d]. It is based on baryonic inputs and is not a measured escape speed. It omits a separate integrated bulge correction and molecular gas. Its normalization is degenerate with the gain and storage time.

| Candidate | Free universal parameters | Out-of-fold log RMSE, dex | Out-of-fold RMSE, km/s |
|---|---:|---:|---:|
| Local production with gain | 1 | 0.12890 | 36.21 |
| External catalog illumination with gain | 1 | 0.09716 | 32.50 |
| Local + external, separate gains | 2 | 0.11124 | 28.16 |
| Local + depth-dependent retention | 2 | 0.12890 | 36.21 |
| External + depth-dependent retention | 2 | 0.09731 | 32.64 |
| Mixed + depth-dependent retention | 3 | 0.11126 | 28.16 |
| Constant external bath control | 1 | 0.09585 | 33.13 |
| Local + constant bath control | 2 | 0.08545 | 28.37 |

The different rankings under log and linear residuals matter: one must not choose the more favorable metric after seeing the results. The predeclared fit metric is log residual. Actual catalog illumination is not better than a constant-bath control under that metric. Therefore these fits do not establish the environmental companion connection. A constant bath is a diagnostic control, not a measured source.

Local and mixed retention fits put m at zero on the 19-object training set; they do not find the hypothesized preferential retention in deeper wells. The external-only fit finds m ~ 0.123 but does not improve its out-of-fold score. Fitted local-only gain is 1.62e9; external-only gain is 7.06e9. Mixed fitted gains are 6.91e8 and 4.55e9. These large independent gains explicitly fail the desired shared energy normalization.

For perspective, baryons alone give 56.01 km/s on the matched sample; ordinary local plus external energy changes that only in negligible decimal places. The earlier empirical luminosity-size core gives 22.10 km/s on those objects, but many helped train it. That is a descriptive reference, not an equivalent out-of-fold comparison.

A paired bootstrap of out-of-fold per-galaxy residuals gives mixed-minus-local log-MSE difference -0.00424, with percentile interval [-0.00823, -0.00102]. This conditional comparison ignores overlapping-fold training uncertainty and multiple candidate development. It supports adding profile flexibility; the constant-bath control prevents interpreting it as evidence for catalogued companion flux.

### Boundary diffusion check

For companions entering at the edge and diffusing inward while leaking, stationary source-free interior transport satisfies

D nabla^2 u - u/tau = 0,  ell = sqrt(D tau).

The regular spherical solution is u(r) proportional to sinh(r/ell)/r. Its enclosed mass fraction is

[y cosh(y) - sinh(y)] / [Y cosh(Y) - sinh(Y)], y=r/ell, Y=R_h/ell.

This becomes uniform as ell/R_h grows and concentrates near the edge when ell is small. Normalize to the same total external supplied-energy template and allow one universal gain. This is a shape test; the gain is not a solved stationary capture normalization.

| Diffusion length / halo radius | Out-of-fold log RMSE | RMSE km/s |
|---|---:|---:|
| 0.1 | 0.17210 | 48.02 |
| 0.3 | 0.11169 | 37.59 |
| 1.0 | 0.09826 | 32.96 |
| Uniform-density limit | 0.09716 | 32.50 |

The narrow edge deposit performs poorly. Redistribution through the interior is preferable within these assumptions. Changing the local model's outer boundary from 10 R_d to 5 or 20 R_d leaves its poor four-object test RMSE at 75.31 or 66.24 km/s (69.02 baseline). This is a sensitivity check, not a calibrated halo boundary.

## 3. A more promising gravitational-response formula

Keep the force-law route explicitly separate from companion energy accounting. Use the same formula at every measured radius:

 g_X = A a_0 (g_b/a_0)^q [sqrt(2 r g_b)/(200 km/s)]^m,
 g_total = g_b + g_X,
 a_0 = 1e-10 m/s^2.

Fit A, q and m on the original 89 SPARC training galaxies only. The pointwise depth factor sqrt(2 r g_b) uses the same definition for galaxies and clusters. It is a baryonic circular-speed proxy, not an exact escape speed or a local covariant field scalar; r and the center must ultimately arise from a field solution.

Best fit:

A = 0.699, q = 0.3131, m = 0.2651.

| Formula | Validation log RMSE | Test log RMSE | Test RMSE km/s | Cluster median required/predicted extra mass |
|---|---:|---:|---:|---:|
| Earlier acceleration power law, no depth term | 0.09743 | 0.07914 | 17.20 | 1.899 |
| Global compactness depth factor | 0.09420 | 0.07487 | 16.56 | 1.140 |
| Consistent pointwise depth factor | 0.09350 | 0.07376 | 16.95 | 0.920 |

The pointwise law's cluster ratio ranges from 0.815 to 1.145. A ratio of 0.920 means predicted extra mass is about 8.6% above required at the median ratio, not that all total cluster masses are accurate to 8%. No cluster-specific gain was fitted. The global compactness variant uses R_d for galaxies and R_500 for clusters, making its transfer less interpretable; the pointwise version removes that particular mismatch.

In 120 galaxy-bootstrap refits, the pointwise depth exponent's percentile interval is approximately 0.032–0.473; the predicted median cluster required/predicted ratio spans 0.605–1.577. These are modest-size exploratory bootstrap intervals with fixed observational nuisance parameters, not complete credible intervals. The assumed cluster stellar fraction 0.01/0.02/0.03 changes the median ratio to 0.959/0.920/0.885.

This is the strongest lead for further development. However, it fits a gravity response directly, has an additional parameter, and does not calculate that response from photons. Cluster gas mass also uses published hydrostatic-related quantities and distances. A new lensing prediction needs both metric potentials; rotation constrains one combination and does not determine the other. No lensing likelihood has been run here.

## 4. Could closed geometry and repeated circulation provide enough energy?

Potentially it increases exposure or lets a continuously replenished reservoir accumulate. It does not, by itself, make the present budget sufficient.

A positively curved closed three-sphere can have returning geodesics; compact flat topologies can also have returns. Curvature alone does not guarantee repeated intersections with galaxies. A three-sphere is an illustrative geometry, not a demonstrated background solution of the candidate theory.

For constant loss per path length, let one loop retain fraction s of a photon's energy. After N loops E_gamma = E_0 s^N. Total energy ever transferred out of that photon is E_0(1-s^N), never more than E_0. Repeated encounters can make eventual capture more probable but cannot be counted as repeated creation of the original energy. If the photon is captured, its trajectory ends unless an explicitly powered reemission process supplies a new photon.

The proposed shared secondary redshift makes the limit tighter for *free companion* supply. With no capture and n/n_emit = S,

E_gamma/E_0 = 1/S,
E_companion/E_0 = log(S)/S,
E_timing/E_0 = 1 - [1+log(S)]/S.

The free companion energy peaks at 1/e of the original energy when S=e; for indefinitely increasing S, it too decays into the timing sector. Using that timing-sector energy as a halo source requires an additional transfer and gravitational-response law.

### A consequence of the currently fitted field history

For stationary spatial geometry, v_light=c/n and dn/dt=gamma n^p, with n_today=1, gamma=7.7315e-11 yr^-1 and p=0.263906,

R_future(n) = [1-n^(-p)]/(k p).

Therefore R_future,max = 1/(k p) = 49,010 million light-years, about 49 billion light-years. This is a conditional future path budget in the fictional model, not an observed cosmological horizon. It follows from extrapolating the fitted law indefinitely, which observations have not justified.

| Closed-path circumference | Maximum future path / circumference | Photon energy left after first return |
|---|---:|---:|
| 1 billion light-years | 49.01 | 92.49% |
| 10 billion light-years | 4.90 | 42.12% |
| 30 billion light-years | 1.63 | 2.76% |
| 100 billion light-years | 0.49 | No completed return |

Thus the current positive-p history does not permit infinitely many future circuits along a fixed nonzero circumference. Returning to p=0, or changing the remote future history, can remove this path bound, but not energy conservation. These circumferences have not been fitted to brightness, repeated images, or other geometry observations.

### An external-background scale

Driver et al. (2016) estimate integrated optical plus infrared background intensity of approximately 24+26 = 50 nW m^-2 sr^-1. This is ordinary background light, not a companion detection. Its equivalent energy density is u_gamma=4 pi I/c = 2.10e-15 J/m^3.

The 26-galaxy median external companion density needed for perfect interception and retention over 10 billion years is 1.01e-10 J/m^3, roughly 48,000 times that benchmark. If a companion bath had only the benchmark energy density, the required accumulation time at fixed interception would be about 4.8e14 years.

In a constant-gamma, spatially mixed stationary illustration,

 du_c/dt = gamma u_gamma - (gamma + lambda_capture + lambda_escape)u_c.

With shared secondary redshift and nonnegative capture/escape, u_c <= u_gamma in steady state. Returning trajectories do not change that balance. If companions instead do not redshift and no capture depletes them, u_c ~ gamma u_gamma t; reaching the required density from this benchmark would take about 6.2e14 years. That is an alternative law, not the shared-redshift model; it also requires sustained luminosity, a reservoir that gravitates before capture, and stable background dynamics. These are optimistic order-of-magnitude illustrations, not upper bounds on every fictional reservoir.

## Pros, cons and next decision

| Candidate | What worked | What remains problematic | Status |
|---|---|---|---|
| Ordinary local photon conversion | Produces a plausible cored shape under simplified illumination | Energy shortfall near a billion at fixed local loss rate | Quantitatively insufficient under tested assumptions |
| Catalogued external illumination | Adds a useful broad spatial component | Billion-scale gain; constant-bath control performs as well or better | No demonstrated causal environmental signal |
| Preferential retention in deep baryonic wells | Clear, bounded shared law can be tested | Local/mixed best fit removes the dependence; normalization degeneracy | Not supported by this restricted test |
| Boundary capture plus diffusion | Broad penetration performs better than a thin edge deposit | Needs capture and redistribution dynamics and energy supply | Useful shape constraint |
| Depth-dependent gravitational response | Better galaxy validation and much closer cluster transfer | Empirical force law; uncertain transfer; no photon normalization or lensing action | Most promising mathematical lead |
| Closed, circulating universe | Can increase encounter probability and retain globally supplied energy | No energy multiplication; current p>0 gives finite future optical path; no static solution | Worth a separate constrained branch, not a budget solution yet |

The focused next milestone is to derive a field response that reduces to the tested pointwise force law while making its coefficient depend on independently calculated companion supply. It must predict both material accelerations and lensing, and specify what gravitates before capture. Do not yet replace the paper's energy budget with the successful empirical response or claim that the circulation option supplies the deficit.

The emitter–propagation–receiver interaction, exact clock-ratio protection, shared photon/tensor dynamics, nonexpanding background stability, and an independent lensing test remain open. This work has tested restrictive submodels and improved the empirical targets; it has not completed the underlying theory.

## Reproduction and sources

Run from repository root with Python 3 and NumPy/SciPy installed:

```sh
python companion_causal_test/run.py
python companion_causal_test/followup.py
python companion_causal_test/circulation.py
```

The follow-up script includes 240 bootstrap optimizations and takes longer. Data are local; scripts do not need network access. Protocol, source hashes, energy budgets, every out-of-fold prediction, fitted parameters, cluster comparisons and sensitivity results are included. The protocol predates the first candidate fits in this checkpoint; follow-up controls were designed after inspecting those fits and are exploratory.

- SPARC: Lelli, McGaugh & Schombert (2016), https://arxiv.org/abs/1606.09251 ; https://astroweb.case.edu/SPARC/ . Original rotation-curve references are in SPARC Table 1.
- Nersesian et al. (2019), A&A 624 A80, https://arxiv.org/abs/1903.05933 ; CDS https://cdsarc.cds.unistra.fr/ftp/J/A+A/624/A80/ . THEMIS and DL14 tables and their ReadMe are included.
- Bianchi et al. (2018), https://arxiv.org/abs/1810.01208 . Bolometric starlight and dust-reprocessed light must not be double counted.
- Ettori et al. (2019), https://arxiv.org/abs/1805.00035 ; Eckert et al. (2019), https://arxiv.org/abs/1805.00034 . Cluster table transcription and prior comparison are in companion_wave_test.
- Driver et al. (2016), https://arxiv.org/abs/1605.01523 . The 50 nW benchmark is a published integrated-light estimate, not a newly fitted background spectrum.
