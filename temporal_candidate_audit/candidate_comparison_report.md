# Temporal-field candidates compared with astronomical observations

The completed comparison retains every candidate as a distinct branch of a fictional, nonexpanding universe. Real observations are used as constraints on the equations, not as a requirement to abandon a world-building premise. Each branch is assessed as currently written; a proposed repair is kept separate from the original branch and must carry its own consequences.

The most useful empirical result is that whole-signal stretching remains consistent with supernova timing, while several distance laws remain almost indistinguishable in the nearby redshift sample. The new density-dependent rotation templates improve substantially on the fixed baryonic baseline, but are less accurate than the earlier acceleration-dependent template. New directional and endpoint-environment redshift extensions do not improve the reserved-sample root-mean-square error. Neither those results nor the numerical field calculation establish a universal law of time.

A connected one-dimensional field calculation now demonstrates simultaneous propagation redshift, pulse stretching and conservative energy exchange under specified artificial conditions. This closes a mathematical gap in the earlier separate-patch examples. Its dependence on epoch, radiation intensity and a finite packet-width prescription shows why it is not yet a cosmological model.

## 1. Evidence and experimental scope

The observational calculations use released, processed measurements. They do not reprocess telescope exposures, detector strain, satellite clock telemetry or Planck time streams. The redshift and rotation splits were examined earlier in this research sequence. Parameters for new variants are fitted on the original training objects only, but their subsequent evaluation is exploratory reuse, not a newly blind experiment. Reserved means excluded from that parameter fit, not previously unseen in the research program.

| Observation | Data used | Calculation completed |
|---|---|---|
| Nearby redshift | 164 Cosmicflows-4 SBF-selected groups; 104 training, 35 validation, 25 reserved | Six temporal/kinematic variants; original reference score retained |
| Supernova event widths | 1,504 DES event averages from previously exported released widths | Duration exponent and fixed-exponent residual comparisons |
| Spectral aging | 35 published supernova aging-rate measurements | Independent-method exponent and residual comparisons |
| Joint wave arrival | Published GW170817/GRB170817A timing summary | Conditional light-only and shared-coupling calculations |
| Local clocks | Published Galileo and atomic-frequency-ratio results | Constraints on explicitly defined extra couplings |
| Microwave spectrum | 43 FIRAS residual channels | Single-temperature and two-temperature mixture fits |
| Distant microwave temperature | Published HFLS3 interval at z=6.34 | Constant-temperature and coherent-scaling predictions |
| Angular microwave structure | 83 Planck TT and 66 EE binned measurements | Strictly homogeneous, initially unpolarized branch versus measured structure |
| Galaxy rotation | 149 eligible SPARC galaxies; 89 training, 29 validation, 31 reserved | Two new surface-density templates and prior dynamical benchmarks |
| Field dynamics | Synthetic periodic one-dimensional system | Connected propagation, finite packet energy exchange, resolution and epoch checks |

Cosmicflows-4, DES, spectral-aging and SPARC products inherit distance calibration, source modeling and selection assumptions. The catalog-derived galaxy distances are not assumed to be geometric distances without a flux law. Planck and FIRAS products inherit instrument calibration and foreground subtraction. These dependencies remain part of every result.[^1][^2][^3][^4][^5][^6]

## 2. Candidate formulas and observable meaning

Let S=1+z, R be geometric path length, D_L the brightness-based distance, and κ the fitted inverse-distance rate. Three phenomenological distance branches are retained:

\[
\mathrm{R1}:\quad S=e^{\kappa R},\qquad D_L=RS;
\]
\[
\mathrm{R2}:\quad S=1+\kappa R,\qquad D_L=RS;
\]
\[
\mathrm{R3}:\quad S=e^{\kappa R},\qquad D_L=R\sqrt S.
\]

R1 describes full-signal stretching with the previously assumed flux law. R2 is the homogeneous rolling-history shape with that same phenomenological flux law. R3 is the energy-loss-only branch, without event-duration stretching. The flux relations specify these branches; they have not been derived for the connected screened field.

For an observable temporal field, retain n=e^χ and the local optical Hamiltonian ω=c|k|e^−χ. Its propagation law gives ln S=∫∂tχ dt along the ray, subject to unchanged endpoint atomic standards and separate treatment of ordinary gravity and motion. A homogeneous χ=γt yields R2 in geometric distance when n at observation is normalized to one. An affine n produces R1 in the homogeneous limit but has the previously identified finite-past positivity problem. These are different temporal histories, not interchangeable derivations of the same equation.

A spatially static optical field can change travel times and ray directions without accumulating this propagation frequency shift. A common homogeneous lapse applied identically to all material and radiation behavior can be a coordinate redefinition. Neither mathematical observation excludes more general temporal theories; both require that a candidate specify an observable relation rather than only relabel time.

## 3. Nearby redshift: distance, direction and environment

The original source catalog contains surface-brightness-fluctuation distance measurements and CMB-frame redshifts. This analysis uses the preserved 164-group selection, with one representative per selected group and sky-tile-separated splits. For each model, the fit minimizes a Gaussian negative log likelihood on the 104 training groups, with 300 km/s added in quadrature to propagated distance uncertainty. The quoted errors below are redshift residuals multiplied by c, not interpreted peculiar velocities caused by the time model.[^1]

Three diagnostic extensions supplement R1. R4 changes the slope with direction, k(n̂)=k0 exp(a·n̂), adding three parameters. R5 changes the slope with standardized ln N, where N is the number of entries in the full Cosmicflows-4 individual-galaxy catalog sharing the source's group identifier. R6 adds a bulk-velocity term V·n̂ to the predicted cz while retaining R1. Here k=cκ is expressed in km/s/Mpc.

| Branch | Fitted parameters | Validation RMSE, km/s | Reserved RMSE, km/s |
|---|---:|---:|---:|
| R1 exponential, whole-signal | 1 | 435.93 | 412.44 |
| R2 rolling-history shape | 1 | 435.90 | 411.54 |
| R3 energy loss only | 1 | 436.06 | 413.47 |
| R4 directional rate | 4 | 454.82 | 419.44 |
| R5 source-group entry proxy | 2 | 475.51 | 460.11 |
| R6 added bulk velocity | 4 | 452.46 | 404.73 |
| Original expanding reference, fixed density parameters | 1 | 435.90 | 411.93 |

R1's fitted k is 75.598 km/s/Mpc; R2 gives 76.125 and R3 75.074. R2 improves the reserved RMSE relative to R1 by only 0.91 km/s. A descriptive sky-tile bootstrap interval for R2 minus R1 is −3.58 to +1.55 km/s. The corresponding interval for R3 is −1.45 to +3.73 km/s. Those data do not discriminate the three laws reliably.

R4 lowers training error but increases both validation and reserved error. Its reserved difference from R1 has a broad descriptive interval, −55.89 to +77.72 km/s. R5 also improves training fit but performs worse on both other subsets; its reserved difference is +47.67 km/s, with a descriptive interval of +3.79 to +95.51. The fitted R5 coefficient is +0.0296 per training-standard-deviation of log catalog membership. This is a property of the chosen catalog proxy and fit, not a measured environmental law.

R6's reserved error is lower, but its validation error is higher, and its descriptive reserved difference spans −83.59 to +64.62 km/s. It is a nuisance comparison, not a detection of a bulk flow or a temporal direction. These intervals resample existing sky tiles at fixed fitted parameters. They do not include training-fit uncertainty, a full peculiar-velocity covariance or the research program's model-selection history.

**Pros.** The simple redshift branches remain useful empirical descriptions. Directional and environmental extensions have explicit formulas that can be reproduced and improved without altering R1. The bulk-motion comparison illustrates that a directional residual need not originate in new timing physics.

**Cons.** No extension supplies a convincing predictive improvement here. Catalog membership is not physical mass density, excludes unseen group members, and inherits selection and group-assignment assumptions. Group identification itself uses redshift information, so this is not a redshift-independent environmental tracer. Most importantly, it measures a source-endpoint proxy, not the matter encountered throughout the photon trajectory. This is not the proposed full void-path test, and does not test screening by the foreground matter field.

**Outcome and goal.** Keep R1–R6 separately. A future environmental claim requires a selection-corrected spatial matter tracer, an explicit field solution and a derived flux law. An arbitrary direction fitted to these nearby galaxies is not an explanation of the CMB's large-angle alignment.

## 4. Event stretching: two observational methods

For light-curve widths, fit w=(1+z)^b. For spectral aging, fit a=(1+z)^−b. The common b=1 branch predicts that events appear longer and spectral evolution correspondingly slower. The b=0 branch predicts no temporal stretching. The two methods probe related physics with different measurements and different template assumptions.[^2][^3]

| Data | Objects | Fitted b, formal statistical error | χ² at b=1 | χ² at b=0 |
|---|---:|---:|---:|---:|
| DES event widths | 1,504 | 1.00288 ± 0.00483 | 2,166.82 | 30,117.37 |
| Spectral aging | 35 | 0.96644 ± 0.10402 | 26.95 | 150.57 |

These reproduce the previous derived-data calculations. The DES width fit includes correlated information through reference curves and shared calibration that is not represented in the diagonal uncertainties. The published width method used here is partly a consistency check because its reference construction includes assumed stretching. The spectral-aging table provides a different-method check, but its low-redshift template treatment also has assumptions. Neither comparison is raw-data or assumption-free evidence.

**Pros.** Whole-signal temporal branches reproduce a necessary observational property without requiring expansion as an input to their propagation equation. The independent-method spectral-aging result is useful corroboration of the observable stretching target.

**Cons.** Energy-loss-only R3 gives much worse timing residuals as written. It can be retained with a separate source-evolution or propagation-duration extension, but that extension must predict the observed relation rather than assigning a correction to each object. Branches sharing b=1 remain indistinguishable in this test. The large χ² differences are diagnostic residual comparisons, not advertised extreme-sigma probabilities from a complete likelihood.

**Outcome and goal.** Signal stretching should remain a primary observable in any new time formula. Agreement with b=1 is not evidence selecting the field mechanism over every other explanation.

## 5. Light, gravitational waves and material clocks

The published gamma-ray minus gravitational-wave arrival delay for GW170817/GRB170817A is 1.74±0.05 seconds. A source-emission delay is part of the interpretation. Instead of importing the standard constant-speed bound into the evolving-field model unchanged, calculate its travel times explicitly.[^7]

For nγ=e^(γt), nγ(0)=1 and photon arrival at t=0, the photon flight time is ln(1+γR/c)/γ. If gravitational waves travel at constant c and were emitted simultaneously, their arrival is later by R/c−ln(1+γR/c)/γ. The table uses the original γ=7.7314966×10⁻¹¹ per year. Distances of 30–50 Mpc are illustrative sensitivity cases; no gravitational-wave distance has been re-inferred in modified gravity.

| Distance | Light-only coupling: gravity arrives later | Shared coupling: propagation difference |
|---|---:|---:|
| 30 Mpc | 368,251 years | 0 seconds |
| 40 Mpc | 653,578 years | 0 seconds |
| 50 Mpc | 1,019,516 years | 0 seconds |

For a generalized nGW=e^(pγt), exact shared propagation is p=1. Around that value, a chosen 10-second propagation-delay allowance at 40 Mpc corresponds to |p−1| of approximately 4.83×10⁻¹³. This is a conditional sensitivity estimate under this particular time history and normalization, not a measured confidence bound. The 10 seconds is an illustrative emission/propagation allowance, not the instrumental timing error.

**Pros.** Shared light–gravity propagation supplies a clear way for a temporal theory to preserve near-coincident arrivals. The light-only branch now has a quantitative consequence rather than a qualitative concern.

**Cons.** The light-only homogeneous branch has an enormous discrepancy under the stated assumptions. Source timing, a different history, or spatial transitions could be explored as distinct extensions; ordinary seconds-scale source delays cannot cancel hundreds of thousands of years. Shared propagation is not a demonstrated covariant gravity theory and does not by itself establish how material clocks remain different.

Local clock observations constrain additional couplings. Galileo's published gravitational-redshift deviation is (0.19±2.48)×10⁻⁵. Its approximate Gaussian 95% interval is −4.67×10⁻⁵ to +5.05×10⁻⁵. That bounds deviations from the ordinary gravitational-redshift coefficient within the measurement model; it is not a direct bound on all field amplitudes.[^8]

Atomic-frequency comparisons give an additional conditional check. Define explicitly α̇/α=q_eff γ. Using the published (1.0±1.1)×10⁻¹⁸ per year summary yields a conservative approximate |q_eff| requirement below 4.08×10⁻⁸.[^9]

| Intrinsic α coupling q | Screening fraction | Predicted α̇/α, per year | Interpretation under this coupling |
|---:|---:|---:|---|
| 0 | 1 | 0 | Compatible; does not test common clock rescaling |
| 1 | 1 | 7.73×10⁻¹¹ | Much larger than the observed drift scale |
| 1 | 10⁻⁶ | 7.73×10⁻¹⁷ | Approximately 69 quoted standard errors from the central value |
| 1 | 10⁻⁸ | 7.73×10⁻¹⁹ | Within the quoted uncertainty |

A millionfold screening example is therefore not automatically adequate if the intrinsic α coupling is unity. Conversely, a field can leave α unchanged through compensating changes of constitutive quantities; q=1 must not be assumed merely because n changes. The common-clock and atomic-ratio questions are distinct.

**Outcome and goal.** Retain light-only, shared-wave and screened-clock branches. The next common requirement is a specified action or equivalent equations for matter, photons and gravitational waves, from which these couplings follow. Choosing them separately is presently additional freedom.

## 6. CMB spectrum, temperature and angular structure

The spectrum test fits the released FIRAS residual table about its reference 2.725 K Planck spectrum. The model adds a temperature correction and the supplied Galactic template, and optionally averages two equal-weight Planck spectra at T(1−δ) and T(1+δ). Quoted channel errors are treated diagonally. This deliberately restricted mixture represents differing *received* thermal temperatures, not every possible emission or thermalization process.[^5]

| Received temperature split | Approximate χ², 43 channels |
|---|---:|
| One temperature | 45.02 |
| ±0.001% | 45.02 |
| ±0.1% | 45.03 |
| ±0.3% | 46.51 |
| ±1% | 235.84 |
| ±5% | 117,314.36 |

The fitted nuisance temperature near 2.7250005 K is relative to the released residual convention. It is not a new absolute temperature measurement and should not be contrasted with the independently calibrated 2.72548 K result as a disagreement. Missing covariance and the restricted mixture prevent interpreting these values as a calibrated general exclusion bound.

**Pros.** Uniform, frequency-independent scaling can preserve an existing thermal occupation distribution. Small differences between received temperatures can evade this approximate monopole-spectrum test. This allows coherent temporal cooling as a branch in a fixed-size universe.

**Cons.** The model does not create thermal radiation from arbitrary starlight. Mixtures separated by percent-level temperatures produce large spectral discrepancies here. Thermalization, a cavity, or energy recycling can be retained as extensions, but need opacity, emissivity, equilibration-time and polarization predictions. Arbitrary source spectra can fit a spectrum without identifying a physical mechanism.

A separate temperature-history comparison uses HFLS3 at z=6.34. Its published molecular-excitation analysis gives 16.4–30.2 K at one sigma. A branch with locally constant 2.72548 K lies outside that range; the relation T(z)=T0(1+z) predicts 20.005 K, inside it. This is a published-summary comparison dependent on the source excitation and atomic/molecular physics, not a new fit of the absorption spectrum.[^10][^11]

The optical model also needs correct mode counting: with ω=ck/n in a fixed volume, photon energy density for a Planck occupation distribution scales as n³T⁴, rather than using the vacuum coefficient unchanged. Under homogeneous adiabatic evolution T∝1/n, this gives uγ∝1/n. Entering a screened detector environment can change the mapping to observed intensity. Therefore neither the simple temperature relation nor the spectrum fit establishes a global screened CMB solution.

### Angular power and polarization

Newly downloaded official Planck products contain 83 binned temperature-power points over effective multipoles 47.7–2499.0 and 66 E-polarization points over 47.7–1988.0. They are processed bandpowers, not sky-map pixels. A strictly homogeneous scalar cooling an initially uniform, unpolarized sky predicts zero TT and EE fluctuations. The following checks quantify the additional structure that branch would need.[^6]

| Quantity | TT temperature | EE polarization |
|---|---:|---:|
| Number of bins | 83 | 66 |
| Bins exceeding zero by five quoted errors | 81 | 46 |
| Zero-signal diagonal χ² diagnostic | 715,128.92 | 54,424.97 |
| Released best-fit reference diagnostic | 65.12 | 62.90 |

The reference column is the best-fit standard-model curve supplied with the data; it is not a newly fitted or blind comparison. Quoted bandpower errors include assumptions about sample variance and processing. These numbers should not be turned into extreme-tail significance for the zero-signal branch. Their practical meaning is that angular fluctuations and polarization are substantial measured targets that a homogeneous temperature law leaves unexplained.

These binned products start well above the quadrupole and octopole and discard orientation information. They cannot test the so-called axis of evil. The Planck large-scale analysis reports temperature anomalies but does not supply an unambiguous matching polarization detection.[^12] In the fictional universe the alignment can be stipulated; the formula still needs a direction, amplitude, odd/even multipole structure and independent consequences before it can be scored as an explanation. A perfectly symmetric sphere selects no axis by itself. No preferred-axis fit or topology constraint is claimed here.

**Outcome and goal.** Preserve coherent cooling, recycling, cavity and directional branches. Only the restricted spectrum, temperature-history and homogeneous-angular limits are numerically tested here. A complete CMB candidate must predict the received spectrum, anisotropy and polarization together.

## 7. Galaxy rotation: density versus acceleration dependence

The SPARC calculation keeps the original quality, inclination and valid-row cuts, yielding 149 galaxies. Stellar mass-to-light ratios are fixed at 0.5 for disks and 0.7 for bulges, and published distances and inclinations are retained. Fits minimize the mean over galaxies of the mean squared logarithmic speed residual, so densely sampled galaxies do not dominate merely through point count. The reserved set has 31 galaxies and 642 radial measurements.[^4]

The earlier acceleration-dependent branch is

\[
g_{\rm extra}=A\,a_K(g_{\rm bar}/a_K)^p,\qquad
v_{\rm pred}^2=r(g_{\rm bar}+g_{\rm extra}),\quad a_K=c^2\kappa.
\]

Its preserved parameters are A=10⁻⁰·⁶¹⁸⁶⁹⁹ and p=0.462459. The two new density-proxy branches instead use

\[
v_{\rm pred}^2=v_{\rm bar}^2
\left[1+\frac{B}{1+(\Sigma_\star/\Sigma_*)^q}\right],
\]

where Σstar=0.5 SBdisk+0.7 SBbulge in solar masses per square parsec. This is a stellar surface-density proxy, not total volume density; gas density and external screening are absent. One branch fixes q=1 and fits B and Σ*. The other also fits q. Parameters are trained on the original 89 galaxies, with no per-galaxy tuning.

| Branch | Reserved log-speed RMSE, dex | Reserved speed RMSE, km/s | Non-Hubble-distance subset RMSE, km/s |
|---|---:|---:|---:|
| Fixed baryonic baseline | 0.25608 | 47.77 | 45.37 |
| Density proxy, q=1 | 0.09339 | 29.25 | 32.95 |
| Density proxy, fitted q | 0.09309 | 28.19 | 31.81 |
| Earlier acceleration power | 0.07914 | 17.20 | 14.98 |
| RAR functional benchmark | 0.07840 | 16.40 | 13.11 |

The density fits have log10 B=0.54884 and log10 Σ*=1.23373 for q=1; the three-parameter version gives log10 B=0.59570, log10 Σ*=1.00714 and q=0.70119. The non-Hubble-distance subset has 11 galaxies and 248 measurements, using the catalog's TRGB, Cepheid or supernova distance-method categories. These are less directly dependent on an assumed Hubble law, not completely free of calibration assumptions.

For the fitted-q density branch minus the acceleration-power branch, the galaxy-resampled descriptive 95% interval for log-RMSE difference is +0.00255 to +0.02678 dex. The power-minus-RAR interval is −0.00248 to +0.00405 dex. Thus density suppression is useful compared with the baryonic baseline but is less successful than acceleration dependence within this restricted comparison; power and RAR remain closely matched. Distance, inclination, mass-to-light and correlated observational uncertainties are not marginalized here.

**Pros.** A simple suppression function captures a meaningful part of the observed pattern, and the improvement persists on galaxies excluded from its parameter fit. Its bounded response is easy to interpret and test. The acceleration branch provides a more accurate empirical target for a future temporal-gravity completion.

**Cons.** These are dynamical enhancement templates. They are not predictions of the optical screened-field equations and do not demonstrate that a clock illusion explains rotation. A common wavelength multiplier on the approaching and receding sides cancels from the normalized splitting, (λrec−λapp)/(λrec+λapp), under the usual symmetric Doppler construction. A pure apparent-time explanation must derive the modified Doppler law and motion consistently rather than multiplying tabulated rotation speeds after the fact.

The apparent use of the redshift scale aK does not yet unify the two datasets: gextra=A aK^(1−p) gbar^p allows a change in κ to be absorbed into the free amplitude A. A genuine cross-scale prediction must remove that degeneracy or independently fix A. No new lensing prediction has been derived.

The previously explored gravity-cliff branch remains archived separately. It improved the earlier 42-galaxy cross-validation score, but its cliff coefficient hit the allowed upper bound in every fold. That boundary behavior and the different sample prevent treating it as an established improvement over the new 31-galaxy comparison. It remains an empirical extension requiring a physical force and clock law.

## 8. Connected temporal fields and conservative energy exchange

This section is a new mathematical experiment, separate from the real-data fits. It evolves two fields across a periodic one-dimensional box containing smooth dense regions around the boundary and an intervening void. Positive quadratic field potentials couple χ to an energy-bearing field ψ, while a density-dependent positive term suppresses χ in dense regions. Both fields have positive gradient and kinetic energies. The mass-weighted lattice contains 81, 161 or 321 sites per field.

In dimensionless units, take ε=10⁻⁴, coupling frequency 10, maximum density-induced frequency 100, and initial ψ rate 0.0077314966. Initial χ velocity follows the local equilibrium suppression fraction; initial field values vanish. These choices are illustrative, not fitted screening parameters. The spatial gradients mean that the separate homogeneous-patch solution is not simply pasted across the box.

A finite-width ray surrogate has Hamiltonian Hγ=p exp(−χbar), where χbar is a Gaussian-weighted field average with width 0.02 box lengths. The derivative of the same Hamiltonian supplies the field's energy-exchange source and the ray's momentum equation. Matched interpolation and deposition conserve the lattice field-plus-ray energy. The finite width is a coarse-graining assumption requiring a microscopic justification; it is not a measured physical photon size.

| Grid sites per field | Negligible-backreaction frequency stretch | Pulse stretch |
|---|---:|---:|
| 81 | 1.003115627 | 1.003115617 |
| 161 | 1.003086497 | 1.003086497 |
| 321 | 1.003080703 | 1.003080705 |

At 321 sites, the pulse and frequency factors agree to approximately 2.0×10⁻⁹ relatively. The redshift itself changes by approximately 0.19% between the 161- and 321-site calculations, so the digits are solver outputs rather than physical precision claims. Frequency-independent propagation and event stretching coexist across the connected inhomogeneous model for this launch epoch.

With a ray packet initially carrying 0.1% of the field energy, the 321-site stretch is 1.00335034. At 10% it is 1.02987139. The corresponding photon-energy changes are −9.98×10⁻⁷ and −8.67×10⁻⁴ in dimensionless energy units. Total relative energy error remains below 10⁻¹² in these runs. The redshift at the larger energy ratio changes by about 2.1% between 161 and 321 sites, so backreaction convergence is less complete than the negligible-backreaction comparison.

An earlier numerical branch tied packet width to grid spacing; refining the grid therefore changed the physical source, and its backreaction changed substantially. Its code and outputs are retained in the archive. The fixed-width branch above separates grid refinement from changing the source prescription. This correction is a numerical-model clarification, not observational evidence favoring a particular packet width.

**Pros.** The connected calculation supplies an explicit conservative place for a ray's lost energy to go. It no longer relies on an externally prescribed field history alone. The negligible-backreaction limit demonstrates the matching signal-stretch observable required by the supernova comparisons.

**Cons.** The field's behavior depends on epoch. For the 161-site model, launching at dimensionless times 0, 10, 20 and 40 gives S=1.00308650, 1.00276697, 1.00073333 and 0.99647533. The last is a blueshift. Positive minimum mode frequency squared is consistent with stable oscillatory behavior, not a perpetual monotonic redshift source. Density is fixed, gravity and atomic clocks are absent, and the ray is not a complete electromagnetic field. Radiation intensity and the undetermined physical field-energy scale can change the shift, creating an additional universality test.

**Outcome and goal.** The finite-epoch energy-conserving mechanism is mathematically feasible. Its amplitude, long-term history, packet prescription and full matter–radiation–gravity completion remain uncalibrated. No claim that the connected model fits the astronomical redshift sample, the CMB or galaxy rotation is warranted yet.

## 9. Retained candidate scorecard

| Candidate | Goal | Observed or computed outcome | Development status |
|---|---|---|---|
| R1 exponential distance | Predict nearby redshift | 412.44 km/s reserved RMSE | Retain as empirical reference |
| R2 rolling history | Obtain redshift from positive evolving timing factor | 411.54 km/s; difference from R1 inconclusive | Retain; derive connected flux law |
| R3 energy loss only | Explain spectral redshift without event stretching | Similar nearby redshift score; larger supernova timing residuals | Retain with explicit duration-extension requirement |
| R4 directional rate | Predict an anisotropic redshift effect | No reserved-score improvement | Retain as exploratory direction branch |
| R5 endpoint environment | Connect timing to source environment | Worse validation and reserved scores for catalog-membership proxy | Retain; distinguish from foreground-path screening |
| Shared wave propagation | Preserve light–gravity arrival coincidence | Equal coupling removes propagation difference by construction | Derive coupling from gravity theory |
| Screened material response | Preserve local clock measurements | Requirements depend on atomic coupling; 10⁻⁶ suppression insufficient for q=1 example | Retain coupling-specific constraints |
| Coherent CMB cooling | Preserve thermal shape and permit higher remote temperature | Restricted spectrum compatible; 20.005 K inside HFLS3 interval | Origin, intensity and boundaries unresolved |
| Thermal mixture/recycling | Generate or maintain a microwave bath | Percent-level received-temperature mixtures fit poorly | Specify thermalization physics |
| Directional/topological CMB | Explain alignment and angular structure | Uniform branch omits TT/EE; no axis prediction scored | Need angular and polarization equations |
| Density-dependent rotation | Explain outer-galaxy enhancement | Improves baseline; less accurate than acceleration power | Empirical branch, not derived time theory |
| Acceleration/cliff rotation | Find a gravity-linked temporal law | Power matches RAR closely; earlier cliff fit boundary-limited | Need unique cross-scale normalization and lensing |
| Connected two-field energy model | Combine screening, redshift and conservation | Finite-epoch example succeeds; later shifts can reverse | Retain with history and backreaction dependencies |

The observational assessment is complete for the specified formulas and restricted comparisons above. A full foreground-path screening fit, an axis-of-evil map likelihood, a modified-gravity lensing calculation and a complete CMB-generation model are not numerical results hidden elsewhere in this report: the candidate equations do not yet define them. Their missing inputs are stated explicitly so future extensions can be compared against the same scorecard.

## 10. Reproduction and data provenance

The accompanying archive contains this report, all input tables needed for the calculations, published-summary values with source URLs, scripts, per-object predictions, numerical results, file hashes and verification checks. Run `analyze.py`, then `uncertainty_checks.py`, and `connected_field.py` with Python, NumPy, SciPy and pandas. Set `OPENBLAS_NUM_THREADS=1` for predictable small-matrix performance. No network access is required after extraction. `prepare_inputs.py` documents acquisition and staging from the original session; it is not required for reproduction from the bundled data.

`results/redshift_predictions.csv` includes all candidate predictions and the endpoint catalog proxy. `results/SPARC_predictions.csv` contains all five models' reserved-galaxy predictions and stellar surface-density values. The Planck comparison CSVs contain bandpowers, quoted errors and the released reference curve. The archive also retains the original SPARC split and parameters, historical cliff results and the expanding-reference redshift scores. Original labels and old frozen objects were not overwritten.

Verification checks confirm the 164 unique selected groups, their 104/35/25 split, sky-tile and group separation, SPARC galaxy separation, finite predictions and successful optimizations. Energy conservation is checked independently of the frequency-stretch comparison. Spatial convergence is reported rather than inferred from solver tolerance alone. The source products remain processed measurements with the limitations stated above.

## Sources

[^1]: Tully, R. B., et al. *Cosmicflows-4* (2023). https://arxiv.org/abs/2209.11238 . Catalog J/ApJ/944/94, table2 and ReadMe: https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJ/944/94 . The selected SBF group data, calibrations and prior split are preserved in the archive.
[^2]: White, R. M. T., et al. *The Dark Energy Survey Supernova Program: Slow supernovae show cosmological time dilation out to z~1* (2024). https://arxiv.org/abs/2406.05050 . Released derived widths: https://github.com/ryanwhite1/DES-Time-Dilation . This report reuses event averages previously exported from those releases.
[^3]: Blondin, S., et al. *Time Dilation in Type Ia Supernova Spectra at High Redshift* (2008), Table 3. https://arxiv.org/abs/0804.3595 . The table includes 22 low-redshift and 13 high-redshift objects.
[^4]: Lelli, F., McGaugh, S. S., and Schombert, J. M. *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves* (2016). https://arxiv.org/abs/1606.09251 . Data release: https://astroweb.case.edu/SPARC/ . Source files SPARC_Lelli2016c.mrt and Rotmod_LTG.zip are bundled.
[^5]: Fixsen, D. J., et al. *The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set* (1996). https://arxiv.org/abs/astro-ph/9605054 . NASA LAMBDA product documentation: https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html . Released residual table: https://lambda.gsfc.nasa.gov/data/cobe/firas/monopole_spec/firas_monopole_spec_v1.txt .
[^6]: Planck Collaboration. *Planck 2018 results. V. CMB power spectra and likelihoods* (2020). https://www.aanda.org/articles/aa/full_html/2020/09/aa36386-19/aa36386-19.html . Official archive: https://pla.esac.esa.int/ . Products COM_PowerSpect_CMB-TT-binned_R3.01.txt and COM_PowerSpect_CMB-EE-binned_R3.02.txt; exact retrieval URLs and hashes appear in input_manifest.json. NASA archive documentation: https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/ .
[^7]: LIGO Scientific Collaboration, Virgo Collaboration, Fermi GBM and INTEGRAL. *Gravitational Waves and Gamma-rays from a Binary Neutron Star Merger: GW170817 and GRB170817A* (2017). https://arxiv.org/abs/1710.05834 . Arrival-delay summary used; detector strain and waveform inference not reprocessed.
[^8]: Delva, P., et al. *Gravitational Redshift Test Using Eccentric Galileo Satellites* (2018). https://doi.org/10.1103/PhysRevLett.121.231101 . Published coefficient and uncertainty used.
[^9]: Lange, R., et al. *Improved Limits for Violations of Local Position Invariance from Atomic Clock Comparisons* (2021). https://doi.org/10.1103/PhysRevLett.126.011102 . Published fine-structure-constant drift summary used.
[^10]: Riechers, D. A., et al. *Microwave Background Temperature at a Redshift of 6.34 from H2O Absorption* (2022). https://arxiv.org/abs/2202.00693 . Published one-sigma excitation-analysis interval used.
[^11]: Fixsen, D. J. *The Temperature of the Cosmic Microwave Background* (2009). https://arxiv.org/abs/0911.1955 . Calibrated temperature used for the remote-temperature comparison, not re-estimated from residuals here.
[^12]: Planck Collaboration. *Planck 2018 results. VII. Isotropy and Statistics of the CMB* (2020). https://arxiv.org/abs/1906.02552 . Published large-angle interpretation; no new orientation-sensitive map reanalysis.
