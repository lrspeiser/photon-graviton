# First-principles temporal physics: hypotheses, tests, and milestones

## Research objective

Develop a physical law that accounts for clock comparisons in gravitational fields, spectral redshift and event-duration stretching in fixed material geometry. The working setting does not assume expansion, a Big Bang, or particulate dark matter. Those premises do not excuse a failure to reproduce the observations: the model must specify material standards, photon propagation, energy exchange, and gravitational effects consistently.

The most useful current results are restrictive. A spatially homogeneous change of the time coordinate cannot create observable cosmological redshift by itself. Ordinary stationary gravitational clock differences do not accumulate through intervening wells. A dynamically evolving photon response can stretch signals, but a passive field pinned around matter does not automatically sustain the required drift. Thermal preservation is possible conditionally; thermal generation and global consistency remain unestablished.

The next milestone is to construct a driven or evolving field with explicit energy accounting and derive its effects on material clocks and light. A successful calculation would establish internal consistency. Evidence identifying that field would require additional predictions that distinguish it from other mechanisms.

## Evidence inventory and limitations

The prior study used 164 Cosmicflows-4 galaxy-group representatives with surface-brightness-fluctuation distances: 104 train, 35 validation, and 25 test. The fit κ=7.7314966×10⁻⁵ per million light-years remains an empirical calibration under stated brightness assumptions. The same test objects have already been examined; all subsequent comparisons on them are exploratory. Distance-calibration and peculiar-motion assumptions remain inherited.

The previous extension reproduced event-duration fits from 1,504 DES supernovae using released per-band widths. Their reference curves include dilation corrections; that reproduction is a consistency check, not independent raw-photometry evidence. The exponent was b=1.002877±0.004825 with formal statistical errors; the original publication gives an additional systematic uncertainty of 0.010.[1]

This report adds a different observation method: Table 3 of Blondin et al. (2008), containing spectral-aging rates for 35 supernovae. All table entries were extracted and checked against the published table. These are derived spectral measurements, not raw spectra. The low-redshift template ages include a small dilation correction; the authors examine its influence. No claim of entirely assumption-free evidence is made.[2]

The microwave calculations use the 43-channel NASA LAMBDA FIRAS residual product, individual channel uncertainties and its Galactic template. It is processed data, including calibration and foreground assumptions. The supplied diagonal errors do not replace a full covariance analysis.[3] Gravitational-clock and distant-temperature checks use published summary measurements rather than new reductions of experimental data.[4–6]

## Hypothesis ledger

| ID | Hypothesis | Goal | Outcome that would support development | Current result and decision |
|---|---|---|---|---|
| H1 | A homogeneous universal clock-rate change causes redshift | Establish a measurable effect independent of coordinate choice | A dimensionless source-to-detector difference that survives consistent clock calibration | The pure lapse-only model is a coordinate reparameterization. Close this restricted branch. |
| H2 | Static gravity-well differences accumulate with traveled distance | Recover local clock tests and obtain a path-dependent redshift | Nonzero accumulated shift after consistent endpoint comparisons | Stationary metric propagation yields endpoint shifts, not accumulation. Keep ordinary gravity as a baseline; close the accumulation claim under these assumptions. |
| H3 | Matter-screened, time-dependent photon response | Produce redshift and whole-signal stretching with stable endpoints | Derive screening from matter and reproduce both observables with one parameter set | Prescribed ray construction works. An undriven pinned field does not sustain linear drift. Continue only with specified dynamics and energy source. |
| H4 | Light loses energy without whole-signal stretching | Explain redshift through energy transfer alone | Spectral-aging rates consistent with no stretching, or independently predicted intrinsic evolution | New spectral-aging reanalysis disfavors the no-stretch prediction under the template assumptions. Do not retain energy loss alone as the mechanism. |
| H5 | The field retains environmental history | Generate a sustained, predictable response from finite memory | One relaxation timescale predicts independent epoch/environment residuals | Constant input generally approaches equilibrium; memory alone does not imply permanent drift. Timescale and driver are not identified by existing data. |
| H6 | Coherent temporal cooling preserves a thermal bath | Match CMB spectral shape and remote excitation temperatures in fixed space | One transport model preserves the spectrum through screened paths and predicts T(z) | Homogeneous thermal preservation is conditional; present spectrum cannot identify the mechanism or its rate. Global screened solution remains missing. |
| H7 | Mixed sources, recycling or a finite cavity generate the CMB | Generate thermal radiation without assuming an initially thermal bath | Explicit equilibration yields the spectrum, photon abundance and energy budget | A tested ±1% two-temperature mixture is strongly disfavored by the approximate spectral fit. Origins requiring thermalization remain open, not confirmed. |
| H8 | Distance and void exposure modify the response | Derive the distance law from measured local conditions | Independent path-density inputs improve frozen predictions beyond calibration and velocity effects | Uniform active fraction is exactly degenerate with field strength. No independent environment reconstruction is in the current dataset. |
| H9 | Global direction or topology produces the large-scale alignment | Predict a preferred axis rather than choose it from the target | A direction fixed elsewhere predicts CMB temperature and polarization consistently | Arbitrary harmonic alignment is not an explanation. No specified model currently supports a meaningful new map likelihood. |
| H10 | The same field explains gravitational clock effects and orbital anomalies | Derive clocks, forces and light bending together | One action passes local tests and predicts orbital and lensing observables without separate rescaling | Local clock departures are tightly constrained. A unified force and lensing model has not been supplied. |

A hoped-for outcome is not a required result. A branch should be rejected when its declared assumptions fail a test; new assumptions create a new branch and must be reported explicitly.

## H1–H2: distinguish physical time from coordinate time

For a spatially homogeneous lapse and fixed spatial metric,

    ds² = −N(t)² c² dt² + dx² + dy² + dz²,
    τ = integral N(t) dt.

The transformation to τ gives a static Minkowski metric. Stationary identical clocks and photons sharing this metric do not acquire a cosmological redshift from N(t) alone. This is an exact statement for this restricted ansatz, not a rejection of all universal scalar or time theories. A scalar that also changes spatial geometry, mass ratios, interactions or gradients is a different physical proposal and must be evaluated on those observables.[7]

For a static position-dependent lapse, stationary gravitational frequency comparisons depend on the endpoint values. A photon can change locally while traversing a well, but intermediate stationary factors cancel in the net comparison. A source and detector with equal lapse do not acquire an extra redshift simply because the path crosses more static wells. Time dependence, moving interactions or explicitly nonmetric transport must be introduced to change that conclusion.[7]

The practical goal is to define every prediction as a dimensionless observable: frequency ratios, interval ratios, angular measurements or quantities compared to explicit material rulers. Distance from an observer cannot be an intrinsic local state unless the model identifies a physical reference structure. Path length can enter through integration of a local response; distance from matter can enter through a sourced field.

## H3: a new consistency obstruction for passive screening

The prescribed screened model previously studied has n=1+γt f(x), where f vanishes at source and detector. Its ray equations give ln S=(γ/c) integral f dx and the same factor for pulse stretching. The numerical agreement establishes the characteristic equations, not that matter generates this response.

Consider the simplest undriven linear field around a stable minimum:

    χ_tt − vχ² χ_xx + m²(x)χ = 0,
    m²(x) ≥ 0, with χ=0 at both ends of a finite interval.

A proposed secular solution χ=t f(x) requires −vχ² f''+m²f=0. Multiplying by f and integrating yields

    integral [vχ²(f')² + m² f²] dx = 0.

The endpoint term vanishes. Both integrands are nonnegative, so the only such solution is f=0. This rules out a nonzero indefinitely linear profile in this finite, undriven, positive-operator model. It does not rule out a rolling nonlinear potential, an external or dynamical driver, changing matter, different boundaries or an additional field. Those additions must specify energy supply and stability.

A numerical illustration used the lowest massless standing mode with pinned endpoints, positive optical response, and photons propagated through the resulting time-dependent medium. The normalization corresponds to κR=0.0077315 for R=100 Mly; the sine spatial profile has average 2/π. If its initial drift were held fixed, the expected z would be about +0.004934. With the field propagation speed equal to c and initial phase zero, the evolving mode instead gives z=−0.000007465: leading-order accumulation cancels as the field turns over during transit. Other phases can give blueshift. Slower modes can remain nearly monotonic over one flight, but still require an explanation for their rate and phase over cosmic times.

These are synthetic field and ray calculations, not astronomical data or a fit to galaxies. They identify a concrete weakness in the proposal that matter pinning alone explains secular time evolution. Screening theories in the literature likewise have nontrivial dynamical and laboratory restrictions; their results cannot simply be imported as validation of this different model.[8]

Development goal: choose one action or explicit coupled dynamical system, solve for its response without prescribing f(x), and show that energy transfer sustains the required evolution with controlled backreaction. The driver may be an evolving reservoir field, but its energy must not disappear from the accounting.

## H4: independent-method spectral-aging test

For whole-signal stretching, the apparent spectral-aging rate is a(z)=(1+z)⁻¹. Without stretching, comparable intrinsic evolution gives a(z)=1. We fitted the general exponent a(z)=(1+z)⁻ᵇ to the 35 published rates using their quoted Gaussian aging-rate errors.[2]

| Sample | Objects | Fitted b, formal 1σ | χ² for b=1 | χ² for b=0 |
|---|---:|---:|---:|---:|
| All spectral-aging measurements | 35 | 0.9664 ± 0.1040 | 26.95 | 150.57 |
| High-redshift subset | 13 | 0.9506 ± 0.1038 | 3.65 | 123.69 |

The fixed b=1 and b=0 comparisons have the same number of fitted parameters: zero. The high-redshift subset is part of the full sample, not independent additional evidence. The results reproduce the published rounded conclusions. Both expansion and our temporal mechanism predict b=1, so these measurements distinguish stretching from energy-only loss, not the cause of stretching. No extreme-significance claim is inferred from these χ² values; template errors and correlations matter.

This is an improvement over relying entirely on DES widths because spectral aging measures changes in spectral features rather than photometric event width. It is not entirely independent of low-redshift template calibration. The appropriate retained requirement is that the final field model reproduce spectral evolution and pulse intervals together without freely adjusting intrinsic evolution with redshift.

## H5: finite memory and the energy reservoir

A simple memory law, χ_dot+χ/τ=J, tends toward χ=Jτ for a constant driver. After equilibration χ_dot approaches zero. This means that memory can store a response, but does not by itself create permanent redshifting through continuing temporal evolution. For intervals much shorter than τ, a response can mimic a linear drift. That observational approximation cannot determine a unique τ.

The next test must specify J in terms of independently constrained matter or radiation and include the reservoir equation. A schematic pair of energy balances is ∂t uγ+div Fγ=−Q and ∂t uχ+div Fχ=+Q. The previous present-day CMB-based estimate Q≈1.02×10⁻³¹ W/m³ assumes all radiation couples at the fitted present fractional rate; it is not a measured energy density of the new field. Increasing field inertia can reduce backreaction but can also increase stored kinetic energy and gravitational stress.

Development goal: derive the evolution of the effective redshift coefficient, then test its curvature or epoch drift with frozen parameters. A memory model must predict a difference from the memory-free model larger than calibrated measurement uncertainty. Existing nearby distance fits do not do this: energy-only, full-signal and alternative positive-history prescriptions differ by less than 2 km/s in test RMSE against scatter near 412 km/s.

## H6–H7: coherent cooling versus mixing thermal histories

In the homogeneous adiabatic optical model, conserved occupation at fixed wave vector and ω proportional to 1/n preserve a Planck distribution when T is proportional to 1/n. The variable mode density gives photon density proportional to n³T³ and radiation energy proportional to n³T⁴. Fixed-volume photon number is then constant and radiation energy decreases as 1/n. An initially thermal bath can cool without increasing material distances. This does not explain how the bath became thermal or determine its present temperature.

Spatial screening complicates the argument. If a detector combines paths with different effective temperatures, the result is generally not exactly Planckian. We therefore tested a concrete family directly against the FIRAS residuals: equal weights of blackbodies at T(1+δ) and T(1−δ), with mean temperature and residual Galactic-template amplitude fitted as nuisance parameters.[3]

| Fractional temperature contrast δ | χ² | Increase above best mixture fit |
|---|---:|---:|
| 0, a single temperature | 45.02 | 0.00066 |
| 0.00001 | 45.02 | 0.00066 |
| 0.001, ±0.1% | 45.03 | 0.0127 |
| 0.003, ±0.3% | 46.51 | 1.49 |
| 0.01, ±1% | 235.84 | 190.82 |
| 0.05, ±5% | 117314.36 | 117269.34 |

For a fixed δ, there are 43 channels and two nuisance parameters. The best-fit contrast is approximately 0.000431, but its improvement over zero is negligible. A diagnostic profile threshold Δχ²=3.84 occurs at δ≈0.003787. This is NOT a calibrated 95% physical bound: the fit omits full covariance and calibration uncertainty, restricts the source mixture, and has a nonnegative contrast parameter. It nevertheless identifies substantial spectral curvature for percent-scale mixing in this model.

The result limits a specific mixture hypothesis, not every recycled-radiation model. Arbitrary frequency-dependent emissivity, absorption or genuine thermalization changes the prediction and requires separate data. A 10⁻⁵ directional temperature variation is far too small for this simple monopole test to select or reject the proposed sky alignment. Different source temperatures may also yield the same received temperature if source evolution and transfer compensate; the calculation concerns the received distribution.

A published water-absorption measurement at z=6.34 gives 16.4–30.2 K at one standard deviation. T0(1+z)=20.005 K for T0=2.72548 K lies inside that interval.[5,6] This is compatible with an evolving bath in fixed space but does not identify its dynamics. Molecular excitation must be recomputed if the field changes local atomic physics. A complete thermal milestone requires a collision model or an explicit initial thermal state, photon-number accounting, screened transport and the observed angular and polarization properties.

## H8–H9: environment, direction and geometry

The current redshift formula measures only the product of local response strength and average active path fraction. No statistical procedure can separate two exactly degenerate parameters without an independent measurement. The prior 8.74% slope change is also equivalent in magnitude to a 0.182 mag distance-scale change. It is not evidence for environmental dependence by itself.

The environmental test should compare independently anchored distances and reconstructed matter along the paths, with shared calibration and peculiar-velocity covariance. Sky coordinates alone are not a void-density measurement. Source and observer screening, field evolution, refraction and path length must be predicted from the same field model. Success means improved frozen prediction on new objects; failure includes improvement that vanishes after calibration or velocity corrections.

A homogeneous field or perfectly symmetric spherical geometry selects no preferred direction. An arbitrary aligned quadrupole and octopole can be fitted to a target, but that construction has no explanatory force without a source law. An antipodally even response cannot generate an odd octopole. A preferred field direction must generate testable polarization or directional redshift consequences.

Existing Planck topology work finds no compact topology in the tested classes below the relevant last-scattering scale and illustrates the importance of polarization checks.[9] Its quoted lengths are tied to its geometry and radiation-origin assumptions and cannot be imported as bounds on a no-Big-Bang cavity. The useful transferable requirement is to calculate repeated patterns, angular correlations and polarization in the candidate's own geometry. The current project supplies neither a radius nor a unique directional transfer function, so no new sky-map likelihood is claimed.

## H10: preserve gravitational clock behavior before extending galaxy dynamics

For the local phenomenological form Δν/ν=(1+εG)ΔΦ/c², the Galileo satellite analysis gives εG=(0.19±2.48)×10⁻⁵.[4] A Gaussian 95% interval calculated from that summary is approximately −4.67×10⁻⁵ to +5.05×10⁻⁵. This is a constraint on departures in that model and environment, not a fractional bound on every possible time field or every gravitational potential.

Atomic comparisons additionally constrain changes in dimensionless constants. If α_dot/α=qγ is postulated, the previously used atomic-clock result gives approximately |q|≲4.1×10⁻⁸ at the fitted γ.[10] This conditional limit cannot be applied to a common rescaling of every clock, and the illustrative optical constitutive relations can cancel in α. The theory must calculate the actual transition ratios.

Development goal: reproduce gravitational and kinematic clock comparisons, then derive both the force on matter and light bending from the same model. Adding an unrelated galaxy acceleration function would not establish a unified law of time. Rotational, vertical-force and lensing tests should follow only once those observables are determined.

## Prioritized milestones and acceptance criteria

1. **Operational definition.** Specify metric or material rulers, atomic clock frequencies, photon dynamics and the measured redshift. Success is a non-removable observable effect. A pure time-coordinate change ends the branch.
2. **Dynamical source.** Replace the prescribed screened profile with a solution. Success is sustained evolution with known initial/boundary conditions and explicit energy supply. A profile that requires an undisclosed growing external source fails.
3. **Emitter-to-detector calculation.** Derive spectral shift, spectral aging, duration stretch, flux and angular transport together. Success is agreement with the stretching observations and local standards without independent per-observable adjustments.
4. **Thermal transport.** Apply the same field to many radiation paths and model their observed spectrum. Success is an acceptable fit with calibration/foreground covariance and a declared thermal origin or initial condition. Percent-level unthermalized mixing like the tested example is a failure.
5. **Parameter identification.** Measure or constrain screening/environment independently of redshift. Success is breaking the active-fraction/rate and distance-zero-point degeneracies. A fit that measures only their product remains phenomenological.
6. **Discriminating prediction.** Before accessing new labels, freeze a prediction for environment, redshift drift, brightness–angular-distance behavior or polarization. Success is a reproducible improvement on independent data with calibrated uncertainty; matching another model's prediction does not identify a cause.
7. **Unified gravity.** Derive orbital, vertical and lensing observables only after the common action is fixed. Success requires transfer across systems with shared parameters, not a new normalization for each dataset.

The immediate research priority is a physically driven screened field, not additional unconstrained curve fitting. The spectral-aging and thermal-mixing results are suitable additions to the scientific evidence base. The passive-pinning obstruction is equally important and should be reported as a limitation. An honest theory matures by narrowing its permitted mechanisms as well as by matching observations.

## Sources

[1] White, R. M. T., et al. The Dark Energy Survey Supernova Program: Slow supernovae show cosmological time dilation out to z~1 (2024). https://arxiv.org/abs/2406.05050 . Derived data: https://github.com/ryanwhite1/DES-Time-Dilation . Prior local reanalysis archived in time_theory_revision_data.zip.

[2] Blondin, S., et al. Time Dilation in Type Ia Supernova Spectra at High Redshift, ApJ 682, 724 (2008), Table 3 and Sections 3–4. https://arxiv.org/abs/0804.3595 . The present report recalculates fits from all 35 published table rows.

[3] NASA LAMBDA. FIRAS CMB Monopole Spectrum, based on Fixsen et al. (1996) and calibration updates. https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html . Data: https://lambda.gsfc.nasa.gov/data/cobe/firas/monopole_spec/firas_monopole_spec_v1.txt . The present report fits the released residual column.

[4] Delva, P., et al. Gravitational Redshift Test Using Eccentric Galileo Satellites, Physical Review Letters 121, 231101 (2018). https://doi.org/10.1103/PhysRevLett.121.231101 . Published summary constraint; no new satellite-data analysis.

[5] Riechers, D. A., et al. Microwave Background Temperature at a Redshift of 6.34 from H2O Absorption (2022). https://arxiv.org/abs/2202.00693 . Published temperature interval only.

[6] Fixsen, D. J. The Temperature of the Cosmic Microwave Background (2009). https://arxiv.org/abs/0911.1955 . Present temperature normalization.

[7] Tong, D. General Relativity, chapter 1, clock comparisons and geodesic dynamics. https://www.davidtong.org/teaching/general-relativity/grhtml/S1 . Coordinate and endpoint arguments applied here to specified restricted ansatzes.

[8] Khoury, J., and Weltman, A. Chameleon Fields: Awaiting Surprises for Tests of Gravity in Space (2004). https://arxiv.org/abs/astro-ph/0309300 . Khoury, J. Chameleon Field Theories (2013). https://arxiv.org/abs/1306.4326 . Screening precedent and limitations, not validation of the present model.

[9] Planck Collaboration. Planck 2015 results. XVIII. Background geometry and topology. https://arxiv.org/abs/1502.01593 . Model-dependent topology and polarization tests; no bounds imported into the fictional geometry.

[10] Lange, R., et al. Improved Limits for Violations of Local Position Invariance from Atomic Clock Comparisons (2021). https://doi.org/10.1103/PhysRevLett.126.011102 . Published dimensionless-constant constraint.

## Reproduction

The accompanying archive includes this report, run_tests.py, the published spectral-aging PDF and extracted text, the parsed table, the NASA residual table and results.json. Run with Python, NumPy and SciPy. Inputs are hashed in results.json. The table comparisons use published errors; the field examples are explicitly synthetic. The previous DES and 164-group calculations are documented in the separately saved earlier archives. No new blind holdout, microscopic field completion, CMB origin, environment reconstruction or sky-map fit is claimed.
