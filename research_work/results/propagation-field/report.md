# Evolving propagation field (PF-1)

**Result.** The proposed wave law is internally consistent. One driven simulation reproduces the redshift, event stretching, photon-number conservation and energy exchange together. PF-1 is not promoted, for two reasons:
- **Light cannot distinguish it from expansion.** Its radiation sector is exactly that of an expanding universe with scale factor a = n; no light-only observation separates the two.
- **Its brightness law is disfavored.** It is the flux law already scored on Pantheon+, and it again loses to a matched FLRW comparator.

Whether the universe is "nonexpanding" is therefore decided entirely by how clocks and rulers couple to n. Every completion of that coupling tested in this repository fails.

**Files.**
- Protocol: [protocol.md](protocol.md), declared before execution at `main` 636a906.
- Scripts: [wave.py](wave.py) (library), [pf1.py](pf1.py) (T1–T3) and [brightness.py](brightness.py) (T4).
- Results: [pf1-results.json](pf1-results.json) and [brightness-results.json](brightness-results.json).

## Provenance

- **Proposed (project hypothesis).** P1–P3:
  - a homogeneous index n(t), measured against fixed material clocks and rulers;
  - one wave law shared by electromagnetic and gravitational waves;
  - a dynamical n with its own inertia.
- **Already in the repository.**
  - [matched-wave/derivation.md](../matched-wave/derivation.md) has the same Lagrangian, its exact optical-time solution, the exchange law and the charged-matter blueshift.
  - Clock and ruler completions are classified in [clock-ruler-completion](../clock-ruler-completion/report.md), [matter-clock-closure](../matter-clock-closure/report.md) and the [shared interaction test](../../../shared_interaction_test/shared_interaction_report.md).
  - The flux law was scored in [brightness.py](../../../shared_interaction_test/brightness.py) and in the [rate report](../brightness-distance-consistency/rate-report.md).
- **Established physics used.** Wave propagation in a time-varying, impedance-matched medium; the conformal invariance of Maxwell's equations; generalized least squares.
- **New in this run.**
  - T1, a driven simulation in physical (material-clock) time.
  - T2's restoring-potential case, which shows exchange in both directions.
  - The identification of the radiation sector with coasting FLRW as an observational degeneracy.
  - T4's comparator at equal freedom.

## Declared tests

**T1. One driven simulation, all observables at once (pass).** An emitter running on material-clock time drives two Gaussian pulses into the wave law: carrier ν = 1, σ = 3, emitted at t = −100 and t = −76.
- **Integration.** Physical time, a staggered leapfrog with a fourth-order Laplacian, and 32 points per emitted wavelength.
- **Index.** Linear, giving z = 1 at the fixed separation D = ln 2/ṅ with ṅ = 0.005.
- **Grid.** Long enough that nothing reaches the outflow edge before the run ends.

| Quantity | Predicted (control) | Linear index | Constant-index control |
|---|---:|---:|---:|
| Carrier ratio, first pulse | 2 (1) | 1.99979 | 0.9999984 |
| Carrier ratio, second pulse | 2 (1) | 1.99982 | 0.9999984 |
| Envelope-width ratio | 2 (1) | 1.99982 | 0.9999984 |
| Pulse-spacing ratio | 2 (1) | 1.99989 | 1.0000000 |
| Largest drift of wave energy × n after emission | 0 | 7.6×10⁻⁵ | 3×10⁻¹³ |

The declared tolerances were 0.5% for the linear index and 0.2% for the control. The wave energy halves (168.4 to 84.3 in code units) while n doubles (0.695 to 1.389). Energy × n stays fixed, so photon number is conserved.

**T2. Coupled index and radiation (pass).** Four modes (k = 1, 2, 3, 5) were coupled to a dynamical n with inertia M = 40. The last four columns are errors.

| Case | Index range | Wave energy | Total energy | Photon number per mode | Exchange law dE_w/dt = −(ṅ/n)E_w | Reduced mechanics |
|---|---|---|---:|---:|---:|---:|
| Free, V = 0 | 1 → 28.3 | 0.201 → 0.0071 | 1.6×10⁻¹² | 1.0×10⁻¹¹ | 8.6×10⁻⁷ | 6.7×10⁻¹² |
| Restoring, V = 0.01(n−1.5)² | 1.0–5.6, oscillating | 0.036–0.201 | 2.3×10⁻¹² | 2.0×10⁻¹¹ | 9.1×10⁻⁷ | 1.4×10⁻¹¹ |

- The exchange residual is limited by the finite-difference derivative, not by the physics.
- Radiation pushes the index, M n̈ = −V′(n) + E_τ/n², so a free index is driven forward and then coasts.
- With a restoring potential, energy returns to the waves (a blueshift) whenever n falls.
- The free case repeats the matched-wave checks; the restoring case is new.

**T3. Closed-form observables (reported).**
- 1+z = n_o/n_e; a linear index gives 1+z = e^{αD} with α = ṅ/c.
- The time-dilation exponent is b = 1 exactly.
- D_L = (1+z) ln(1+z)/α.

Against flat FLRW with Ω_m = 0.3 and the same low-redshift slope (H0 = cα):

| z | 0.1 | 0.5 | 1 | 1.5 |
|---|---:|---:|---:|---:|
| μ_PF1 − μ_FLRW (mag) | −0.054 | −0.182 | −0.232 | −0.231 |

**T4. Exposed Pantheon+ brightness (reported; not blind).**
- **Data.** 960 light curves at z ≥ 0.1 with the full released covariance, and the absolute magnitude from 77 Cepheid-host light curves.
- **Freedom.** Each model has one scale parameter. It enters as a redshift-independent magnitude shift, so the best fit equals the free-offset fit.
- **Bin contrasts.** Mean residual of each bin minus the 0.1–0.3 bin.

| Model | χ² (960 rows) | Scale | 0.3–0.6 minus first | 0.6–1 minus first | 1–3 minus first |
|---|---:|---:|---:|---:|---:|
| PF-1 | 871.55 | cα = 69.76 ± 0.99 km/s/Mpc | +0.062 ± 0.013 | +0.080 ± 0.024 | +0.225 ± 0.065 |
| FLRW, Ω_m = 0.3 | 836.51 | H0 = 74.10 ± 1.06 km/s/Mpc | +0.005 ± 0.013 | −0.036 ± 0.024 | +0.091 ± 0.065 |

- **Δχ² = +35.0 at equal freedom.** Both χ² values reproduce the archived shape-only scores (871.552 and 836.512).
- **Distant supernovae are fainter than PF-1 predicts.** No value of α changes the bin contrasts.
- **The rates disagree, as the archive already found.** The fitted cα differs from the galaxy-group value of 74.62.
- **Time dilation is consistent.** The declared comparison is the archived DES width test ([time_revision/analyze.py](../../../time_revision/analyze.py)). It gives b = 1.0029 ± 0.0048, consistent with b = 1, but its reference light curves are already dilation-corrected.

## Why light alone cannot decide

With τ = ∫dt/n, the wave law becomes the flat wave equation in (τ, x). Maxwell's equations are conformally invariant, so the medium ε = μ = n(t) is equivalent to the metric ds² = −dt² + n(t)² dx². That is FLRW with scale factor a = n, conformal time τ and comoving coordinate x. It is the r = 2 member of the [clock-ruler family](../clock-ruler-completion/report.md) (N = 1, A = n). A linear index is the coasting universe a ∝ t, which is why D_L = (c/H0)(1+z) ln(1+z) is the known coasting-cosmology distance law.

Every observable carried by light or gravitational waves alone is therefore identical in PF-1 and in that expanding universe:
- redshift and event stretch;
- photon number and a blackbody temperature falling as 1/n;
- flux;
- the arrival delay between electromagnetic and gravitational waves.

Only the matter standards separate the two.

- **Rulers.** Rulers fixed in x give D_A = D, so D_L = (1+z)D_A. The distance-duality ratio η = D_L/[(1+z)²D_A] is then 1/(1+z), where expansion gives 1. Tolman surface brightness falls as (1+z)^−2 instead of (1+z)^−4. The published early-type exponents to z ≈ 0.9 are 2.6 (R band) and 3.4 (I band) (Lubin & Sandage 2001). Read under PF-1, they would require galaxies that were fainter in the past, the opposite of passive stellar evolution. That reading is model-inferred, because the test is degenerate with luminosity evolution.
- **Clocks and atoms.** This decides the branch.

| Completion (where tested) | Measured redshift | Problem |
|---|---|---|
| Ordinary charges and masses in the ε = μ = n field ([matched-wave](../matched-wave/derivation.md) §4, [atomic-line-response](../atomic-line-response/report.md)) | 1+z = n_e/n_o, a blueshift | Atomic lines scale as n^−2, and different transitions shift differently |
| Universal clock, dτ = dt/n ([matter-clock-closure](../matter-clock-closure/report.md)) | None (below 3×10^−16) | Signal and clock change together |
| Partial clock response, dτ = dt/n^β (same) | (n_o/n_e)^(1−β) | Local light speed changes by 0.76% at 100 Mly |
| Universal metric for clocks and rulers ([clock-ruler-completion](../clock-ruler-completion/report.md)) | Full shift at r = 2 | Physical separations grow as A = n; this is expansion |
| Fixed electrostatic binding ([shared interaction test](../../../shared_interaction_test/shared_interaction_report.md)) | Full shift in gross structure | Magnetic hyperfine scales as n², so optical/Cs ratios drift at about −2ṅ/n ≈ −1.5×10^−10 per year. Lange et al. (2021) measure the Yb⁺ E3/Cs drift as −3.1(3.4)×10^−17 per year |
| Spatially varying n(x,t) ([inhomogeneous-clock-completion](../inhomogeneous-clock-completion/)) | Kinematic only | Needs a holding force c²∇ln(1/n); no dynamics |

P1 leaves matter uncoupled from n. That is the "nonuniversal light interaction with unchanged local matter standards" that clock-ruler-completion left open. The first and fifth rows show why matter cannot simply be left uncoupled: atoms are held together by the same electromagnetic field whose propagation n changes, since static and radiative fields are one field. A completion must name a matter action in which gross, fine and hyperfine atomic structure are all independent of n while free radiation sees ε = μ = n. None has been found.

## Common electromagnetic and gravitational propagation

Both kinds of wave obey the same law, so they share the τ-mapping.
- **Arrival times.** For any n(t) they arrive together, which satisfies GW170817 by construction. The [joint propagation audit](../joint-propagation-audit/report.md) shows that photon-only transport would leave light at least 674,646 years behind.
- **Chirp.** A binary's chirp is stretched by 1+z, so a detector infers a chirp mass (1+z) times the source value. Its amplitude, which falls as 1/D in the flat τ frame, then implies a siren distance D_GW = D(1+z) = D_L.
- **Siren caveat.** This holds only if gravitational-wave generation has no dependence on n. Any dependence would appear as D_GW ≠ D_L, an observable handle on a completion.

## Limits

- **Linear n is an interval model.** n = 1 + ṅt reaches zero 1/ṅ = 13.1 Gyr before the present (cα = 74.62 km/s/Mpc), the counterpart of a coasting universe's age. Light emitted near that time has unbounded redshift. The history of n outside the interval is not modeled.
- **The field channel does not feed galaxies.** Q = (ṅ/n)u_γ feeds a homogeneous field, not galaxy reservoirs. The same supply rate was found about 6×10^5 short of reservoir growth ([optical-growth report](../brightness-distance-consistency/optical-growth-report.md)).
- **No thermal origin.** PF-1 preserves a thermal spectrum (T ∝ 1/n) but does not create one.
- **Homogeneous n only.** Lensing and local gravity are outside this branch.

## Decision

PF-1 is kept as an internally consistent propagation law and is not promoted:
- on exposed data its brightness law is disfavored at equal freedom (Δχ² = +35);
- its claim to nonexpansion needs a matter completion that no tested coupling supplies.

For the other branches it supplies two things:
- **A rate.** c|d ln n/dt| = c²α ≈ 7.25×10^−10 m/s² at the present epoch; RPG-1 tests it as the origin of a*.
- **An exchange channel.** A conserved exchange between radiation and the index.

Reproduce with `python research_work/results/propagation-field/pf1.py` and `python research_work/results/propagation-field/brightness.py`. Both regenerate into a fresh directory and compare their numbers with the archived results. Both are also jobs in `research_work/run_checks.py`.
