# Spherical spreading: startup redshift approaches a fixed delay

The same optional photon/companion interaction behaves differently when waves spread spherically in three dimensions. A steady localized light supply builds a field that approaches a static profile. New signals then retain an extra travel delay but acquire very little additional redshift. The one-dimensional sustained rolling result cannot be assumed to describe a localized three-dimensional source.

This is a limitation of the tested scalar candidate and source geometry. It does not reject the user's cumulative-time premise, exclude a distributed cosmic background, or test nonlinear graviton self-binding. Earlier redshift is not erased when a photon leaves a region or when the field later settles.

## Model and energy supply

The [derivation](derivation.md) obtains the spherical equations from the energy-weighted radial Hamiltonian. The field is w=r phi, with regularity at the origin, and every radiation shell carries a stated total energy over all directions. The radial area factors are retained in both production and field energy. These are spherical symmetry calculations in three spatial dimensions, not general three-dimensional structure simulations.

The [protocol](protocol.json) fixes K=1, field speed v=0.5, Gaussian smoothing width 0.25 and cosine-squared production half-width 0.5 in reference c=1 units. The compact cases have centers at radii 3, 6 and 9, a source at 0.5, and a detector at 12. The extended case adds centers 12, 15 and 18 and uses a detector at 21. No physical distance or stellar luminosity calibration is assigned.

The compact source emits for 64 time units; its finite initial fuel energy is 19.2 at power 0.3 and 38.4 at power 0.6. The extended source emits for 112 units, with initial energy 33.6. Energy is transferred from fuel into each newly emitted shell with its current optical factor included. Nothing is supplied after the finite fuel budget is exhausted. Emitter microphysics, atomic clock response and material recoil are not derived by this prescription.

Two grids, with radial spacings 0.0625 and 0.03125, are run for each of four cases. The compact domain ends at radius 56 and time 80; the extended one ends at radius 100 and time 140. The exterior boundary is beyond significant wave propagation for the recorded intervals. The origin implements regular spherical propagation. Photon shells may leave the field domain after leaving all production support; their energy remains in the total ledger.

The compact cases compare packet intervals 0.5 and 0.25 at fixed power 0.3, then double power at interval 0.25. The extended case uses power 0.3 and interval 0.25. Nine or ten nonuniform launch times give 74 independent frequency/event comparisons across eight backgrounds. Later launch times are 44.31, 52.07 and 60.19 in the compact geometry, and 92.07, 100.19 and 108.13 in the extended geometry. Their passage through the production regions precedes the local end of illumination.

## Result

Fine-grid measurements:

| Case | Largest sampled startup/transient z | Largest absolute later z | Last sampled z | Later mean excess travel time |
|---|---:|---:|---:|---:|
| Compact, power 0.3, interval 0.5 | 0.00170592 | 1.48047e-6 | 2.65972e-7 | 0.0328210 |
| Compact, power 0.3, interval 0.25 | 0.00168519 | 1.50460e-6 | 2.66470e-7 | 0.0328209 |
| Compact, power 0.6, interval 0.25 | 0.00337132 | 5.99375e-6 | 1.06746e-6 | 0.0656284 |
| Extended, power 0.3, interval 0.25 | 0.00205889 | 1.25284e-6 | 4.88374e-7 | 0.0756272 |

The table reports finite samples, not a proof of the infinite-time limit or an exhaustive maximum over every launch phase. It shows a large reduction of the redshift increment while travel delay remains. Doubling source power increases the delay and transient redshift, but does not preserve the earlier rolling behavior. Increasing the radial extent lengthens the transient and increases delay; it does not establish a permanently rolling localized field.

An independently derived static radial profile predicts the field without fitting its amplitude. In the late field windows, the maximum profile difference is below 0.073% of the predicted maximum amplitude for all fine-grid cases. Small field velocities and residual redshifts remain, so the finite results are described as approaching the static profile, not as an exact stationary solution.

All raw probes, finite event intervals, static profiles, emitter-energy transfers and sampled energy ledgers are in [spherical-propagation-results.json](spherical-propagation-results.json). The final field energy includes turn-on/turn-off radiation and field energy still present at the chosen final time. It must not be relabelled as permanently deposited gravitational mass.

## Verification

The total ledger is emitter fuel plus photon energy plus companion-field energy. Every photon-energy loss is paired with field work, including during the source's response to earlier radiation. Its maximum relative residual across the eight runs is below 5.1e-14. This near-roundoff residual checks the implemented Hamiltonian; it does not establish fundamental physics.

Grid refinement changes final field energy by less than 0.00446%, an individual probe stretch by less than 6.94e-6 absolutely, and predicted static delay by less than 2.69e-6 in reference-time units. Independent frequency and event-time integrations use the same declared standards and pass the protocol's 1e-7 relative agreement and 1e-7 arrival-difference gates. A separate 1e-7 absolute refinement gate covers the much smaller late signal; its largest measured grid difference is below 6.97e-9. Across all runs, frequency/event disagreement is below 9.23e-11 relatively and independent arrival difference below 2.99e-13. The saved output retains the actual residuals. Outer boundary energy stays below the frozen threshold.

For efficiency, emitted shells more than twelve Gaussian widths outside the production support are omitted from the source kernel. Their neglected Gaussian factors are below exp(-72); the probe uses the full kernel. The total energy and grid checks cover the resulting calculation. This is a numerical tail cutoff, not an extra physical capture or attenuation law.

Whether the redshift becomes small is a measured result, not a requirement for passing the numerical checks. No unfavorable case has been discarded to obtain a static outcome. Angular structure, gravitational backreaction, interactions between companions and cosmic illumination remain absent.

## Consequence for the research plan

The current free scalar field driven by a steady compact radiation source is not a demonstrated sustained cosmic redshift mechanism. It can produce transient cumulative stretching, but a static field mainly supplies a common delay. The required low-gravity time law must therefore explain what maintains its evolution in three dimensions.

The user's [self-binding and time-formula question](../../../research_plan/time-field-and-self-binding.md) is recorded as a possible extension. Self-binding could address retention, but a static bound store alone would not maintain partial_t n. A next calculation must derive the interaction and its energy, capture and stability consequences and determine whether it supports continuing phase/event stretching. A distributed finite source population is a separate unfinished route. All 20 tasks and 32 requirements remain in scope.

## Reproduction

```sh
python research_work/results/spherical-propagation/check.py
```

NumPy and SciPy are required. The standard runner includes this calculation and writes fresh generated results while preserving saved evidence. No expansion, dark-matter component or Big-Bang source history is imported.
