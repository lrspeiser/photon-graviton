# Longer illumination approaches a derived rolling branch

## Result in plain language

In this optional one-dimensional model, the time/companion field can keep growing while successive faint signals receive almost the same stretch. Longer, smoother illumination brings both the redshift and the event-duration stretch close to a value derived from the same interaction and source power. The earlier large drift is substantially reduced in the tested later windows.

This is a conditional model result, not an astronomical fit or a proof that the state is stable forever. In particular, waves spreading in three dimensions may behave differently. Low-gravity coupling, matter clocks, actual companion capture and extra gravity remain unproved.

## What was derived

The [derivation](derivation.md) starts from the previously tested Hamiltonian, including the radiation's effect on the field. For a rolling field phi=u t+B(x), constant incoming photon power Q, field inertia K, companion speed v, and total smoothed coupling G, it predicts

```
4 K v u^2 = Q [1-exp(-2uG)],   u>0,
frequency stretch = local event-duration stretch = exp(uG),
field-energy growth = in-flight photon-energy growth = 2Kvu^2.
```

The positive rolling rate is solved from this equation; it is not fitted to the simulation's redshift. In these tests K=1, v=0.5 and G=1.5. All units are model reference units, not calibrated cosmic times, distances or luminosities.

The last equality has an important bookkeeping consequence. The difference between incoming and outgoing photon power includes both energy transferred to companions and energy remaining in photons whose arrival is delayed. All energy actually lost by an individual photon still goes to the field. Counting the entire detector brightness deficit as new companion energy would overcount the supply in this branch.

## Finite histories and controls

The final [protocol](protocol.json) places all photons in the initial state; total initial radiation energy is Q times 48, equal to 0.144 or 0.288. The field starts with zero energy. Packet intervals 0.5 and 0.25 at Q=0.003 test radiation sampling at fixed power; Q=0.006 and interval 0.25 test source-power dependence. There are 96 or 192 packets. No infinite energy source or background reset is supplied.

Each case runs on 3,072 and 6,144 cells over -56<x<40 to time 60. The three conversion profiles, widths and Gaussian smoothing are unchanged. The initial source is at x=-2 and the detector at x=9. Ten nonuniform probe launches between 6 and 42.07 give 60 independent carrier/event comparisons across the six backgrounds. Six launches from 30.19 to 42.07 form the later summary. These are deterministic samples, not a random source population or an observational uncertainty distribution.

Energy rates are measured in the fixed time window 40 to 46, before the finite train ends. A smooth control window between x=-2 and x=9, with Gaussian edge width 0.25, records in-flight photon energy, incoming and outgoing powers, and conversion work. Its exact identity is

```
d E_photon,window/dt = incoming - outgoing - conversion,window.
```

The conversion term uses the same window weights, so small tails are accounted for rather than silently dropped. Both grids use identical time-window endpoints. Cubic interpolation integrates the sampled fluxes; Hermite interpolation of energy uses its independently specified derivative. This avoids accidentally changing the measurement window with the solver's adaptive time grid.

The [shorter-history protocol](shorter-history-protocol.json), [results](shorter-history-results.json), and [exact earlier source text](shorter-history-check.py.txt) preserve the preceding duration-32 comparison. That run approached the predicted redshift and field-energy growth, but its time-24-to-30 in-flight growth rate was still about 17-18% below the rolling prediction in the interval-0.25 cases. It motivated the longer history. Its windows use accepted solver nodes and report the actual endpoints; it is not overwritten with the final interpolation method. The archived script and protocol hashes correspond to their historical names inside the archived result. They are evidence snapshots, not additional default runner jobs.

## Measured redshift and timing

Final fine-grid results:

| Incoming power | Packet interval | Predicted z | Later mean z | Mean difference from prediction | Later linear dz/dt |
|---:|---:|---:|---:|---:|---:|
| 0.003 | 0.5 | 0.00672747 | 0.00672804 | +0.00840% | -1.13e-7 |
| 0.003 | 0.25 | 0.00672747 | 0.00672804 | +0.00852% | -1.13e-7 |
| 0.006 | 0.25 | 0.01341087 | 0.01341302 | +0.01601% | -4.33e-7 |

These percentages are relative to z, not the full duration factor 1+z. The later sampled z ranges are approximately 0.00672746 to 0.00672854 at Q=0.003 and 0.01341078 to 0.01341491 at Q=0.006. The remaining small negative slopes are retained, not set to zero or removed by detrending. Earlier probes and all finite event intervals remain in the [saved results](sustained-illumination-results.json).

This evidence supports approach toward the rolling branch in these cases; it does not prove asymptotic convergence for all initial states or source histories. The background still grows, and signals emitted later still spend longer traveling. Nearly steady stretching does not mean n is static or that a source of finite fuel can maintain the state indefinitely.

## Energy result

At Q=0.003 the rolling prediction for each of field-energy growth and in-flight photon-energy growth is about 1.99806e-5. At Q=0.006 it is about 7.88750e-5.

| Incoming power | Packet interval | Field growth versus prediction | In-flight photon growth versus prediction |
|---:|---:|---:|---:|
| 0.003 | 0.5 | -0.000515% | -0.92375% |
| 0.003 | 0.25 | -0.000314% | +0.03148% |
| 0.006 | 0.25 | -0.001273% | +0.06046% |

Finer packet spacing improves the finite-window photon-inventory comparison. The more widely spaced train retains discrete arrival effects even after the mean stretch becomes close to the continuum prediction. It would be incorrect to claim every finite packet statistic equals the continuum branch.

For the smoother Q=0.003 case, the measured incoming power is approximately 0.003, outgoing power 0.00296003265, conversion work 1.99805e-5 and photon-inventory growth 1.99869e-5. These close the control ledger. No deposited reservoir or new gravitational mass is inferred from the result.

## Numerical verification and limits

- The compact field-history interpolation reproduces a saved earlier background's final energy exactly and four reference probe stretches within 8.6e-11 absolutely. The longer solver uses the same Hamiltonian and numerical method, retaining only interaction-region field histories for probe interpolation.
- All 60 independent frequency/event comparisons agree within 9.9e-12 relatively, and independent arrival times within 1.2e-13 in reference units. These use declared common standards, not derived atom clocks.
- Across all runs, relative total-energy drift is below 1.6e-15 and the momentum-rate ledger residual below 5.6e-16. The latter includes the fixed profile/grid reaction; it is not a complete dynamical matter momentum theory.
- Grid refinement changes final photon loss by less than 0.000311%, later mean redshift by less than 0.0000033%, and measured energy-growth rates by less than 0.000175%. The largest individual probe-stretch change, including early probes, is below 9.42e-6 absolutely.
- Integrated control-volume energy residual is below 2.6e-7 of conversion power. Boundary-strip energy checks exclude relevant periodic wave return during the test.

The numerical gates test the calculation and its accounting. Agreement with the rolling prediction is measured, not imposed as a gate that discards disagreeing cases. No real observation has been fitted or passed by these checks.

## Next requirement

Test dimensionality before extrapolating this branch to a universe. Derive how a companion field driven by sustained radiation behaves when it can spread through three-dimensional space, and compare localized sources with extended illumination under finite causal energy budgets. Determine whether persistent time stretching survives or whether it depends on the one-dimensional propagation geometry. Keep the same energy recipient, include radiation backreaction, and predict frequency and event timing together.

Source evolution, clock response and an environmental dependence on weak gravity must still be derived. Entry into wells, irreversible capture, support against escape and a shared lensing/stellar-motion response remain separate unfinished requirements. The complete 20-task, 32-area goal remains active.

## Reproduction

```sh
python research_work/results/sustained-illumination/check.py
```

The diagnostic uses NumPy and SciPy and writes to the ignored generated-results directory. The standard runner also includes it. Source hashes, both grids, raw probes, finite event intervals, numerical checks and sampled energy histories are retained in the result JSON. The paper's version-0.1 manuscript remains a dated evidence snapshot; its README points to this later result.
