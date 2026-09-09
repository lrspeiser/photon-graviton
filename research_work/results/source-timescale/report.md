# Source size changes both companion delivery and photon loss

## What this answers

The previous joint source/barrier calculation mostly reflected its generated companion waves. Was that inevitable, or did the assumed conversion-region size make the waves oscillate too slowly? This comparison changes that size in the same Hamiltonian and measures the absolute energy delivered, not just the transmission percentage.

**A narrower tested conversion region delivers more energy through the unchanged barrier, despite transferring less energy out of the photon.** Reducing its half-width from 2 to 0.5 increases delivered energy by about 47 times for restoring frequency m=1 and 932 times for m=3. This demonstrates a source-timescale tradeoff within the candidate. It does not identify the actual size or structure of low-gravity regions, preserve a fitted redshift automatically, or supply permanent capture.

## Fixed assumptions and changed parameter

The [previous solver](../generated-wave-access/check.py) is reused without modification. Its Hamiltonian contains photon energy P/n_bar, a positive companion-field kinetic and gradient energy, and the same lossless restoring slab. The receiving companion sector is unchanged. Initially all field energy is zero; it is supplied by the photon interaction. Static spatial profiles exchange momentum but supply no energy in this approximation.

The conversion envelope is

```
g(x) = cos^2[pi (x-2)/(2 a)] for |x-2|<a, and zero otherwise.
```

Only half-width a changes among 0.5, 1 and 2. Initial photon energy stays 0.03, smoothing width stays 0.25, field stiffness stays K=1, and companion speed stays v=0.5 in c=1 units. The slab remains -6<x<-4, with m=0, 1 or 3. Source and slab remain separated, so the photon does not directly inject energy inside the monitored slab region. All results use the same end time 50 and the same boundary and energy checks as the original calculation.

These are dimensionless model lengths. A finite envelope is an optional diagnostic device, not a derived distribution of voids. In particular, we have not shown why the time/companion interaction would occur in patches of the favorable size. The smoothing width is a coarse-grained coupling parameter, not an identified photon wavelength or particle radius.

## Results

The table uses 6,144 cells. Delivered energy means energy that has crossed the far monitor by time 50; it does not mean captured energy. Percent transmission is relative to companion energy incident on the barrier, while photon loss is relative to the initial photon energy. The same photon loss applies to the three barrier choices within the numerical accuracy of these causally separated runs.

| Source half-width a | Photon energy lost | Transmission, m=1 | Delivered energy, m=1 | Transmission, m=3 | Delivered energy, m=3 |
|---:|---:|---:|---:|---:|---:|
| 2 | 2.8267% | 0.17129% | 3.74094e-7 | 0.000061971% | 1.35342e-10 |
| 1 | 1.2995% | 7.3476% | 8.00164e-6 | 0.0020611% | 2.24455e-9 |
| 0.5 | 0.5116% | 34.2493% | 1.76320e-5 | 0.24515% | 1.26204e-7 |

The smallest source produces only about 18.1% as much total photon energy loss as the widest one. Its leftward incident companion energy is about 23.6% as large. Nevertheless, the much higher transmitted fraction increases absolute delivery. For m=1, the smallest source delivers about 0.0588% of the initial photon energy across the barrier by time 50. This is still a small amount, and no gravitational response is calculated here.

The m=0 controls transmit essentially all of the leftward pulse and therefore reveal the opposite part of the tradeoff: without a barrier, narrowing the source reduces delivered energy. The advantage of narrowing arises from matching the generated frequencies to the receiver's propagation band; it is not extra energy creation.

The generated spectrum also changes. Its median angular frequency rises from approximately 0.209 at a=2 to 0.410 at a=1 and 0.773 at a=0.5. The incident pulse is produced by the common photon/field equations, not manually replaced with higher-frequency waves.

For a=0.5, m=3, about 0.0380% of incident energy remains temporarily between the monitors at time 50. The spectrum-weighted continuum calculation predicts eventual transmission about 0.2658%, compared with 0.2451% already transmitted. Some remaining ringing energy can leave on either side. Reported delivery factors are finite-time comparisons; they are not exact asymptotic efficiencies. The output retains the source-spectrum predictions, temporary energy and both grid resolutions.

## Why a smaller source changes the spectrum and energy

An analytic limit makes the tradeoff intelligible. Neglect photon feedback and smoothing, take a point source moving at c=1 through g, and retain v<1. The retarded one-dimensional wave equation gives a leftward wave velocity amplitude proportional to g evaluated at the emission time:

```
u_left = P g(t_emission) / [2 K (1+v)],
dt_arrival/dt_emission = (1+v)/v.
```

The leftward energy flux is K v u_left^2, so integrating the pulse gives

```
E_left = P^2 / [4 K (1+v)] integral g(t)^2 dt
       = 3 P^2 a / [16 K (1+v)].
```

The corresponding rightward energy is 3 P^2 a/[16 K (1-v)] for this super-companion-speed source. The c=1 source speed and the v<1 condition matter; these expressions must not be extended through v=1. The code independently integrates g squared and checks its value 3a/4. These point-source, weak-loading formulas are limiting derivations, not exact predictions for the finite-smoothing, reacting numerical packet.

The characteristic pulse duration scales with a, and therefore its characteristic frequency scales approximately as 1/a until smoothing or other dynamics intervene. The available energy scales roughly with a in this limit. Narrowing the region thus trades generated energy for a broader, higher-frequency spectrum. Transmission can more than compensate for the lost production in some ranges, as the numerical examples show, but three widths do not establish an optimum or arbitrary improvement as a tends to zero.

There is also a scale requirement. If the slab half-width is b and its static suppression parameter is B=mb/v, the characteristic source frequency divided by the barrier cutoff scales as

```
omega_source/m proportional to [b/a] / [(1+v) B]  (c=1).
```

The proportionality depends on how pulse frequency is defined and on smoothing. It shows why a very extended slowly varying source and a strongly suppressing small region can be poorly matched. The favorable examples do not establish this scale relation for actual voids and galaxies.

## Numerical verification and retained failure

The initial 1,536/3,072-cell pilot failed the predeclared 2% relative-energy convergence gate: the a=0.5, m=3 transmitted-energy comparison differed by about 3.76%. That failure is retained in [pilot-verification.json](pilot-verification.json). The gate was not relaxed.

The final 18 runs use 3,072 and 6,144 cells for all nine source/barrier combinations. The largest relative change among photon loss, incident energy and transmitted energy is about 0.950%, meeting the unchanged 2% criterion. Both resolutions retain the full energy and temporary-storage ledgers, numerical momentum-reaction check and boundary-contamination check. No observational success is inferred from these numerical passes.

## Consequences and next requirement

The previous poor-transmission result is not a general prohibition on energy delivery. It applies to its generated source timescale. The present comparison finds a route to better delivery within the same candidate, but also changes the amount of redshift-producing photon loss. It would be incorrect to retain the old redshift normalization while silently replacing the production region with a smaller one.

A next calculation must derive or explicitly model the distribution and repeated action of conversion regions, including their evolving fields, and predict accumulated frequency shift, whole-event timing and fluctuations together. Simply multiplying independent single-region losses can miss interference, backreaction and source history. A favorable region size must ultimately follow from a physical environmental law and remain compatible with matter clocks. Capture, stable spatial support and a common motion/lensing response remain unresolved even after successful entry.

The full 20-task, 32-area goal remains active. The current academic manuscript is a dated evidence snapshot; this result and the preceding generated-wave access result are inputs to its next revision.

## Reproduction

```sh
python research_work/results/source-timescale/check.py
```

Requires NumPy and SciPy and imports the existing generated-wave solver. The [protocol](protocol.json) records the scan and unchanged convergence gate. The [saved output](source-timescale-results.json) includes all runs, controls, spectra, source hashes and comparisons. The standard diagnostic suite includes this calculation.
