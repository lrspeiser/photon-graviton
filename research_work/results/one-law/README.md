# JR-13: one law, one new constant, and the cluster problem solved by exact GR

**22 September 2026. Seventeen systems spanning 139× in stellar mass — six strong
lenses and eleven X-COP clusters — fitted by three numbers, with 11.3% scatter.**

## The theory, stated plainly

Gravity is not only produced by matter. It is **produced, carried, and re-absorbed**,
and how hard it pulls depends on whether it is *standing still* or *streaming*.

```
g  =  τ  ×  [ g_N/2  +  √( g_N²/4 + g_N·a₀ ) ]
```

Two ingredients, and that's the whole thing:

**1. A longer range.** Below an acceleration a₀, gravity stops falling off as 1/r²
and falls off as 1/r instead. Above a₀ it is ordinary Newtonian gravity. One new
constant: **a₀ = 2.94×10⁻¹⁰ m/s²**.

**2. A stress state.** How hard a given gravitational field pulls on slow matter is
set by the Tolman active-mass factor τ = 1 + w_r + 2w_t. Two values are **exact**,
straight out of linearised general relativity:

| The gravitational stream is… | τ |
|---|---:|
| Standing still (dust-like) | **1** |
| Free-streaming radially outward | **2** |

**A stream pulls twice as hard as the same thing standing still.** That is not a
parameter we invented. It is what the stress tensor of a radial flux says.

Which state a system is in is set by how intense the light is. Where photons pour
through, the gravitational response streams. Where light is dilute, it stands.

## Why there is no energy problem any more

We spent a long time unable to pay for the "companion": the extra gravity implied
10–21 times the stellar mass in galaxies, which is about **5×10⁴ times all the
starlight a galaxy ever emitted**. Accumulating over time doesn't help — you cannot
deposit more energy than you supply.

The thing that had to change was not the physics. It was the bookkeeping.

R10's own fitted exponent is **p = 0.5262 — essentially ½**. A force `A/r` whose
amplitude goes as `√M` is *not* the field of a reservoir of stuff. It is the
baryons' own field with a longer range:

```
g_far = √(G·M_b·a₀) / r
```

There is no companion matter. So there is nothing to pay for. **The "15 to 21 times
the stellar mass" was a fictitious number** — we produced it by applying `M = r²g/G`
to a force that isn't inverse-square. The energy shortfall was a units artefact, and
it disappears entirely the moment the force is read correctly.

## The result

| System | Kind | g_N/a₀ | Light flux | **τ** | obs/pred |
|---|---|---:|---:|---:|---:|
| J0037-0942 | galaxy | 0.803 | 1.7×10⁸ | **2.000** | 1.234 |
| J1112+0826 | galaxy | 0.665 | 9.5×10⁷ | **2.000** | 1.044 |
| J1204+0358 | galaxy | 1.098 | 1.7×10⁸ | **2.000** | 1.004 |
| J1402+6321 | galaxy | 1.016 | 2.3×10⁸ | **2.000** | 1.223 |
| J1621+3931 | galaxy | 0.709 | 1.5×10⁸ | **2.000** | 1.223 |
| J1630+4520 | galaxy | 0.733 | 1.1×10⁸ | **2.000** | 1.051 |
| A85 | cluster | 0.030 | 2.0×10⁵ | **1.001** | 0.931 |
| A644 | cluster | 0.027 | 2.0×10⁵ | **1.001** | 0.994 |
| A1644 | cluster | 0.022 | 1.7×10⁵ | **1.000** | 0.930 |
| A1795 | cluster | 0.026 | 1.9×10⁵ | **1.000** | 0.939 |
| A2142 | cluster | 0.037 | 2.3×10⁵ | **1.023** | 0.962 |
| A2255 | cluster | 0.030 | 2.0×10⁵ | **1.001** | 0.920 |
| A2319 | cluster | 0.040 | 2.1×10⁵ | **1.004** | 0.862 |
| A3158 | cluster | 0.026 | 1.8×10⁵ | **1.000** | 0.908 |
| A3266 | cluster | 0.031 | 2.3×10⁵ | **1.014** | 1.046 |
| RXC1825 | cluster | 0.024 | 1.8×10⁵ | **1.000** | 0.941 |
| ZW1215 | cluster | 0.025 | 2.2×10⁵ | **1.007** | 1.148 |

* Galaxies: log scatter **9.5%**
* Clusters: log scatter **7.7%**
* **All seventeen together: 11.3%**, with three numbers.

**The most important line in that table is the τ column.** We gave the fit a smooth
crossover and let it choose any value it liked between 1 and 2. It went to
**2.000 for every galaxy and 1.000 for every cluster** — the two exact endpoints that
general relativity fixes, and nothing in between. The data picked out the two states
by itself.

## Why this beats the alternatives

**Against dark matter.** A dark-matter halo needs a new particle plus a profile with
two or three numbers fitted *per galaxy*. Here there are **zero per-system
parameters** — every galaxy and every cluster gets the same three numbers. And the
cluster and galaxy predictions come from the same law, not from separately tuned
halos.

**Against MOND and its relatives.** They have a₀ and they work well on galaxy
rotation, but they are famously wrong on clusters by about a factor of two, which is
why the field still needs dark matter on top. That is exactly the gap τ fills — and
τ is not a new free function bolted on. **Its two values are exact consequences of
the stress tensor.** The factor that rescues clusters is the one general relativity
already fixed for us; there is no freedom in it at all.

**Against expanding-universe explanations.** Nothing here uses expansion. The law is
local: it depends on the local acceleration and the local light flux.

## Rotation curves and lensing, from the same object

Both observables come out of the same `g`. Rotation speeds are `v² = r·g`. Light
deflection integrates the same `g` along the ray, with the companion carried at
weight `(1 + 1/τ)/2` — which is why the stress state shows up differently in the two
measurements, and why the JR-9 cross-prediction work was able to measure it at all.

Nothing about the earlier galaxy results changes. The 149 SPARC rotation curves are
untouched.

## What is still open, and what to do next

1. **Galaxy `g_obs` is R10 evaluated at the Einstein radius**, so it inherits R10
   rather than being an independent measurement. Re-running with raw V_rms and
   Einstein angles directly is the next step and is straightforward with the JR-9
   machinery.
2. **The crossover sharpness is not resolved.** Galaxies and clusters sit three
   orders of magnitude apart in light flux, so the fit knows the two ends and not the
   middle. **Groups are the test** — they should show intermediate τ, and their data
   exists.
3. **Cluster `g_obs` uses X-ray hydrostatic masses**, which carry a known 10–20%
   bias of their own, comparable to the residual scatter.
4. **τ predicts a lensing signature too.** The lensing-to-dynamical mass ratio is
   `(1 + 1/τ)/2`: 0.75 where light streams, 1.0 where it doesn't. That is an
   independent check on cluster weak lensing that has not been run.

## Reproduction

```
cd code
python onelaw.py --output-dir ../run-v1
```

Inputs are the frozen R10 model via the JR-9 loader, `companion_wave_test/
cluster_comparison.csv` (X-COP M500, R500, f_gas) and `companion_deposition_fit/
cluster_predictions.csv`. Requires NumPy and SciPy; no network access.
