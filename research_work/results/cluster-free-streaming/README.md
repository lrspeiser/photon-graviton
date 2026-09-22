# JR-11: clusters are short by exactly the free-streaming factor, and that number is not fitted

**22 September 2026. Eleven X-COP clusters, the archived galaxy-calibrated transfer,
and one number that comes from general relativity rather than from a fit.**

## The problem in one sentence

Take the companion law calibrated on galaxy rotation curves, apply it to clusters,
and it predicts about **half** the extra mass their X-ray profiles require.

| Cluster | Required / predicted |
|---|---:|
| A2319 | 1.704 |
| A3158 | 1.824 |
| A2255 | 1.849 |
| A85 | 1.879 |
| A1644 | 1.879 |
| A1795 | 1.899 |
| RXC1825 | 1.904 |
| A2142 | 2.001 |
| A644 | 2.032 |
| A3266 | 2.191 |
| ZW1215 | 2.416 |

Median **1.899**, mean **1.962 ± 0.059**. Notice how *tight* that is. It isn't a mess
that needs eleven explanations. It's one number, nearly the same everywhere.

## The formula

How hard a given amount of energy pulls on slow matter is set by its internal
stresses, through the Tolman active-mass factor:

```
τ  =  1 + w_r + 2 w_t          w = stress ÷ energy density
```

Two states are exact, no fitting anywhere:

| What the companion is | w_r | w_t | **τ** |
|---|---:|---:|---:|
| Sitting still (static dust) | 0 | 0 | **1** |
| **Streaming outward (free-streaming radial flux)** | 1 | 0 | **2** |

**A stream pulls twice as hard as the same energy sitting still.** That is ordinary
general relativity, not a new assumption.

Now the key step. Rotation curves are a *dynamical* measurement, so a galaxy-calibrated
amplitude already has τ_galaxy baked into it. Transfer it to clusters and you predict
their dynamical mass short by exactly:

```
shortfall  =  τ_cluster / τ_galaxy
```

If galaxies sit in the dust state and clusters in the streaming state, that is **2**.

## The test

| Model | Parameters fitted | Scatter |
|---|---:|---:|
| Static dust, τ = 1 *(what the transfer assumed)* | 0 | **98.0%** |
| **Free-streaming flux, τ = 2** | **0** | **9.5%** |
| Best fitted constant (1.962) | 1 | 9.5% |
| Crossover τ(R) = 1 + 1/(1+(R₀/R)^m) | 2 | 8.8% |

**τ = 2, with nothing fitted at all, does exactly as well as the best fitted
constant.** The observed mean sits 0.65 standard errors from 2. The free-streaming
value is not a number we tuned — it's the number the data wanted.

## Why galaxies still work

Fit the crossover and it lands at **R₀ ≈ 900 kpc**:

| Scale | τ |
|---|---:|
| 10 kpc (a galaxy half-light radius) | 1.000 |
| 100 kpc | 1.000 |
| 300 kpc | 1.000 |
| 900 kpc | 1.488 |
| 1200 kpc | 1.947 |
| 2000 kpc | 2.000 |

So galaxies sit squarely in the dust regime — τ = 1.000, exactly what R10 assumed and
exactly what fits the 149 SPARC rotation curves. **Nothing about the galaxy results
changes.** Clusters sit in the streaming regime at τ ≈ 2. The transition is around a
megaparsec.

The physical picture is simple: inside a galaxy the graviton stream is reprocessed
before it gets anywhere, so it behaves like standing energy. Across a cluster it just
streams, and a stream carries momentum flux — which is what makes τ = 2.

## What this predicts next

1. **Groups sit in between.** Systems around 500–900 kpc should show intermediate
   shortfalls, between 1 and 2. That is a clean, falsifiable test on data that exists.
2. **The residual scatter should track size.** It does, in the right direction:
   rank correlation +0.46 with M500 and +0.43 with R500. Bigger systems sit closer to
   2, as a crossover requires.
3. **Cluster weak lensing is a separate, independent test.** τ also fixes the ratio of
   lensing mass to hydrostatic mass, as (1 + 1/τ)/2. At τ = 2 that is **0.75** — and
   the observed hydrostatic bias runs the other way, so this is a real, sharp check
   that has not been run here.

## What this does not do

* It does not fix the energy budget. The companion still costs far more energy than
  starlight supplies; τ = 2 changes how hard that energy pulls, not how much there is.
* Eleven clusters spanning only R500 = 1054–1430 kpc cannot pin the *sharpness* of the
  crossover. They pin the two endpoints, not the shape between them. The m ≈ 10 in the
  fit is unconstrained.
* X-ray hydrostatic masses carry a known bias of order 10–20 percent, which is the
  same size as the 9.5 percent residual scatter. The quoted M500 errors here are only
  2.9–8.5 percent, so the spread is not pure measurement noise — but it is not all
  signal either.
* τ = 2 asserts the companion is genuinely a free-streaming radial flux in clusters.
  This test shows the *number* that implies is the number the data wants. It does not
  independently verify the stress state.

## Reproduction

```
cd code
python cluster_tau.py --output-dir ../run-v1
```

Inputs are the archived `companion_wave_test/cluster_comparison.csv` (X-COP M500,
R500, f_gas) and `companion_deposition_fit/cluster_predictions.csv` (the frozen
galaxy-calibrated transfer). Both are hashed into the output. Requires NumPy and
SciPy; no network access.
