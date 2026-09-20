# RW-1: the whirlpool beyond the rut, backed into regime by regime
20 September 2026. [Protocol](protocol-rw1.md) 17466e7, declared before any calculation; [amendment 1](protocol-rw1-amendment-1.md) beb41d9 (the shell series carries 40 terms, after the synthetic gates and before any data fit). Archive `rw1-results.json`, runtime 644 s; suite job `rw1_checks.py`. Inputs: CL-2 stage 2's caches and archive, the X-COP extract and release correlations, CR-2's lens interface under the static geometry, the Milky Way baselines, Coma's six exposed shear bins.

**The question, in plain language.** The owner asked to let the rut reach: the track an orbit writes keeps pulling far beyond the orbit, catches light, and the pull of all the ruts inside a radius accumulates there. The family carries that reach as an inward pull from every element of ordinary matter, f(s) = s (s^2 + w^2)^(-(p+1)/2), with a reach p (2 is Newton with a soft core, 1 is the whirlpool geometry in which a ring pulls as m/R outside and not at all inside, below 1 the accumulated pull keeps growing outward), a core w (the rut's own width) and a strength lambda (the pull relative to Newton's at 1 kpc). Each regime was fitted on its own, as many formulas as needed, and then the formulas were compared.

**Outcome.** c: not matched even by its own formula: lenses.

| status | finding |
|---|---|
| **reproduction** | **passes.** lambda = 0 with the tabulated force reproduces stage 1's archived baryons-alone galaxies, 52.5609 km/s, to 2.2e-16; stage 2's cluster references under the correlations, 11.6382 and 248.12 per point, to 0.0e+00; NL-1's Milky Way baryons-alone speeds to 0.0e+00; CL-F1's Coma baryons-only shape chi-square, 4.284, to 2.2e-16; stage 1's stars-only lens benchmarks under their own geometry to 1.9e-13. |
| **numerical verification** | **passes every gate.** Rings: the two-dimensional Gauss law at p = 1 to 2.5e-11, the elliptic ring force at p = 2 to 1.2e-10; shells: the shell theorem to 2.5e-11, the closed form at p = 1 to 5.3e-12, adaptive quadrature at 525 points to 2.1e-08 with the two branches agreeing to 2.5e-10; the far field to 1.5e-04. Galaxy columns change by 1.1e-10 under node doubling and 9.5e-06 under source-grid doubling; cluster columns by 3.2e-06 between the 1200-point evaluation and the full grid and 3.4e-06 under source-grid doubling; lens rows by 3.1e-09. Every one-parameter fit's two solvers agree to 1.1e-11 with a worst derivative ratio of 3.0e-13; every NNLS has a KKT residual below 1.5e-15; separability holds to 4.1e-16. |
| **scientific outcome** | **c: not matched even by its own formula: lenses** Exposed data, static geometry for the lenses, mean-field kinematics, the hot-gas field standing in for cluster lensing (no matched shear catalogue): a description within the declared family, not a validated theory. |

## What each galaxy regime asks for (89 training galaxies, reconstructed source)
| regime | one kernel, individual strengths (F2) | speed error | universal (F3) | speed error | baryons alone | individual best reach / core, the mode | at each galaxy's own best kernel |
|---|---|---|---|---|---|---|---|
| G-in, the stars closer in (R < R_half) | p = 2, w = 1 kpc | **16.53 km/s** (Newton kernel 16.53) | p = 1.25, w = 8 kpc, lambda = 0.513 | **29.55 km/s** (Newton kernel 29.60) | 35.6 | 2.0 (36 of 84) / 0.25 (44 of 84) | 15.43 |
| G-out, the stars beyond | p = 1, w = 2 kpc | **8.50 km/s** (Newton kernel 11.50) | p = 1.25, w = 0.25 kpc, lambda = 0.199 | **34.64 km/s** (Newton kernel 35.74) | 62.1 | 1.25 (26 of 85) / 0.25 (19 of 85) | 5.18 |
| G-edge, the last five points | p = 1, w = 0.25 kpc | **4.18 km/s** (Newton kernel 7.14) | p = 1.5, w = 0.25 kpc, lambda = 0.539 | **33.78 km/s** (Newton kernel 35.97) | 72.4 | 0.5 (17 of 89) / 0.25 (26 of 89) | 2.04 |
| G-whole, every point | p = 1, w = 0.25 kpc | **13.94 km/s** (Newton kernel 16.33) | p = 1.25, w = 2 kpc, lambda = 0.208 | **31.49 km/s** (Newton kernel 32.21) | 53.7 | 1.5 (26 of 89) / 0.25 (36 of 89) | 10.32 |

Acceptance 21.9 km/s. "Newton kernel" is the p = 2 member, a mass rescaling, the N1 control. The point-mass whirlpool (N2, the same kernels sourced by each galaxy's mass at its centre) gives 12.67 km/s at p = 0.75, w = 2 kpc against 13.94 for the matter-sourced whirlpool, a difference of -1.27 km/s. With the tabulated force (the sensitivity): F2 p = 1, w = 0.25 kpc 13.03 km/s, F3 p = 1.25, w = 1 kpc lambda = 0.198 30.63 km/s on the whole curves. Galaxies excluded from a regime for having fewer than two points there: inner NGC2683, NGC3992, NGC4013, NGC4138, UGC06614; outer UGC05918, UGC06667, UGC07577, UGC07866.

## The strength against the mass (S2), at each regime's F2 kernel
| regime | kernel | log lambda against log M_b | log lambda_h (at the half-mass radius) against log M_b |
|---|---|---|---|
| inner | p = 2, w = 1 kpc | slope -0.39 [-0.46, -0.32], scatter 0.37 dex, 73 systems (11 at zero) | slope -0.39 [-0.46, -0.32], scatter 0.37 dex, 73 systems (11 at zero) |
| outer | p = 1, w = 2 kpc | slope -0.61 [-0.65, -0.59], scatter 0.19 dex, 85 systems (0 at zero) | slope -0.43 [-0.48, -0.39], scatter 0.30 dex, 85 systems (0 at zero) |
| edge | p = 1, w = 0.25 kpc | slope -0.51 [-0.55, -0.47], scatter 0.22 dex, 89 systems (0 at zero) | slope -0.33 [-0.39, -0.28], scatter 0.34 dex, 89 systems (0 at zero) |
| whole | p = 1, w = 0.25 kpc | slope -0.51 [-0.57, -0.47], scatter 0.29 dex, 89 systems (0 at zero) | slope -0.33 [-0.40, -0.27], scatter 0.40 dex, 89 systems (0 at zero) |
| clusters | p = 2, w = 1 kpc | slope -0.35 [-0.44, -0.29], scatter 0.05 dex, 12 systems (0 at zero) | slope -0.35 [-0.44, -0.29], scatter 0.05 dex, 12 systems (0 at zero) |

## The clusters (twelve X-COP fields under the SZ correlations)
One kernel with individual strengths (F2): p = 2, w = 1 kpc, **19.02 per point** against the release reference's 11.64 (Newton kernel 19.02, baryons alone 248.1); universal (F3): p = 2, w = 1 kpc, lambda = 5.97, 24.52 per point (Newton kernel 24.52). The individual best reach: 2.0 (12 of 12).

| cluster | F2 chi-square per point | F2 strength | own best kernel | own best per point |
|---|---|---|---|---|
| A1644 | 12.8 | 7.54 | p = 2, w = 16 kpc | 12.5 |
| A1795 | 8.3 | 6.93 | p = 2, w = 16 kpc | 7.9 |
| A2029 | 16.7 | 5.78 | p = 2, w = 16 kpc | 16.6 |
| A2142 | 13.7 | 4.85 | p = 2, w = 64 kpc | 11.6 |
| A2255 | 33.0 | 4.54 | p = 2, w = 1 kpc | 33.0 |
| A2319 | 44.5 | 5.2 | p = 2, w = 1 kpc | 44.5 |
| A3158 | 15.2 | 6.05 | p = 2, w = 1 kpc | 15.2 |
| A3266 | 15.1 | 7 | p = 2, w = 1 kpc | 15.1 |
| A644 | 15.0 | 6.82 | p = 2, w = 1 kpc | 15.0 |
| A85 | 12.1 | 6.62 | p = 2, w = 1 kpc | 12.1 |
| RXC1825 | 25.1 | 7.62 | p = 2, w = 1 kpc | 25.1 |
| ZW1215 | 11.6 | 7.19 | p = 2, w = 1 kpc | 11.6 |

At the galaxies' whole-curve kernel p = 1, w = 0.25 kpc the clusters with their own strengths score 64.01 per point; their strengths there are a median 3.3 times the galaxy line extrapolated to their masses; at their own kernel p = 2, w = 1 kpc the ratio to the same line is 3165.3 (different kernels, a stated comparison of magnitudes).

**The lensing curves the fields imply** (alpha(b) divided by the release NFW reference's, median over the twelve clusters; the repository holds no matched cluster shear catalogue, so these are the curves the fitted fields imply, not a fit to lensing): at b = 100 / 300 / 1000 kpc the baryons alone give 0.10 / 0.13 / 0.15; the F2 fields 0.78 / 0.89 / 1.11; U1 2.79 / 4.21 / 8.22; U2 0.26 / 0.37 / 0.72; U3 0.76 / 0.88 / 1.11.

## The lenses (six, static geometry, Chabrier bracket)
At the galaxies' whole-curve kernel, the strength that meets each Einstein radius exactly and the stellar motions it then gives (total kinematics 3965.6 against 128 allowed, not matched); and each lens's own best kernel on the grid.

| lens | Einstein strength at the galaxy kernel | kinematics | anisotropy | own best kernel | its kinematics | stellar mass (Msun) |
|---|---|---|---|---|---|---|
| J0037-0942 | 0.288 | 1000.8 | 0.45 | p = 1.75, w = 0.25 kpc | 19.6 | 2.43e+11 |
| J1112+0826 | 0.28 | 270.8 | 0.45 | p = 1.25, w = 0.25 kpc | 55.4 | 2.37e+11 |
| J1204+0358 | 0.295 | 930.0 | 0.45 | p = 2, w = 0.25 kpc | 7.4 | 1.29e+11 |
| J1402+6321 | 0.249 | 927.8 | 0.45 | p = 1.75, w = 0.25 kpc | 46.4 | 2.85e+11 |
| J1621+3931 | 0.3 | 308.4 | 0.45 | p = 1.75, w = 0.25 kpc | 52.4 | 2.23e+11 |
| J1630+4520 | 0.19 | 527.8 | 0.45 | p = 1.5, w = 0.5 kpc | 89.6 | 3.22e+11 |

Across the six lenses the Einstein strengths against the stellar mass: slope -0.37 [-1.19, -0.14], scatter 0.05 dex, 6 systems (0 at zero).

## The universal candidates, scored everywhere without retuning
| candidate | kernel | strength law | galaxies train / validation / test (km/s) | inner / outer / edge | clusters per point | lenses | Milky Way I / II (km/s) | Coma shape chi-square, brackets 1 / 2 | passes |
|---|---|---|---|---|---|---|---|---|---|
| U1 | p = 1.25, w = 2 kpc | lambda = 0.208 | 31.49 / 34.29 / 28.81 | 30.10 / 34.77 / 35.43 | 6548.70 | 14696.5 and 51.8% | 29.0 / 15.0 | 5.92 / 5.92 | none |
| U2 | p = 1, w = 0.25 kpc | lambda = 10^4.59 M^-0.51 | 21.32 / 29.68 / 17.28 | 26.10 / 19.47 / 19.42 | 153.23 | 21007.5 and 62.9% | 19.4 / 11.9 | 6.96 / 6.59 | galaxies |
| U3 | p = 2, w = 1 kpc | lambda = 10^5.69 M^-0.35 | 726.57 / 784.39 / 772.81 | 681.76 / 770.32 / 621.83 | 20.64 | 1008130.7 and 1973.8% | 1212.2 / 1155.2 | 4.28 / 4.28 | clusters |

Acceptances: galaxies 21.9 km/s on training; clusters 23.28 per point; lenses 128 and 3%. Milky Way baryons alone 52.6 / 62.3 km/s; Coma baryons alone 4.28 on both brackets, CL-F1's written field 9.90, CL-1's NFW reference 3.86. Regimes matched by their own formula: inner yes, outer yes, edge yes, whole yes, clusters yes, lenses no.

## What is shared and what is not
**What each regime backs into.** The stars at the edges and the outer halves of the galaxies ask for the whirlpool geometry exactly: a reach of 1, the pull falling as one over the distance, so that a ring of written matter pulls as m/R outside itself and not at all inside, with a small core (0.25 kpc at the edge, 2 kpc for the outer half) and each galaxy its own strength. The stars closer in ask for nothing beyond Newton: a reach of 2 with a 1 kpc core, which is a mass rescaling galaxy by galaxy, and their individual best reaches split between the two ends of the grid, so the inner curves do not measure the reach at all. The whole curves take the edge's answer, a reach of 1 with a 0.25 kpc core. The clusters ask for the opposite: a reach of 2, Newton's pull with about six times the baryonic mass under one strength, and they cannot use the galaxies' reach (64 per point at the galaxies' kernel even with their own strengths). The lenses are matched by no member: at the galaxies' kernel the strengths that meet the Einstein radii ruin the stellar motions, and each lens's own best kernel, a reach of 1.25 to 2 with a core of 0.25 to 0.5 kpc, still totals above the limit. The Milky Way, never fitted, takes the galaxies' law with its mass scaling and lands at 19 and 12 km/s on the two baselines against 53 and 62 for its baryons alone. Coma's shear shape gets worse with every candidate that adds anything, as CL-F1's written field did.

**What is shared.** Wherever the galaxies carry an extra pull, the common kernel's reach is 1: the outer halves, the edges and the whole curves, with the reconstructed and with the tabulated force alike (galaxy by galaxy the best reaches scatter across the grid, with modes of 1.25, 0.5 and 1.5, because a single curve constrains the reach weakly; the regime as a whole does not). In every regime that measures the reach, the strength fitted galaxy by galaxy falls with the galaxy's baryonic mass as one over its square root: the slope of log strength against log mass is -0.51 for the whole curves and the edges and -0.61 for the outer halves, each with a bootstrap interval of a few hundredths and a scatter of 0.2 to 0.3 dex. Put together, a reach of 1 with a strength falling as the root of the mass is a pull proportional to the root of the mass divided by the distance, which is the deep form of the root law that NL-1 fitted as a response and PM-1 fitted as a local law. With the line of the whole curves, the far pull of a galaxy of mass M is 10^4.59 G M^(0.486)/(R r_*): at the slope -0.51 this is the root law sqrt(a_0 G M)/R with a_0 = 10^(2b) G M^(2q+1)/r_*^2, whose residual mass dependence is M^-0.027; at 10^10 Msun a_0 = 3408 (km/s)^2/kpc, 1.1e-10 m/s^2 (the outer line gives 1.5e-10 m/s^2, the edge line 1.2e-10 m/s^2), against PM-1's fitted a* of 2104 (km/s)^2/kpc, 6.8e-11 m/s^2, and the intercept's 68% interval spans a factor 10.2. The whirlpool premise therefore backs into the same formula from a new side, and this time it carries its own scale.

**What is not shared.** A universal strength fails the galaxies (31.5 km/s, as every universal linear family has), and the point-mass control fits the whole curves 1.3 km/s better than the matter-sourced whirlpool, so the rotation curves constrain the reach and the strength but not where the ruts are written. The clusters share neither the reach nor the line: their strengths at their own kernel sit three orders of magnitude off the galaxy line extrapolated to their masses, and the cluster line applied to the galaxies is absurd (727 km/s). The lenses share the clusters' preference for a Newton-like reach and fail the stellar motions anyway. The lensing curves the cluster fields imply follow the fits: between 100 kpc and 1 Mpc the clusters' own whirlpool fields bend light at 0.78 to 1.11 of the release reference's deflection (the baryons alone 0.10 to 0.15), the galaxies' law at 0.26 to 0.72, the universal galaxy strength at 2.8 to 8.2; no matched shear catalogue exists to judge these, and they are what the hot-gas field implies under the static coupling. The one law that carries anything across systems, the galaxies' kernel with the root-of-mass strength (U2), passes the training galaxies by 0.6 km/s, misses the validation galaxies as NL-1 did, fits the Milky Way without being asked, and fails the clusters and the lenses.

## What is not shown
The reach and core are read off a declared grid, not fitted continuously; the strengths are fitted on exposed data with no energy budget, formation or propagation behind them; the sideways part of a swirl is TF-1's vector law and is not carried; the cluster lensing curves are implied by the hot-gas field under the static well coupling, not measured; the galaxy sources are stage 2's razor-thin reconstruction; the lens kinematics are mean-field with constant anisotropy; the Milky Way and Coma are transfers without acceptance; no held-out system.

## Next decision
The galaxies have now backed into the same formula three times, from three premises: a local root law (PM-1), a root response of the written field (NL-1) and a whirlpool of reach 1 with a root-of-mass strength (this experiment), the last carrying an acceleration scale that agrees with the fitted local law within the line's uncertainty. No premise has moved the clusters or the lenses onto that formula, and the whirlpool moves the clusters away from it. The next question is not another kernel. It is what sets the strength of the pull by the system it acts in, so that galaxies, the Milky Way and clusters read off one rule: an energy budget of the written field is the natural candidate, and the six lenses, which no static family has passed, are the test that decides whether the obstruction is in the family or in the lens model (its stellar mass bracket, its constant anisotropy, the static geometry).
