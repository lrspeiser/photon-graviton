# RW-1 protocol: the whirlpool beyond the rut, backed into regime by regime
Declared 20 September 2026, before any calculation on data. Owner's request (20 September): assume that the rut
written by an orbiting body not only deepens the well along its path but also creates a swirling whirlpool of gravity
far beyond it, which catches light and lenses it; that each rut adds to the speed, so that the entrained gravity
increases with distance from the centre; back into the formula from the lensing curves of clusters, from the stars at
the edges of galaxies and from the stars closer in, as many formulas as needed to match each individually, and then
find the behaviours the formulas share.

## The question, in plain language
Every field tested so far was written where the matter is (a blurred copy of it, CL-2), or grew as the root of the
ordinary pull (NL-1), or was written a fixed distance away (NK-1). This experiment lets the rut reach: the track an
orbit writes keeps pulling far outside the orbit, and the pull from all the ruts inside a radius accumulates there.
Whether the far reach is a swirl or a plain pull, what it can add to lensing and to orbital speeds is its inward
component, so that is what the family carries; the sideways part of a swirl is TF-1's vector law, which stays
obstructed for light by the factor c/v it measured, and a steady sideways push on stars has no static counterpart.
The question is then: what reach, what core and what strength do the inner stars, the outer stars, the cluster
fields and the lenses each ask for, and do they ask for the same thing?

## The family (formula provenance)
Each element of ordinary matter dm at r' pulls a body at r toward r' with

    g_whirl(r) = lambda G int dm f_{p,w}(|r - r'|),   f_{p,w}(s) = s (s^2 + w^2)^(-(p+1)/2) r_*^(p-2),   r_* = 1 kpc.

- p is the reach: p = 2 is Newton's pull with a Plummer core (established, the control, a pure mass rescaling);
  p = 1 is the force of the two-dimensional logarithmic potential, the whirlpool geometry in which a ring's pull
  outside it is m/R and inside it nothing (established); p < 1 makes the accumulated pull grow with radius without
  bound (the owner's "entrained gravity increases as you move out"); the intermediate values are proposals.
- w is the rut's own core: inside w the pull is softened (the deepened well along the path); the grid of w is the
  RUT-1 footprint scale.
- lambda is the strength, dimensionless: the whirlpool pull relative to Newton's at 1 kpc from a point mass. For each
  system the size-referred strength lambda_h = lambda (R_h/r_*)^(2-p), the pull relative to Newton's at the system's
  own half-mass radius R_h, is also reported.
- Circular speeds: v^2 = R (g_N + g_whirl) with g_N the reconstructed Newtonian force of stage 2 (the same components
  that write the ruts; the tabulated force is the sensitivity). Light: the static Euclidean route PF1 with well coupling,
  alpha(b) = (4/c^2) int (g_N + g_whirl)(b/r) dz, the three-dimensional quadrature.
- The ring and shell averages of f_{p,w} are new elementary integrals: the ring average by 128-node Gauss-Legendre
  quadrature in the azimuth; the shell average in closed form (two power integrals in u = r^2 + r'^2 + w^2 - 2 r r' mu,
  differences taken through expm1 and log1p) where x = 2 r r'/(r^2 + r'^2 + w^2) > 0.2 and by the binomial series in x
  (20 terms) elsewhere. Kernels are declared on the grid p in {0.5, 0.75, 1, 1.25, 1.5, 1.75, 2}; galaxy and lens cores
  w in {0.25, 0.5, 1, 2, 4, 8, 16} kpc; cluster cores w in {1, 2, 4, 8, 16, 32, 64, 128, 256} kpc.

## Data and regimes, each backed into on its own
- G-in, the stars closer in: the points of each of the 89 training galaxies with R < R_h, the radius enclosing half
  the galaxy's baryonic mass on the stage-2 reconstruction (median 2.2 disk scale lengths); galaxies with fewer than
  two such points are excluded from this regime and listed.
- G-out, the stars beyond: R >= R_h, the same rule.
- G-edge, the stars at the edge: the last five points of each galaxy (all points when there are fewer).
- G-whole: every point (the standing acceptance metric).
- C, the cluster fields: the twelve X-COP clusters' pressure profiles under the release's SZ correlations with the
  corrected boundary fit, exactly stage 2's block. Under the static well coupling the same field bends light, so the
  cluster's lensing curve alpha(b) is what its fitted field implies; the repository holds no matched cluster shear
  catalogue, and this is stated wherever those curves appear. Coma's six exposed, figure-reconstructed shear bins are
  a shape test (one nonnegative nuisance, CL-F1's construction), transfer only.
- L, the six lens galaxies (static geometry, Chabrier bracket, constant anisotropy in [-2, 0.45]): transfer of the
  galaxy formulas, and per lens the strength that meets its Einstein radius exactly, read against its stellar motions.
- MW, the Milky Way baselines I and II (38 bins): transfer only.
Every fit below uses stage 2's galaxy loss (the equal-galaxy speed error, convex in lambda on its domain), the
cluster whitened chi-square with the boundary nuisance, and the lens rows of stage 2 with the three-dimensional bend.

## The fits
- F1, individual: for every system, regime and kernel (p, w), lambda >= 0 by a certified one-parameter solve
  (bisection on the analytic derivative, checked by a bounded Brent minimisation of the value; convexity is analytic,
  F'' = sum w v R^2 W^2/(2 u^(3/2)) > 0). The best kernel per system is read off the grid: (p_i, w_i, lambda_i).
- F2, one kernel per regime with individual strengths: (p, w) minimising the equal-galaxy loss with each galaxy's own
  lambda_i; the problem is separable, so the optimum is the sum of the individual optima (gate C2).
- F3, universal: one (p, w, lambda) for the whole regime, the certified one-parameter solve at every kernel.
- Clusters: F1 per cluster per kernel by nonnegative least squares on the whitened rows with the boundary nuisance;
  F2 the common kernel with per-cluster strengths; F3 the universal strength (13-variable NNLS: one lambda, twelve
  nuisances). KKT residuals recorded.
- Lenses: per lens and kernel, lambda_l from the Einstein condition (linear; zero if negative), then the kinematics
  chi-square at lambda_l with the anisotropy fitted; and the transfers of the universal candidates.
- Sensitivity: F2 and F3 for the galaxy regimes with the tabulated force.

## The shared behaviours (declared statistics, no thresholds)
- S1: the distribution of the individual best (p_i, w_i) in each regime, and each regime's F2 kernel.
- S2: in each galaxy regime at its F2 kernel, the slope and intercept of log lambda_i against log M_i (baryonic mass)
  by least squares, with a 68% bootstrap interval over galaxies (1000 resamples, seed 7), and the scatter in dex; the
  same for lambda_h.
- S3: the clusters' lambda_c at their F2 kernel against the galaxy line of G-whole extrapolated to each cluster's
  baryonic mass: the ratio per cluster and its median; the same at the galaxies' F2 kernel applied to the clusters.
- S4: the universal candidates. U1: F3 of G-whole. U2: the F2 kernel of G-whole with lambda set by the S2 line,
  lambda = lambda_0 (M/M_0)^q (no fit to any system beyond the line). U3: the clusters' F2 kernel with the clusters' own
  strength line. Each is scored on every regime, the validation and test galaxies, the Milky Way, the six lenses and
  Coma, without retuning.
- The cluster lensing curves: for each cluster, alpha(b) at b = 100, 300 and 1000 kpc from its F2 field and from
  U1, U2, U3, divided by the release NFW reference's alpha(b); the Newtonian parts are truncated at the cluster's grid
  edge r_out for model and reference alike, the whirlpool part is extended analytically.

## Negative controls
- N1, p = 2: the Newtonian kernel is a mass rescaling; its results are reported beside every regime's. A whirlpool
  geometry is claimed for a regime only if its F2 reach differs from 2 by more than one grid step.
- N2, the point-mass whirlpool: the same kernels sourced by each galaxy's total baryonic mass placed at its centre,
  F2 on G-whole; if it matches within 1 km/s of the matter-sourced F2, the rotation curves do not test where the ruts
  are written, and the report says so.
- N3, lambda = 0: the baryons alone in every regime.

## Numerical gates (thresholds set from synthetic prototypes, before any data fit)
- K1 ring, p = 1, w = 1e-6: two-dimensional Gauss (m/R outside, 0 inside) at R/R' in {0.3, 0.7, 0.9, 1.1, 1.5, 3, 30},
  absolute error below 1e-9 in units of 1/R'.
- K2 ring, p = 2, w = 1e-6, against the exact elliptic ring force: relative error below 1e-8 at the same points.
- K3 shell, p = 2, w = 1e-6: the shell theorem, error below 1e-8 in units of 1/R'^2.
- K4 shell, p = 1, w = 1e-6, against the closed form 1/(2r) + (r^2 - r'^2)/(4 r^2 r') ln((r + r')/|r - r'|): relative
  error below 1e-9.
- K5 the two-branch shell kernel against adaptive quadrature (relative tolerance 2e-14) at 7 reaches x 5 cores x 15
  radii: error below 1e-6 relative to max(|K|, C^(-a) min(r/r', 1)); the two branches agree to 1e-8 where 0.15 < x < 0.3.
- K6 far field at 100 core radii, rings and shells, every p: |K r^p - 1| below 1e-3.
- K7 galaxy columns: 128 against 256 nodes on NGC2403, NGC3198 and DDO154, relative change below 1e-8 for every kernel;
  3000 against 6000 source points, below 1e-3.
- K8 cluster columns: the 1200-point evaluation with the cubic spline in ln r against the full 6000-point evaluation on
  A1795 and A2255 at w in {1, 16, 256} kpc and every p, relative change of the pressure columns below 1e-5; 6000 against
  12000 source points on the same, below 1e-3.
- K9 lens columns: the 1000-point evaluation with the spline against the full 4001-point evaluation on J0037-0942 at
  w = 1 kpc and every p, relative change below 1e-5.
- C1 every one-parameter fit: the two solvers agree in the loss to 1e-9 relative; at an interior optimum the derivative
  is below 1e-8 of its value at lambda = 0; every NNLS has a KKT residual below 1e-10 relative.
- C2 separability: F2's loss at the assembled strengths equals the sum of the individual minima to 1e-12.
- R1 galaxies: lambda = 0 with the tabulated force reproduces stage 1's archived baryons-alone reference,
  52.5609 km/s, to 1e-9; stage 2's caches load by hash.
- R2 clusters: lambda = 0 reproduces stage 2's newtonian_baryons, 248.1196 per point, and the release reference,
  11.6382, under the correlations to 1e-9.
- R3 Milky Way: lambda = 0 reproduces NL-1's archived baryons-only speeds, 52.5749 and 62.3402 km/s, to 1e-9.
- R4 Coma: the baryons-only shape chi-square reproduces CL-F1's 4.284 on both brackets to 1e-6.
- R5 lenses: the stars-only lens rows reproduce stage 2's archived stars-only kinematics and Einstein residuals per
  lens to 1e-9 where the archive holds them.
A failed gate is recorded with its magnitude and is never re-thresholded; an amendment may add a gate or a replacement
run, never remove a failure.

## Acceptance and decision, fixed beforehand
Galaxies: the equal-galaxy speed error on the 89 training galaxies at or below 21.9 km/s (validation and test
reported); a regime is matched when its own points give at or below 21.9 under the same metric. Clusters: within a
factor two of the release reference, 2 x 11.638 per point. Lenses: total kinematics chi-square at or below 128 and
worst Einstein residual at or below 3%. Milky Way and Coma are reported without acceptance.
- (a) a universal candidate (U1, U2 or U3) passes the galaxies, the clusters and the lenses at once;
- (b) every regime is matched by its own formula (F2 in the galaxy regimes, F2 in the clusters, the Einstein-fixed
  strengths in the lenses with the kinematics limit), but no universal candidate passes all three; the shared
  behaviours are reported as measured, with the trade-offs;
- (c) some regime is not matched even by its own formula on the declared grid: named, with its best value.
Sub-labels record which acceptance fails for each candidate.

## Outputs and reproduction
Archive `rw1-results.json`; report `report-rw1.md`; library `rw1_lib.py`; driver `rw1.py`; caches under
`research_work/generated/routes/rw1_*.npz`; suite job `rw1_checks.py` registered before `rut8_checks.py`. Exposed data
throughout; static Euclidean geometry for the lenses; mean-field kinematics; no dark matter, no halo fitted; a
description within a declared family, not a validated theory.
