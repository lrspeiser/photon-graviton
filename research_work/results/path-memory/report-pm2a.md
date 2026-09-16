# PM-2A stages A and B: the source survives its audit, PM-1's mass is not the source's mass, and the readout is the accuracy bottleneck

Stage A of [protocol-pm2a.md](protocol-pm2a.md), declared in d28b67a; the two source-geometry corrections and stage B declared in fa8b407, both before the runs they govern. It builds one frozen ordinary-matter source per SPARC galaxy, independently of the rotation speeds being fitted, and gates it before any field is solved. No fit is performed here and no force law is scored.

**All four gates pass, on every galaxy each one applies to, with no construction failures.** The source builds for all 149; the row-removal gate covers the 131 with at least eight velocity rows.

| Gate | Requirement | Measured |
|---|---|---|
| Mass-preserving refinement | conserve the source mass to 10⁻¹⁰, no clipping | **2.2×10⁻¹⁶** |
| Row-removal invariance, end to end | predicted force unchanged to 10⁻¹⁰, source hash unchanged | **0.0**, 131 galaxies, every hash unchanged |
| — its positive control | a source rebuilt from shortened rows must be *caught* | **1.4×10⁻⁴ to 2.3×10⁻²**, every hash changed |
| Integral agreement | closed form vs an independently converged quadrature, 10⁻⁸ | **6.6×10⁻¹³** (quadrature self-convergence 1.8×10⁻¹⁵) |
| Collective tail algebra | removing the tail leaves M_cyl(<r) alone and scales the term by 1/√(1−f) | **2.2×10⁻¹⁶** |

## What the source is

RPG-1's archived construction, kept intact and integrated exactly. ln Σ is interpolated linearly in radius over the retained SBdisk samples, so the stellar cylindrical mass is analytic segment by segment; beyond `R_star_source_end` the endpoint-anchored continuation closes to 2πΣ(R_end)·rd·(R_end + rd); the gas disk's total is its own normalization; the bulge is the existing prescription. The primary model keeps the continuation, the declared sensitivity removes only the extrapolated stellar tail, and **neither is selected by which fits better** — that choice was made before this run and is not revisited by it.

This is a source independent of the observed rotation speeds, not one free of force-model inputs: the gas scale is fitted to the tabulated gas contribution and the bulge is reconstructed from the tabulated bulge contribution.

**How much of the source is extrapolated.** The stellar tail is a median **1.06%** of the total source mass, 8.0% at the 90th percentile, and at most 28.1% (UGC06667, whose SBdisk stops at 7.85 kpc against rd = 5.15 kpc). As a fraction of the stellar mass alone the median is 2.77%, the 90th percentile 15.9% and the maximum 57.1%. Any stage C conclusion that moves between the two prescriptions will be reported as outer-source dependent.

## The gate that needed a positive control

Row-removal invariance is trivially satisfied by a pipeline that never consults the mask — which is exactly what a correct pipeline looks like, and also what a vacuous test looks like. So the gate carries a positive control: for eight galaxies the source is deliberately rebuilt from a shortened rotmod array, the mistake the gate exists to catch. Every rebuilt source changes its hash, and its predictions move by 1.4×10⁻⁴ to 2.3×10⁻², six to eight orders of magnitude above the 10⁻¹⁰ tolerance and above the measured 0.0. The gate can see the failure it is meant to see.

The reason the rebuild is wrong, rather than merely different, is that the rotmod array sets three things at once: the stellar profile, the gas scale and the bulge. Truncating it moves the source itself — one galaxy's total mass changes by 8.6% — which is not what "remove some velocity measurements" should mean.

## Exact integrals, and what the archived quadrature actually carries

Across all 149 galaxies the closed forms agree with an independently converged Gauss–Legendre integration of the same profile to 6.6×10⁻¹³, and that quadrature's own 8-node-to-16-node refinement moves by at most 1.8×10⁻¹⁵, so its convergence is shown rather than asserted.

Against those exact totals, RPG-1's `_cylinder_mass` — 20,001 geometrically spaced points under a trapezoidal sum — differs by a median **2.4×10⁻⁷** and at most 4.5×10⁻⁷ across all 149 galaxies, consistent with the 1.2×10⁻⁸–1.4×10⁻⁷ it carries on a unit-total exponential. It is kept as a regression diagnostic and is not used as a reference: two quantities produced by that route agreeing cannot establish 10⁻⁸ accuracy.

The AQUAL cell-mass construction's total differs from the exact total by 0.9–3.9×10⁻⁵ on four galaxies. That comparison is **reported, not gated** — it decomposes a thickened disk over spherical shells, so only its total is commensurable with a cylindrical integral. It is *not* mostly quadrature and outer truncation, as I had implied: correction 9 below accounts for most of it. Replacing the inherited inner mass with the spherical integral drops the worst discrepancy from 3.9×10⁻⁵ to 5.3×10⁻⁶, and ESO444-G084's from 1.4×10⁻⁵ to 4.2×10⁻⁸.

## Two source-geometry corrections

The owner's review of e9e931d found two places where the same declared source was integrated over the wrong volume. Neither changes a fitted parameter or the outer continuation.

**Correction 8: `m_cyl` added a spherical bulge mass to a cylindrical integral.** `m_bulge(R)` is an enclosed mass inside a *sphere*; the disk terms are cylinder integrals. The correct projection adds the portions of exterior shells that fall inside the cylinder, each shell at r > R contributing the fraction 1 − √(1 − R²/r²). On a Plummer sphere, where the closed form is M R²/(R²+b²), the implemented projection reproduces it to **1.1×10⁻¹⁶** at R/b = 0.5, 1, 2 and 5, and refining 16 → 32 nodes moves nothing. Substituting the spherical mass understates the cylinder by **29.2893%** at R = b and 55.3% at R = b/2.

It also corrects a false sentence in stage A's archive. In spherical symmetry M_force equals M_**sph**, by Newton's theorems — not M_cyl, which also contains material outside the sphere, so M_cyl ≥ M_sph whenever the source extends past the evaluation radius.

**Measured impact on the numbers already published: none.** In all 31 bulge galaxies the bulge's tabulated support ends exactly at the outermost rotmod radius, so at R_last the cylinder already contains the whole bulge. At half the bulge support the median change is 0.000%, the 90th percentile 0.005%, the maximum 0.120% (UGC11914) — `np.maximum.accumulate` saturates the adopted bulge mass well inside those radii. The fix is made because the expression was wrong, not because it moved a number.

**Correction 9: the thick-disk builder hands the solver a thin-disk inner mass.** `exponential_disk_cell_masses` returns Σ(0)πr_in² for the unresolved core — a cylinder's mass at constant Σ — while the solver consumes it as the mass inside the inner *sphere*. The exact spherical integral for a varying Σ is 2π∫₀^{r_in} R·Σ(R)[1 − e^(−√(r_in²−R²)/h)] dR, smaller by a leading ratio 3h/(2r_in). Measured on four galaxies at h = 0.1·rd: the inherited value is **15.57 to 15.60×** the spherical one against a leading 15.0, and that inner region is only 0.9–3.8×10⁻⁵ of the total — which is precisely why it accounted for most of the cell-total discrepancy above. Stage B's spherical tests use `spherical_cell_masses` and are unaffected.

## The measurement stage A exists to make

PM-1 fitted a rule in `M_force = R²g_N/G`. Against the source's own cylindrical mass at the last sampled radius:

| Ratio | 10th | median | 90th | range |
|---|---|---|---|---|
| M_force / M_cyl(<R_last) | 1.126 | **1.381** | 1.763 | 1.040 – 3.174 |
| M_force / M_total | 0.915 | **1.177** | 1.397 | 0.670 – 3.009 |

M_force exceeds the cylindrical source mass at every one of the 149 galaxies. That is the expected sign for a disk — a flattened source pulls harder at a given radius than a sphere of the same enclosed mass — and it is the quantitative statement that PM-1's fitted rule cannot be read as "a force sourced by the matter inside r". It is a statement about geometry and about whether the force law holds, not a defect in the source model, and it is what stages B and C have to account for by solving a field rather than substituting a mass.

**But the ratio is not pure geometry, and should not be quoted as if it were.** `v_bar` is SPARC's tabulated ordinary-matter contribution, while M_cyl is the mass of the source reconstructed here; RPG-1 already documents that its reconstructed density and the archived rotation contributions differ materially at some radii. The 1.381 therefore mixes disk flattening with that reconstruction difference, in unknown proportions. Separating them needs the source's own Newtonian field, computed on the same grid — which is exactly why stage C recomputes its baseline inside the experiment instead of comparing against PM-1's archived scores. Until then the honest reading is: the two masses differ by tens of percent, with the right sign for a disk, and the split between the two causes is not yet measured.

**M_sph(<r) is deliberately not quoted.** The disks' spherical enclosed mass is not defined by a surface density alone; it depends on the declared vertical profile. Stage B verifies the spherical law on a genuinely spherical source instead of reporting a cylindrical integral under a spherical label.

## Stage B: the field equations, verified before any galaxy

Four equations, each with its own μ, spherical inverse, log-conductance relaxation weight and field functional, in [fields.py](fields.py). RPG-1's `aqual.py` is subclassed, never modified.

| Equation | μ(x) | spherical x(y) | weight | F(u) |
|---|---|---|---|---|
| Newtonian | 1 | y | 1 | u |
| **Completion I** | 4x/[1+√(1+4x)]² | y + √y | √(1+4x)/[1+√(1+4x)] | t³(t + 2/3) |
| **Completion II** (auxiliary ψ) | x | √y | 1/2 | (2/3)u^{3/2} |
| simple-μ comparison | x/(1+x) | ½(y + √(y²+4y)) | (1+x)/(2+x) | u − 2√u + 2ln(1+√u) |

**The algebra, before any solve.** μ(x)·x = y holds to ≤ 8.9×10⁻¹⁶ for every equation over 10⁻¹⁰ ≤ y ≤ 10¹⁰, ν is consistent with x(y) to ≤ 4.4×10⁻¹⁶, and each weight matches a numerical d ln μ/d ln x to ≤ 3.7×10⁻¹⁰. Each F′ matches 2xμ(x) to ≤ 8.5×10⁻⁹, all finite-difference limited. The two completions' functional densities agree identically on the matched spherical solution, F_I(y+√y) = y² + (2/3)y^{3/2}, to **1.3×10⁻¹⁵** over twenty decades — an energy check independent of the acceleration.

**Cancellation.** Completion I's t-form subtracts no nearly equal numbers at all: at x = 10⁻¹² it returns 6.6666666667×10⁻³⁷ against an exact (2/3)x³. The inherited simple-μ functional, written as `y − 2√y + 2·log1p(√y)`, returns **exactly 0.0** at x = 10⁻⁸ where the truth is 6.667×10⁻²⁵; the series branch returns it correctly.

**The solves.** Each equation solved on a Plummer sphere at compactness GM/(a\*b²) = 10⁻², 1 and 10², so weak, transition and strong regimes are all exercised:

| Check | Tolerance | Worst measured |
|---|---|---|
| Interior flux identity r²μ(g/a\*)g = GM(<r), imposed faces excluded | 10⁻⁶ | **4.4×10⁻⁷** |
| Against the analytic law g_N + √(a\*g_N) | 10⁻⁶ | **2.3×10⁻⁷** |
| The two completions agreeing with each other | — | 2.3×10⁻⁷ |
| Equation residual off the gauge cell | 10⁻⁸ | 1.8×10⁻¹² |
| Iteration change | 10⁻¹¹ | 3.4×10⁻¹² |

Convergence requires **both** the residual and the iteration change; `aqual.py` sets its flag from the change alone. Completion II is two separate solves whose accelerations are added as vectors, and its auxiliary μ is never evaluated on the total gradient. At g_N = a\* the completions give 2.000000 and the simple-μ law 1.618034 = (1+√5)/2 — they must agree with each other there, and a test requiring them to differ would itself be wrong.

**The finding: the readout, not the solve, is the accuracy bottleneck.** The protocol required off-grid readout to be tested separately from face flux, because they are different numerical paths. They are, by three to four orders of magnitude:

| Readout path | error at nr=300 | nr=600 | nr=1200 | order |
|---|---|---|---|---|
| the face gradients themselves | 1.1×10⁻⁷ | 1.2×10⁻⁷ | 3.1×10⁻⁷ | — |
| inherited `gradient_at` (cell-centre averages) | 1.65×10⁻⁴ | 4.11×10⁻⁵ | 1.02×10⁻⁵ | 2 |
| inherited `midplane_speed` (faces, linear in ln r) | 1.61×10⁻⁴ | 4.02×10⁻⁵ | 1.01×10⁻⁵ | 2 |
| **face gradients, cubic spline in ln r** | **1.26×10⁻⁸** | 9.3×10⁻¹⁰ | 7.3×10⁻¹¹ | ~4 |

Both inherited paths average or linearly interpolate, costing O(h²) on a solution whose faces are already accurate to 10⁻⁷. Only the face spline is gated (worst 2.2×10⁻⁸); the inherited numbers are published as a measurement, not used. This matters directly for stage C: a difference between the two completions smaller than the readout error is unresolved, and at 1.6×10⁻⁴ the inherited readout would have set that floor a hundred times higher than it needs to be. It is not a defect in RPG-1's published rotation curves — 1.6×10⁻⁴ on g is far below their observational errors — but it is not good enough for differencing two completions.

Stage B establishes an equation implementation on a spherical source. It is not path memory, and the solver's relaxation time is a numerical quantity that is never a physical memory time.

## What this does and does not establish

**Shown.** One frozen, independently specified source per galaxy, integrated exactly, that is invariant to which velocity rows are kept — with a control proving the invariance test has teeth. The size of the extrapolated outer mass, per galaxy. The size of the error in the quadrature the archive relies on, measured on every galaxy rather than a sample. The gap between PM-1's force-equivalent mass and a real source integral. Four field equations implemented consistently down to their functionals and verified on a spherical source to 4.4×10⁻⁷, with the two completions agreeing where they must and the simple-μ law distinct where it must be. And the size and order of the readout error that would otherwise have limited stage C.

**Not shown.** No fit, no galaxy, no force law scored: stage A audits a source and stage B verifies equations on a sphere. The archived quadrature's 2.4×10⁻⁷ is a diagnostic of that routine, not a revision of RPG-1's published results, which are untouched. And the gates are internal-consistency gates; passing them says the source is well defined, not that it is right.

## Reproduce

```sh
python research_work/results/path-memory/pm2a.py
```
