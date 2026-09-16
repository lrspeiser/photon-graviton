# PM-2A stage A: the source survives its audit, and PM-1's mass is not the source's mass

Stage A of [protocol-pm2a.md](protocol-pm2a.md), declared in d28b67a before this run. It builds one frozen ordinary-matter source per SPARC galaxy, independently of the rotation speeds being fitted, and gates it before any field is solved. No fit is performed here and no force law is scored.

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

The AQUAL cell-mass construction's total differs from the exact total by 0.9–3.9×10⁻⁵ on four galaxies. That comparison is **reported, not gated** — it decomposes a thickened disk over spherical shells, so only its total is commensurable with a cylindrical integral, and it carries its own quadrature and outer-truncation error. Any later step that needs the cells to be accurate must establish that separately, and stage C's convergence gate is where that gets tested.

## The measurement stage A exists to make

PM-1 fitted a rule in `M_force = R²g_N/G`. Against the source's own cylindrical mass at the last sampled radius:

| Ratio | 10th | median | 90th | range |
|---|---|---|---|---|
| M_force / M_cyl(<R_last) | 1.126 | **1.381** | 1.763 | 1.040 – 3.174 |
| M_force / M_total | 0.915 | **1.177** | 1.397 | 0.670 – 3.009 |

M_force exceeds the cylindrical source mass at every one of the 149 galaxies. That is the expected sign for a disk — a flattened source pulls harder at a given radius than a sphere of the same enclosed mass — and it is the quantitative statement that PM-1's fitted rule cannot be read as "a force sourced by the matter inside r". It is a statement about geometry and about whether the force law holds, not a defect in the source model, and it is what stages B and C have to account for by solving a field rather than substituting a mass.

**But the ratio is not pure geometry, and should not be quoted as if it were.** `v_bar` is SPARC's tabulated ordinary-matter contribution, while M_cyl is the mass of the source reconstructed here; RPG-1 already documents that its reconstructed density and the archived rotation contributions differ materially at some radii. The 1.381 therefore mixes disk flattening with that reconstruction difference, in unknown proportions. Separating them needs the source's own Newtonian field, computed on the same grid — which is exactly why stage C recomputes its baseline inside the experiment instead of comparing against PM-1's archived scores. Until then the honest reading is: the two masses differ by tens of percent, with the right sign for a disk, and the split between the two causes is not yet measured.

**M_sph(<r) is deliberately not quoted.** The disks' spherical enclosed mass is not defined by a surface density alone; it depends on the declared vertical profile. Stage B verifies the spherical law on a genuinely spherical source instead of reporting a cylindrical integral under a spherical label.

## What this does and does not establish

**Shown.** One frozen, independently specified source per galaxy, integrated exactly, that is invariant to which velocity rows are kept — with a control proving the invariance test has teeth. The size of the extrapolated outer mass, per galaxy. The size of the error in the quadrature the archive relies on, measured on every galaxy rather than a sample. The gap between PM-1's force-equivalent mass and a real source integral.

**Not shown.** Nothing about any field, any fit, or any force law: stage A solves nothing and scores nothing. The archived quadrature's 2.4×10⁻⁷ is a diagnostic of that routine, not a revision of RPG-1's published results, which are untouched. And the gates are internal-consistency gates; passing them says the source is well defined, not that it is right.

## Reproduce

```sh
python research_work/results/path-memory/pm2a.py
```
