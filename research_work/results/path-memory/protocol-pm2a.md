# PM-2A: which shared field reproduces PM-1's acceleration law, from an independently specified source

Declared before execution, 15 September 2026. Baseline: `main` at 89ef51d.

**The question.** Starting from one independently specified ordinary-matter distribution, which shared-field completion reproduces the empirical acceleration law PM-1 fitted, and what further predictions follow from its geometry?

**Why it is run.** PM-1's correction 7 established that its fitted rule is g = g_N + √(a\*·g_N), a local acceleration law built on a force-equivalent mass, not a force sourced by matter inside r. Two things follow. The source has to be specified independently of the force it is meant to explain, and the field equation has to be solved rather than assumed. The owner's review supplies both the completions to test and the audit that must precede them.

**What this is not.** Neither completion derives a spectrum, a decay law, phase coherence or path memory. They are controlled static benchmarks against which a later dynamical mechanism can be judged. Memory is tested only after a candidate survives the source and field stages, and RC-2b and the optical candidates stay separate obligations.

## Stage A: separate source mass from force-equivalent mass

Four quantities stay distinct throughout, and every result names which one it used:

    M_force(R) = R² g_N/G        reconstructed from a force
    M_cyl(<R)                    density integrated over a cylinder
    M_sph(<r)                    density integrated over a sphere
    M_total                      including the source's stated outer continuation

The planar ring calculation takes the cylindrical mass, the spherical analytic test takes the spherical mass, and the disk solver takes the density itself rather than either enclosed-mass surrogate. The source model is RPG-1's `baryons.py`: `sparc_components` for the stellar surface density, gas and bulge, `sparc_total_mass` for the independently integrated total, and `milky_way_masses` for the Milky Way. These are model-inferred quantities with stated assumptions, not assumption-free measurements, and the report says so.

**Recomputed with the physical source:** the cumulative rule, the collective rule and the mass–speed slopes. **Preserved unchanged:** PM-1's force-proxy versions, under that label.

### The outer-source convention, after the owner's ruling

**Primary.** PM-2A's primary source model retains RPG-1's archived stellar-profile interpolation and endpoint-anchored exponential continuation, its separately normalized gas disk, and its existing bulge prescription. The continuation parameters are fixed from the archived source inputs and are **not** adjusted using rotation-curve residuals. The stellar endpoint is the final retained positive SBdisk sample in the archived source construction — `R_star_source_end` — and it is not assumed to equal the final original photometric measurement.

**Sensitivity.** A second source model removes only the extrapolated stellar tail, without renormalizing the interior stars, the gas or the bulge. Both source models stay fixed under velocity-row removal, and the same frozen endpoint is used in both. The extrapolated stellar-mass fraction and the extrapolated fraction of the total source mass are reported per galaxy, and every compared force law uses matched source models.

**Interpretation.** Continued mass is model-inferred, not directly observed. Tail removal is a source-assumption sensitivity, not an observational confidence bound. Any conclusion that changes between the two prescriptions is reported as outer-source dependent, and **neither prescription is selected by which fits better**.

**What "independent" means here.** The gas scale is fitted to the tabulated gas contribution and the bulge is reconstructed from the tabulated bulge contribution, so this is a source model independent of the observed rotation speeds being fitted — not one that uses no force-model inputs at all. The gas is not truncated at the stellar endpoint: its normalization comes from a separate input, and cutting it there would change two assumptions at once.

**Exact integrals, and the archived quadrature as a diagnostic.** Because the archived construction interpolates ln Σ linearly in radius, the stellar cylindrical mass is analytic piece by piece, and the continuation closes analytically to 2πΣ(R_end)·rd·(R_end + rd); the gas disk's total is its own normalization. Stage A therefore uses those exact expressions. RPG-1's `_cylinder_mass` — 20,001 geometrically spaced points under a trapezoidal sum — is measured against them and reported as a regression diagnostic: on a unit-total exponential it carries relative errors of 1.2×10⁻⁸ to 1.4×10⁻⁷, so agreement between two quantities produced by that route cannot establish 10⁻⁸ accuracy. Conservation and quadrature accuracy are checked separately: telescoping annular masses conserve by construction, which is not evidence that the underlying cumulative mass is accurate.

**Three gates before any fit, run first on analytic sources and then on the galaxies:**
1. **Mass-preserving refinement.** Refining the annuli must conserve the source mass to 10⁻¹⁰ and must not clip it. The per-ring √N growth is then measured on a genuine density rather than on differenced force-equivalent mass.
2. **Row-removal invariance, end to end.** The requirement is invariance of the *predicted force*, not of a stored mass function. The source is built once from the full frozen source inputs and passed separately from the velocity mask; predictions at the full and retained radii are compared at their common radii with the parameters held fixed, and the source's parameters and hash must be unchanged. `sparc_components` is never rebuilt from a shortened rotmod array, because that array also sets the stellar profile, the gas scale and the bulge. The force-proxy contrast measured earlier stays as a separate diagnostic: a model can show that contrast while its prediction path still hides a dependence on the last evaluation radius, and only the end-to-end test catches that.
3. **Integral agreement.** The closed-form masses must agree to 10⁻⁸ with an independently converged Gauss–Legendre integration of the same interpolated profile, the quadrature's own node refinement being reported so its convergence is shown rather than asserted. The grid's cell-mass construction is then compared against those totals and **reported, not gated**: it spreads a thickened disk over spherical shells, so only its total is commensurable with a cylindrical integral, and it carries its own quadrature and outer-truncation error. Gating that comparison at 10⁻⁸ would be a tolerance no construction in this pipeline can meet — the same point the ruling makes about `_cylinder_mass` — so the measured difference is published instead, and any later step that needs the cells to be accurate must establish that separately. This replaces the earlier wording, which applied 10⁻⁸ to the cell comparison itself.

## Stage B: verify the field equations before any galaxy

**Completion I, one total potential matched to PM-1 in spherical symmetry.** Requiring μ(x)·x = y with x = y + √y, where x = g/a\* and y = g_N/a\*, gives

    μ_PM(x) = (√(1+4x) − 1)/(√(1+4x) + 1) = 4x/(√(1+4x) + 1)²,

the second form for small x. The field equation is ∇·[μ_PM(|∇Φ|/a\*)∇Φ] = 4πGρ_b with a = −∇Φ. The identity is verified over 10⁻¹⁰ ≤ y ≤ 10¹⁰ before use.

**Completion II, Newtonian gravity plus a separate nonlinear field.** ∇²Φ_N = 4πGρ_b together with ∇·[(|∇ψ|/a\*)∇ψ] = 4πGρ_b and a = −∇Φ_N − ∇ψ. In spherical symmetry g_ψ = √(a\*g_N), so the total matches Completion I; in a disk they need not agree, and that difference is the measurement.

**The existing simple-μ solver stays a separate comparison.** RPG-1's `aqual.py` implements μ = x/(1+x), whose spherical solution is ½[g_N + √(g_N² + 4a₀g_N)] — a different law: at g_N = a\* it gives 1.618a\* against PM-1's 2.000a\*. Running it unchanged and comparing with PM-1 would mix a change of force law with a change of geometry.

**Implementation requirement.** Changing `mu()` alone is insufficient: `nu()`, the boundary gradients, the iteration weights and the field-energy function must all be made consistent with whichever equation is being solved, and a test must fail if the three implementations are treated as the same equation.

**Verifications:** both completions reproduce g_N + √(a\*g_N) on a smooth spherical source to 10⁻⁶ relative; the simple-μ solver reproduces its own analytic solution to the same tolerance; the existing disk analytic controls still pass; and the finite-band ring calculation is checked against independently integrated source rings, never against the on-ring expression.

## Stage C: a matched disk comparison

Every model uses the same density grid, observation radii, thickness assumption, gas prescription and bulge prescription:

| Model | Purpose |
|---|---|
| Newtonian field of that density | the common ordinary-matter baseline |
| PM-1's local formula on that baseline | the empirical acceleration prescription alone |
| Completion I | one shared-field realization |
| Completion II | another |
| Existing simple-μ AQUAL | the established comparison, kept distinct |

Two declared passes. **First**, a common frozen a\*, to isolate the effect of the field equation. **Second**, each one-parameter candidate refits its own scale on the 89 training galaxies only, with the same objective and masks; validation and test are then reported without further adjustment. Model-score differences use paired galaxy resampling, not per-radius counting. The split is the archived one, so this is validation, not fresh confirmation.

**Not permitted:** comparing the new field scores directly against PM-1's archived scores and attributing the difference to geometry. RPG-1 already documents that its reconstructed density and the archived rotation contributions differ materially at some radii, so the baseline must be recomputed inside this experiment.

## Stage D: more than midplane rotation

- **The Milky Way vertical force** enters the first comparison, not a later one. RPG-1 found its own response overpredicts K_z at 1.1 kpc while improving rotation, which makes the vertical force a discriminator between the completions rather than an afterthought.
- **A two-concentration source**, with the field determined by the density and no declaration of whether it is one system or two. The current solver assumes axial and equatorial symmetry, so the control is a configuration compatible with those symmetries and is reported as such, not as a general two-galaxy calculation.

## Convergence is a gate

RPG-1's declared convergence gate failed, with 12 galaxies moving by more than 0.5 km/s under doubled resolution and 15 by more than 0.1 km/s when the outer domain was expanded. PM-2A keeps at least those gates, reports unconverged galaxies separately rather than letting an aggregate hide them, and **labels as unresolved any difference between completions smaller than the estimated numerical uncertainty**.

## What would make this informative either way

A completion that reproduces the spherical law, converges, and predicts the disk and vertical force better than the alternatives would be a viable shared-field realization — not a derivation of path memory. A completion that fails would be rejected as a realization, leaving memory itself untouched. Both outcomes are worth the run, and neither is a fit to be tuned.

## Files

`pm2a.py` runs the stages and writes `pm2a-results.json`; the field solvers extend RPG-1's machinery under this directory rather than modifying it in place, so RPG-1's own archive stays reproducible.

The stage A verifications are registered as their own suite job rather than added to PM-1's `checks.py`, so a PM-2A change cannot make PM-1's job fail and the two archives stay independently reproducible. This is where the job lives, not what it tests: `pm2a.py` regenerates its archive, compares every number, and exits nonzero if any gate fails or any number moves. Stage B's verifications join it when they exist.
