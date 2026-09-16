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

### Two source-geometry corrections, after the owner's review of e9e931d

Neither changes a fitted parameter or the outer continuation. Both concern integrating the same declared source over the correct volume, and both are fixed in PM-2A's own implementation with RPG-1's archive left intact.

**Correction 8: `m_cyl` added a spherical bulge mass to a cylindrical integral.** The disk terms integrate a surface density over a cylinder; the bulge term added `m_bulge(R)`, which is an enclosed mass inside a *sphere*. Their sum is neither. For a spherical component with cumulative mass M(r), the mass inside an infinite cylinder of radius R is

    M_cyl(R) = M_sph(R) + ∫_R^∞ M′(r)·[1 − √(1 − R²/r²)] dr,

the extra term being the portions of exterior shells that project inside the cylinder (a shell at r > R contributes the fraction 1 − √(1 − R²/r²)). The bracket is evaluated in the cancellation-safe form (R²/r²)/(1 + √(1 − R²/r²)), and the integral under the substitution u = √(r² − R²), which removes the square-root singularity of the integrand's derivative at r = R.

*The analytic test.* For a Plummer sphere M(r) = M r³/(r²+b²)^{3/2}, the projected surface density gives M_cyl(R) = M R²/(R²+b²) in closed form, so M_cyl/M_sph = √(R²+b²)/R. At R = b that is √2: substituting the spherical mass **underestimates the cylindrical mass by 29.2893%**, reproducing the owner's figure. The shell integral above returns the closed form. This test is added before the component-wise agreement test is extended, because the existing integral-agreement gate checks only the stellar interior and cannot see this.

*What it changes here.* A cylinder includes material outside the sphere whenever the source extends past the evaluation radius, so in spherical symmetry M_cyl ≥ M_sph rather than M_cyl = M_sph. The `four_masses` note asserting that the force-equivalent and cylindrical masses coincide in spherical symmetry is therefore wrong and is replaced: it is M_sph that M_force equals in spherical symmetry, by Newton's theorems. **Measured impact on the archived numbers: none at the radii reported.** The adopted bulge's tabulated support ends exactly at the outermost rotmod radius in all 31 bulge galaxies, so at R_last the cylinder already contains the whole bulge and M_cyl = M_sph = M_bulge. At half the bulge support the median change is 0.000%, the 90th percentile 0.005% and the maximum 0.120% (UGC11914), because `np.maximum.accumulate` saturates the adopted bulge mass well inside those radii. The correction is made because the expression is wrong, not because it moved a number; a different bulge prescription, or an inner-radius evaluation, would not be so forgiving.

**Correction 9: the thick-disk builder assigns a thin-disk mass to its inner sphere.** `exponential_disk_cell_masses` builds cells for ρ(R,z) = Σ(R)e^(−|z|/h)/(2h) but returns Σ(0)·πr_in² as the unresolved core, which is the mass inside a *cylinder* of radius r_in at constant Σ — the razor-thin disk's answer. The solver consumes it as the mass inside the inner *sphere*. Doing the z integral first at each cylindrical radius gives the exact spherical inner mass for a varying Σ,

    M_inner = 2π·∫₀^{r_in} R·Σ(R)·[1 − e^(−√(r_in² − R²)/h)] dR  →  2πΣ(0)r_in³/(3h)  for r_in ≪ h,

so the code exceeds it by a leading ratio 3h/(2r_in). At the grid's own scales — h = rd and r_in = 0.01·rd — the exact ratio is **150.5631**, reproducing the owner's figure; at h = 0.1·rd it is 15.57. That is a large relative error in an unresolved core whose absolute contribution can still be small, but it feeds the inner boundary flux and the solver's total mass, so it is **not** to be folded into the published cell-total discrepancy as though that were all quadrature and outer truncation. Stage B's spherical tests use `spherical_cell_masses`, which differences the supplied spherical cumulative mass and is unaffected. Before any disk comparison the thick disk's inner mass is replaced by the spherical integral and its convergence checked as r_in varies.

**Two gate-hardening changes in the same pass.** The row-removal gate recomputes a fingerprint of the live source specification before and after prediction instead of comparing a cached string with itself, and that fingerprint covers the bulge's full radial profile rather than only its total, plus the vertical prescription once one exists. And `gates_passed` additionally requires that the construction-failure list is empty and that the expected number of sources was built, so a future run cannot pass on a reduced sample.

## Stage B: verify the field equations before any galaxy

Let x = |∇Φ|/a\* for whichever potential the equation solves, and y = g_N/a\*. Every equation below is ∇·[μ(|∇Φ|/a\*)∇Φ] = 4πGρ_b for its own μ.

| Equation | μ(x) | spherical solution | ν(y) = x/y | relaxation weight w = 1/(1+η), η = dlnμ/dlnx |
|---|---|---|---|---|
| Newtonian | 1 | x = y | 1 | 1 |
| **Completion I**, total potential | 4x/[1 + √(1+4x)]² | x = y + √y | 1 + 1/√y | √(1+4x)/[1 + √(1+4x)], from η = 1/√(1+4x) |
| **Completion II**, auxiliary ψ only | x | x = √y | 1/√y | 1/2 |
| Existing simple-μ comparison | x/(1+x) | x = ½(y + √(y²+4y)) | ½ + √(¼ + 1/y) | (1+x)/(2+x) |

The μ_PM identity μ_PM(x)·x = y at x = y + √y is verified over 10⁻¹⁰ ≤ y ≤ 10¹⁰ before use. The weights are each the locally matched spherical weight that cancels the first-order iteration error, derived rather than inherited; `aqual.py`'s hardcoded (1+x)/(2+x) is the simple-μ member of that family and is correct only for simple μ. They are numerical iteration weights: not a proof of convergence for a nonspherical solve, and never a physical memory time.

**Completion II is two solves, added as vectors.** ∇²Φ_N = 4πGρ_b and ∇·[(|∇ψ|/a\*)∇ψ] = 4πGρ_b are solved separately and the accelerations added as vectors, a = −∇Φ_N − ∇ψ, not by summing their magnitudes. The auxiliary μ is never evaluated on the total gradient, and the total boundary gradient is never supplied to the auxiliary solve; either mistake changes the equation while possibly leaving a plausible-looking rotation curve. In spherical symmetry g_ψ = √(a\*g_N), so the total matches Completion I; in a disk they need not agree, and that difference is the measurement.

**Boundary gradients** use the selected equation's own spherical inverse, computed directly from the enclosed mass through y = GM/(r²a\*) rather than as y·ν(y), since ν diverges as y → 0 while x does not.

**The field functionals are derived, not chosen.** The gradient term in the action and the equation are not independent: with u = |∇U|²/a\*² and x = √u, the requirement is dF/du = μ(√u), and the density is a\*²F(u)/(8πG). This is `aqual.py`'s existing normalization. These follow, in forms free of cancellation:

| Equation | F(u), F(0) = 0 | small-x behaviour |
|---|---|---|
| Newtonian | u | exact |
| Completion I | t³(t + 2/3), with t = 2x/[1 + √(1+4x)] so that x = t² + t | exact, no subtraction; → (2/3)x³ |
| Completion II | F_N(u_N) + (2/3)u_ψ^{3/2}, evaluated on each field separately | exact |
| simple-μ | u − 2√u + 2·ln(1 + √u) | series (2/3)x³ − ½x⁴ + (2/5)x⁵ − ⅓x⁶ + … below x = 0.1 |

Here t = √y on the spherical solution, which is why Completion I's form has no subtraction of nearly equal numbers: at x = 10⁻¹² it returns 6.667×10⁻³⁷ = (2/3)x³ exactly. `aqual.py`'s `field_energy_function(y) = y − 2√y + 2·log1p(√y)` is the simple-μ F and *is* cancellation-prone despite `log1p`, because its leading terms cancel to O(x³): at x = 10⁻⁸ it returns exactly 0.0 against a true 6.667×10⁻²⁵, and at x = 10⁻⁶ it is wrong by 7.7×10⁻⁵ relative. PM-2A's implementation takes the series branch below x = 0.1, checked against a high-precision reference at the transition; RPG-1's archived energy results are not silently changed.

**The functional densities and their cross-check.** Completion I's is a\*²F_PM(|∇Φ|²/a\*²)/(8πG); Completion II's is |∇Φ_N|²/(8πG) + |∇ψ|³/(12πGa\*), the two fields entered separately. On the matched spherical solution the two must agree identically,

    E_I = E_II = [a\*²/(8πG)]·[y² + (2/3)y^{3/2}],

which follows from t = √y. This tests the energy implementation independently of checking the final acceleration. These are field-functional contributions, not the complete conserved matter–field energy: calling them an energy reservoir would require the source coupling, gauge convention and boundary terms as well, and none of them explains how matter powers a persistent disturbance.

**What counts as passing.** The source is the same Plummer profile used for the bulge test, M(<r) = M r³/(r²+b²)^{3/2} with g_N = GMr/(r²+b²)^{3/2}, its dimensionless compactness GM/(a\*b²) varied so the weak, transition and strong regimes are all exercised. The 10⁻⁶ relative tolerance applies to identified numerical outputs, not to a mixture of boundary conditions, scalar identities and solutions:

1. **Interior field and flux.** At interior faces, excluding the faces where the answer was imposed as a boundary condition, the solved field must satisfy that equation's own spherical flux identity r²·μ(g/a\*)·g = G·M(<r), with angular symmetry checked too. The auxiliary solve is checked on g_ψ, Completion I on the total solved field.
2. **Off-grid readout.** Interpolated gradients are tested at radii that are not grid faces. Passing a face-flux test does not establish that the readout used for galaxy observations meets the same accuracy — `midplane_speed` interpolates face quantities in log r while `gradient_at` interpolates cell centres, so they are different numerical paths.
3. **Actual convergence.** An acceptable equation residual *and* a small iteration change, both bound to explicit criteria. `aqual.py` computes a residual but sets `converged` from the iteration change alone; both measurements are preserved and both are gated.
4. **Equation distinction.** At g_N = a\*, Completions I and II must give x = 2.000000 and the simple-μ solver 1.618034, and a test must fail if the three are treated as one equation. The two completions **must agree** on the spherical benchmark — requiring them to differ there would itself be a wrong test.
5. **Gauge.** These isolated models have a logarithmic far field, so Φ(∞) = 0 is not a valid normalization; the finite-radius outer pin is kept and every energy comparison uses the same convention.

A static solve that passes all of this establishes an equation implementation, not path memory, and the solver's relaxation time is a numerical quantity that must never be read as a physical memory time.

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
