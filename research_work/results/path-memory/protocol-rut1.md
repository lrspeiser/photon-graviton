# RUT-1: can a moving body write an attractive track, and survive writing it?

Declared before execution, 16 September 2026. Baseline: `main` at d1dc7e0.

**The question.** PM-2A asks whether a specified shared-field equation predicts the field of an observed matter distribution. PM-3 showed that a *relaxational* memory reduces to a static law whenever the source is steady. RUT-1 asks the question neither of those answers: **can moving matter write a persistent field that supplies additional attraction — and survive the process of creating it?**

**Why it is a distinct experiment.** The footprint field retains a nonzero potential after its averaged amplitude stops changing. Writing at rate `q` with retention `τ`,

    ∂Φ_mem/∂t = −Φ_mem/τ − Σᵢ qᵢ·exp[−|x − Xᵢ(t)|²/(2w²)]

so a time-independent writing pattern leaves Φ_steady = −τ·(writing pattern), not zero. That is structurally different from a correction proportional to (p − Q), which vanishes when Q = p, and it is why PM-3's theorem does not dispose of this model. **But a nonzero stored field does not by itself distinguish rotation from an equivalent stationary source, and that is a control here, not an assumption.**

**What this is not.** RUT-1 is a low-speed effective model, not an equation derived from general relativity. The writing coefficient's dependence on source mass and the field's energy accounting are not specified, so the early stages are **kinematic feasibility tests** and are labelled as such: computing work and torque is informative, but assigning a field-energy expression afterwards would not establish conservation. RUT-1 also uses the instantaneous-Gaussian kernel throughout; a causal version, in which a disturbance written at X(s) cannot reach x before |x − X(s)|/c, is a separate obligation and is not implemented here. Propagation speed and retention time are distinct properties.

## Stage 0: the normalization check, first

Before any formation study. For the circular-history field the on-path potential depth is D = −Φ_mem(R), and the extra inward acceleration is ≈ D/(2R). **D is not automatically constant across orbits.** In the long-memory, narrow-track limit, with A = qτ,

    D ≈ A·w/(√(2π)·R)     ⇒     extra inward acceleration ≈ q·τ·w/(2√(2π)·R²).

At fixed writing rate, retention and physical width, a family of self-written circular tracks therefore approaches an **inverse-square** scaling, not inverse-radius. This is verified against the exact ring expression at several radii before anything else runs. It does not invalidate the rut's attraction or its stability; it means **"can create a supporting track" and "can explain flat rotation curves" are separate claims**, and holding D fixed independently at every radius and presenting the resulting 1/R force as a prediction is not permitted. A different width scaling, collective writing process or footprint shape may change this, and would be a declared mechanism to test.

## Stage 1: verify the kernel, then move the probe and not the source

1. **The closed form.** Φ_ring(r) = −A·exp[−(r−R)²/(2w²)]·I₀e(rR/w²) equals the declared angular integral to 10⁻¹²; the scaled Bessel form is used because I₀(rR/w²) overflows at the illustrative parameters.
2. **The mature-ring control**, reproduced as a *reported result to check*, not as new evidence of formation: at GM = 1, R = 1, w = 0.1, D = 0.2 the extra inward acceleration is ≈10.03% of Newtonian, the supported speed ≈4.89% above circular, and κ² ≈ 21.15 against 1 without the track.
3. **Probe scan with the source held fixed.** The same source ring is evaluated at many probe radii *and off its plane*. The source ring is never moved to the probe radius. This is the PM-1 on-ring error, and it is re-tested here rather than assumed retired.
4. **A second, negligible-mass test body** must feel the stored track. If only the writer responds, the implementation has produced a private trajectory constraint rather than a shared field, and that is a failure.
5. **The trough bottom** lies inside the source path, which is why a body on its own track feels an inward force; a symmetric trough centred on R would give exactly zero there. Verified, not asserted.

## Stage 2: prescribe the writer's orbit and measure the complete finite-memory force

A source on a prescribed circular orbit, its field building from zero. For X(s) on a circle of radius R at rate Ω, the self-force separates exactly, with φ = Ω(t−s):

    a_radial(t)    = −q ∫₀^t e^{−u/τ}·[R(1−cos φ)/w²]·exp[−R²(1−cos φ)/w²] du
    a_tangential(t)= −q ∫₀^t e^{−u/τ}·[R·sin φ/w²]·exp[−R²(1−cos φ)/w²] du

so no field grid is needed: both are one-dimensional quadratures over the body's own past. The radial integrand is even about each passage and accumulates; the tangential integrand is **odd** about each passage, so it nearly cancels within old revolutions but not across the most recent one, which is the drag the averaged ring removes by construction.

**Scanned:** w/R and τ/T_orbit. **Recorded, as a map and not a pass/fail label:** mean additional inward acceleration, mean tangential acceleration, the oscillating components, the time to build the field, and the driver work needed to hold the prescribed orbit. The recently written, uneven part of the trail is retained; replacing the actual history with a uniform ring is run alongside *as the control that shows what the averaging hides*.

The decisive diagnostic is

    angular-momentum-change time = |L / mean torque| = v/|a_tangential|

reported against both the orbital period and the field-build time: does the model build useful inward attraction before its own drag substantially changes the orbit?

### Corrections after the owner's review of 0c34b78

Five, all conceded, none of which changes a computed force on main — the shipped `self_force` agrees with an independent integration-by-parts route to 10⁻¹⁴, and its measured C converges exactly to the closed form below. What changes is what the map is called and what may be concluded from it.

**Correction 1: the published map is a *mature-field* calculation, not a build-from-zero one.** `stage2_prescribed_orbit` called `self_force` without `n_rev`, selecting S = 1/(1−e^{−T/τ}), the infinite-past finite-retention solution. Every row therefore describes a writer that has already been circling indefinitely. The rows are kept and relabelled, and the finite-age force is added, using the exact decomposition for t = nT + u:

    a(t) = a(T)·(1 − e^{−nT/τ})/(1 − e^{−T/τ}) + e^{−nT/τ}·a(u),

with a(T) and a(u) the fresh-history integrals over those intervals, applied separately to the radial and tangential components. Substituting a fractional revolution count into the geometric factor alone is **not** equivalent.

**Correction 2: "the drag saturates after one revolution" is wrong.** At *integer* revolutions both components carry the same geometric factor, so the drag follows the same build law as the support. Within an orbit it does not: at w/R = 0.1, τ/T = 10 the drag relative to its mature value is 1.03956 at half an orbit, 0.09516 at one, 1.03580 at one and a half, 0.63212 at ten. The correct statement is that the mature drag approaches a **retention-independent limit** as τ grows — asymptotically, not exactly, and not after one revolution.

**Correction 3: an exact drag identity replaces the cancellation-prone route.** With E(u) = exp[−R²(1−cos Ωu)/w²] and H(t) = ∫₀^t e^{−u/τ}E(u)du, integration by parts gives

    drag(t) = (q/ΩR)·[1 − e^{−t/τ}E(t) − H(t)/τ],     drag(∞) = (q/ΩR)·[1 − H(∞)/τ] → (q/ΩR)·[1 − I₀e(α)],

with α = R²/w², and the accumulated backward impulse J_drag(t) = (q/ΩR)·{t − ∫₀^t e^{−u/τ}E(u)[1 + (t−u)/τ]du}, whose product with R is the angular-momentum impulse the driver must supply. Both routes are computed and compared; H(∞) uses the same exact periodic decomposition, since truncating its revolution sum is what makes a naive check disagree.

**Correction 4: `w/R > 4.4·f` is not a survival theorem and C is not a constant.** The Boolean tested was `L_time > τ`, a chosen timescale diagnostic that no freely moving orbit was subjected to; it is renamed `mature_drag_timescale_exceeds_one_retention_time`. An orbit may migrate substantially without being destroyed, and may lose an unacceptable fraction of its angular momentum while passing the inequality. And C varies with width in closed form,

    C_long(b) = b³/(2π)·[1 − I₀e(α)]/[I₀e(α) − I₁e(α)],   b = w/R, α = 1/b²,   C → √(2/π) = 0.797885,  2πC → 5.013257,

so the threshold is solved from the actual force integrals rather than from a representative C. At τ/T = 10 that gives critical w/R = **0.04907** at 1% mature support and **0.39343** at 10% — both timescale crossings, not measured survival thresholds. The accumulated cost is the sharper statement: holding a prescribed orbit to t = τ at 10% mature support requires driver angular momentum of 4.79, 2.27 and 0.98 times the body's own, at w/R = 0.1, 0.2 and 0.4, by which time the support has reached only 63.2% of mature.

**Correction 5: PM-3 does not exclude a broad track, and the inference is withdrawn.** Width does not change the storage equation: Φ_steady = −τ·(steady writing source) is nonzero whether the source is narrow or broad. That a broad collective field might admit a static description raises a question about telling mechanisms apart observationally; it does not make the stored field vanish or bear on whether it forms. "Not a narrow rut" describes geometry, not grounds for discarding a response. (Separately: I₀(100) = 1.07×10⁴² does **not** overflow double precision. The scaled form is still the right implementation — I₀(1000) does overflow while I₀e(1000) is finite — but the stated rationale was wrong.)

### Stage 2C: a collective control, declared

Nothing in the model requires each body to take its support from its own private trail; the potential is a shared function of position, which stage 1 already verifies. With N evenly spaced writers on the same prescribed orbit at **fixed total** writing rate q_total (each writing q_total/N), the mature kernels become

    radial_kernel(u) = mean_j[(1 − cos φ_j)·E_j],   tangential_kernel(u) = mean_j[sin φ_j·E_j],   φ_j = Ωu + 2πj/N,

and the radial contributions largely add while the tangential ones increasingly cancel. Two controls must stay distinct and both are run: **dividing one writer into coincident copies** whose rates sum to the original must leave the field *exactly* unchanged — this is the numerical-subdivision loophole that killed PM-1's per-ring saturation, re-tested here — whereas **placing distinct writers at different physical positions** changes the source distribution, so the field may legitimately change. A single-writer bound therefore cannot be applied to a collective source without this test.

Evenly spaced, held-on-orbit writers are a deliberately favourable symmetry, so the collective arrangement is **not** thereby shown to form, keep its spacing, tolerate phase disturbances, survive differential motion, or obey a completed matter–field energy law. Stage 3 compares one writer with several physically distinct writers at the same total writing rate, both from an empty field, and perturbs the phases rather than testing perfect symmetry alone.

## Stage 3: release the body and let it write its own trajectory

Declared here and run only after stage 2's map exists. Start from Φ_mem = 0, evolve trajectory and field together, with no circular path imposed and no target speed supplied. Three outcomes are distinguished and **not conflated**:

- **Support** — the written field supplies additional inward force.
- **Stability** — disturbances stay bounded.
- **Settling** — disturbances decay by transferring energy elsewhere.

A stable oscillating orbit is not a failure merely for not circularizing, and a circularizing orbit is not successful support if it loses most of its angular momentum and spirals in. The same declared parameters carry over from stage 2, and an already-built-track control runs beside the empty-field formation so a formation problem can be told from a mature-field problem.

**The numerical method, declared.** The history integral is not evaluated by storing the trajectory — that costs the whole past at every step. Instead the field's *gradient* is carried on a fixed Cartesian grid, because differentiating the stored potential numerically would reintroduce exactly the readout error PM-2A stage B measured. From ∂Φ/∂t = −Φ/τ − q·Σᵢ exp[−|x−Xᵢ|²/(2w²)] it follows that

    ∂g/∂t = −g/τ + q·Σᵢ (x − Xᵢ)/w²·exp[−|x − Xᵢ|²/(2w²)],   g ≡ ∇Φ,   a_mem = −g,

so each component has an analytic source and the acceleration is read by interpolating g rather than by differencing Φ. Motion is planar and the sources lie in the plane, where the three-dimensional Gaussian factorizes exactly, so a two-dimensional grid is not an approximation; a vertical extension belongs to goal 2's robustness sequence. The decay is integrated exactly over a step and the source deposited at the midpoint,

    g(t+h) = e^{−h/τ}·g(t) + τ(1 − e^{−h/τ})·S(x_mid),

whose weight tends to h as h → 0, so **halving the timestep does not double the writing rate** — that invariance is a declared gate, not an assumption, and is tested directly by comparing the field built at h and at h/2.

### Corrections after the owner's review of ce5b0d7, and the two-stage candidate

**Correction 10: "the obstacle is continued writing, not formation or field shape" is withdrawn.** The primed control establishes something narrower — that *a lone writer, continuously writing under the same Gaussian rule, does badly even when it starts in a mature ring.* It does not separate field shape from writing, because it lets the old field decay while the body adds new uneven structure, and it was never run on a primed *disturbed collective*. The replacement statement is that **every fresh footprint immediately becomes an attractive disturbance behind the moving body while older structure persists elsewhere**, which is a property of the *response process* and not of maintaining a field as such. The primed controls are therefore split three ways — a frozen field, a decaying field with no new writing, and a decaying field with continued writing — and repeated for several writers.

**Correction 11: "escaped" is not "unbound", and the close passages were not resolved.** The `bounded` flag meant only that no coordinate crossed 0.9 of the domain half-width. Worse, the timestep was fixed while the local dynamical time collapses inward, and the no-memory control was a *circular* orbit, which velocity-Verlet integrates almost exactly and which therefore cannot detect a mishandled close passage: on a bound eccentric Kepler orbit the same scheme conserves angular momentum to 10⁻¹⁴ while the energy is badly wrong. Measured at the radii RUT-1's plunges actually reached, fixed-step h = 0.01 gives a relative energy error of 4.3×10⁻⁹ at pericentre 0.11, **4.4×10⁻²** at 0.063, and at 0.03 throws a bound orbit to r = 247 and reports it unbound — all converging away by h = 0.002. **The inward migration is credible, because it is driven at radii where the integration is accurate; the ejections and the precise termination times are not.** Three termination labels are now distinguished and never interchanged: *completed*, *left the computational domain*, and *entered an unresolved central region*. Unbinding is asserted only from a computed total energy including the memory potential, which requires carrying the scalar field and not only its gradient. The timestep becomes adaptive, set by the shortest local dynamical time, and the field is advanced consistently with it — which the exponential deposit weight τ(1−e^{−h/τ}) already permits, since it is a function of the step actually taken. Convergence is re-run on the **disturbed** cases that carry the breakup claim, not on the symmetric control.

**The two-stage candidate, declared as a separate model.** The owner's exploratory runs found that letting a disturbance take time to become force-producing performs substantially better. With S the writing pattern, E an excitation field and C the matured, force-producing depth,

    τ_form·∂E/∂t = S − E,      ∂C/∂t = E − C/τ_keep,      Φ_mem = −C,      a_mem = ∇C,

both starting from zero. The one-stage model is exactly τ_form = 0, and the steady state is unchanged: constant forcing gives E = S and C = τ_keep·S, so **the equilibrium field strength is not reduced by weakening gravity** — the comparison holds the total writing rate, retention time and footprint width fixed and changes only the maturation. Per harmonic, H_new(ω) = H_old(ω)/(1 + iωτ_form): the zero-frequency response is untouched while rapid variation is attenuated. That ratio is a statement about one harmonic, **not** a formula for the drag of a freely evolving collective, where transient torques of either sign must stay in the accounting. Equivalently τ_form·C̈ + (1 + τ_form/τ_keep)·Ċ + C/τ_keep = S, which is a starting point for a driven, damped field model whose energy and momentum budget is still owed; adding a formation delay is **not** spatial causality.

The declared comparison is matched: same bodies, same total writing rate, same retention time, same footprint width, same initial positions and velocities, **only the maturation process changes**. The one-stage Gaussian is preserved unchanged as the baseline. Reported alongside: per-body angular momentum distributions rather than only the mean, since a reasonable mean can hide one subset collapsing while another gains; and the radial-velocity dispersion, because a run that has not crossed a boundary but is heating is not a success.

**Acceptance criteria.** Convergence is a gate: every reported run is repeated at half the timestep and half the grid spacing, and a conclusion that moves under refinement is reported as unresolved rather than as a result. The writing rate is set from the mature-ring prediction as a *label* only; the support actually achieved is **measured from the evolved field** and reported as such. No body's own footprint is excluded from its own force — the self-force is finite because the Gaussian's gradient vanishes at zero separation, and it is the drag mechanism under test. A body leaving the grid is flagged and its run reported as terminated, never silently truncated.

## Stage 4: the instrumented hundred-period campaign, after the owner's review of 78b172a

The two-stage response stays the leading candidate. Stage 4 does not rewrite it; it asks whether the radial motion left in its runs is **initial adjustment, coherent oscillation, or a slow instability**, and it is declared here in full before any of its runs.

### Corrections 12–14

**Correction 12: the close-passage gate measured the endpoint, not the passage.** It read the energy error only after the integration finished, and a symplectic step returns close to its starting energy at the same orbital phase while carrying appreciable error through pericentre. At pericentre 0.063 with fixed h = 0.002 the final relative energy error is 4.5×10⁻¹⁴ and the **maximum 6.05×10⁻²**; so the statement in 454155d that those errors "converge away by h = 0.002" was an endpoint statement and is withdrawn. The gate now records the maximum energy error along the orbit, the maximum position error against the analytic Kepler solution, and the minimum radius, and it binds at the edge of the resolved domain, r = 0.25, where the formation runs stop: with the runs' η = 0.01 the maximum energy error must stay below 10⁻³ and the maximum position error below 5×10⁻³ over ten orbits (measured before declaration at 9.4×10⁻⁵ and 6.7×10⁻⁴). Deeper passages are reported, not gated. None of this invalidates the two-stage results, whose runs never go below 0.25.

**Correction 13: a positive specific orbital energy is not an escape test in an evolving field.** The quantity carried as "total energy" is a body's instantaneous specific orbital energy ε_i = |v_i|²/2 − GM/r_i − C(X_i,t). Along the declared equations dε_i/dt = −∂C/∂t at X_i exactly — for the two-stage model C/τ_keep − E — so a positive value at one instant does not establish permanent escape and a negative one does not guarantee boundness. The fields are renamed and the note corrected. The same identity supplies a much stronger check: R_i(t) = ε_i(t) − ε_i(0) + ∫₀ᵗ ∂C/∂t(X_i(s),s) ds must vanish, testing the moving body, the sampled force, the sampled potential and the field update together. It is not the matter–field energy budget.

**Correction 14: the two-stage field is less mature at 20 T0 than reported.** "86% of mature" was the one-stage constant-source figure. For the two-stage model, C(t)/C_∞ = 1 − [τ_keep e^{−t/τ_keep} − τ_form e^{−t/τ_form}]/(τ_keep − τ_form): **80.72%** at 20 T0, 99.04% at 50 T0, 99.9935% at 100 T0 — reference values for a constant source, not measurements of the moving system's field. Separately, rut3's sixteen-writer comparison ended at 10.0 T0 for the one-stage model and 20 T0 for the two-stage one: the completion outcomes are directly comparable, the final-window heating and migration numbers are **not** measurements over identical intervals, and stage 4 therefore records every model at common times. The transfer-function gate is also rebuilt to drive the shipped update with a sinusoid and fit amplitude and phase, rather than evaluating H(ω) algebraically; it realizes the quadrature-ratio crossover the owner derived, [1 + τ_form/τ_keep]/[1 + (ωτ_form)²] = 1 at a drive period of 2π√(τ_form τ_keep) = 34.4 T0, which suggests maturation suppresses the fast wake while modestly **strengthening** the phase-lagged response to slow variation — a hypothesis to test, not a proof of slow instability.

### The frozen equation and the declared runs

Frozen exactly as in rut3.py: the two-stage model with τ_keep = 10 T0 and τ_form = 3 T0, the one-stage baseline (τ_form = 0), the Gaussian footprint, the 10% writing label, 2% phase and speed jitter, adaptive stepping with η = 0.01 and h_max = 0.01, r_min = 0.25. Horizon **100 T0**. Twenty runs, all declared now:

| group | runs |
|---|---|
| Formation from an empty field | two-stage, one-stage and **no-memory** for each of 16 writers w/R 0.1 seed 1; 32 writers w/R 0.2 seeds 1 and 2 (nine runs) |
| Mature matched challenge | two-stage and one-stage, both initialized with the same supported ring field and the same perturbed bodies launched in equilibrium with it, the two-stage excitation set consistently (E = C/τ_keep); for the 16/0.1 and 32/0.2 seed 1 configurations (four runs) |
| Grid-rotation control | two-stage 32/0.2 seed 1 with the whole initial condition turned 0.3 rad on the Cartesian grid (one run) |
| Refinement over the full horizon | **every** two-stage formation configuration at half the timestep and step parameter, and at half the grid spacing (six runs) |

Refining all three configurations removes any choice of which case to refine after seeing the outcomes. The primed challenge is a diagnostic, not proof of formation: it separates "the two-stage response genuinely handles perturbations better" from "much of the benefit is that its force-producing field develops more slowly".

### What is measured during the integration

- **The ε-consistency residual** R_i above, accumulated every step.
- **The torque identity**: L_i(t) − L_i(0) = ∫(x × ∇C)_z dt = ∫∂C/∂θ dt, every step — exact for a central Newtonian force.
- **Mode-resolved memory power and torque on every body**, by angular band m ∈ {0, 1, 2, 3–4, 5–8, 9–16, 17–32, 33–64}: C sampled on five rings around each body's own radius, Fourier-decomposed in angle and reassembled at the body; specific torque ∂C_b/∂θ and power v·∇C_b, the radial derivative by a fourth-order stencil. These are accumulated at a cadence of 0.0025 T0, which calibration showed is needed to resolve the reversible epicyclic power exchange in the narrow-track case, and the bands must sum to the every-step totals. Bands extend past the writer-spacing scales of 16 and 32, where fresh fine structure lives. The axisymmetric band carries no torque but can still do work.
- **Radial mean flow and dispersion** in radial bins of 0.1: the mean radial velocity, and the radial and azimuthal dispersions about each bin's own mean — so a coherent contraction or breathing is not counted as heating, which the bare radial-velocity RMS reported in rut3 could not distinguish.
- **Angular-mode band powers** of the writing pattern S, the excitation E and the matured field C on fixed rings r ∈ {0.8, 0.9, 1.0, 1.1}.
- Scalars every 0.05 T0, spectra every 0.5 T0, full per-body states and band accumulations at 5, 10, 20, 50 and 100 T0, and the exact termination state separately from the last scheduled output.

### Analysis, declared

**Windows** early (0, 20], mid (20, 50] and late (50, 100] T0 — the late window beginning where the constant-source two-stage maturity passes 99%. Within each window, block means over 5 T0 and an ordinary least-squares trend on the block means with its standard error, which tames the autocorrelation of dense samples.

**Numerical uncertainty** for each two-stage configuration and each window quantity is the larger of its differences from the half-step and half-spacing runs, for both the window mean and the window change. If the dynamics are chaotic this includes trajectory divergence, which is the honest uncertainty of a single-trajectory statement.

**A window change is resolved** only if it exceeds three standard errors **and** three times its numerical uncertainty.

**Classification of the late window**, with thresholds fixed now:

| verdict | declared test |
|---|---|
| candidate heating or instability | any of: secular migration (mean radius change resolved and above 0.01); angular-momentum drift (resolved and above 0.01); dispersion growth (radial dispersion change resolved, positive, above 10% of its window mean and above 0.001); excursion growth (radius spread change resolved, positive, above 10% and above 0.005) |
| bounded coherent oscillation | none of the above, and the standard deviation of the mean radial flow exceeds the mean radial dispersion |
| formation adjustment, then stationary | none of the above, and the late-window mean radius differs by more than 0.01 from the **matched no-memory control's** late-window mean |
| stationary | none of the above |

The adjustment test is referenced to the no-memory control rather than to the starting radius because a smoke test of this pipeline, run at a compressed horizon before declaration, classified a *no-memory* ring as adjusted: every jittered body starts at r = 1 on a slightly eccentric orbit, so the ring's instantaneous mean radius breathes coherently at the orbital frequency, and comparing two instants measures epicyclic phase. Positive support is claimed only if the late-window mean support is positive and exceeds three times its numerical uncertainty. Several flags may hold at once and all are reported; the verdict is the first row that applies. **A contraction that stops is not a failure, and neither is a bounded breathing mode** — the target is to eliminate continuing secular migration and growing random motion, not every radius change.

### Gates

| gate | requirement |
|---|---|
| ε-consistency | worst relative residual over all memory runs and bodies below 10⁻³, and the half-step run's residual at least twice smaller than its base run's |
| torque identity | worst \|ΔL − ∫∂C/∂θ dt\| below 10⁻¹⁰ |
| band attribution | at the last checkpoint, bands reproduce the every-step torque to 2% and work to 5% of the summed per-body magnitudes (work is cadence-limited) |
| grid rotation | rotated and unrotated late-window means of radius, dispersion, support and angular momentum agree within the larger of three times the numerical uncertainty and the seed-1-to-seed-2 spread |

### What follows the diagnosis

Only what the diagnosis identifies changes next: bounded breathing → assess it and test slower preparation or broader initial orbital distributions before adding damping; fine-scale angular structure carrying torque or heating → more distributed writers, finite-width annuli and a declared maturation-time scan (τ_form/T0 = 1, 3, 6 at fixed writing rate and retention, run long enough that success cannot mean delayed onset); a growing slow global mode → compare its frequency with the field response, since longer maturation may not help; a work-balance residual that grows under refinement → fix deposition, readout or integration before interpreting physics; a stable field with decaying support → storage and source normalization, not orbital damping. **No third force term, hand-damped radial velocity or imposed symmetry** is added before that diagnosis.

### Stage 4R, in parallel: a reciprocal, energy-accounted form of the same response

A proposed completion from the owner's review, not a change to the model under test. With writing q_i = α m_i and ρ = Σ m_i δ(x − X_i), the two-stage equation is τ_form C_tt + γ C_t + C/τ_keep = α K∗ρ with γ = 1 + τ_form/τ_keep. Writing the Gaussian kernel as a symmetric convolution square, W∗W = K, and C = W∗h, the same dynamics follow from τ_form h_tt + γ h_t + h/τ_keep = α W∗ρ with matter acceleration −∇Φ_N + ∇(W∗h), so source coupling and force come from one interaction term. For a time-independent external potential,

    H = Σ m_i[|v_i|²/2 + Φ_N(X_i) − C(X_i)] + τ_form/(2α)∫h_t² + 1/(2ατ_keep)∫h²,   dH/dt = −(γ/α)∫h_t² ≤ 0,

a derived dissipation rate rather than a variable added afterwards to absorb missing energy; the cancellation was checked in a finite-mode representation to 10⁻¹⁰ before declaration. The benchmark: a spectral periodic box in which deposition and force readout use the same finite mode set, first reproducing the grid two-stage C at prescribed body trajectories, then free motion with the discrete balance H(t) − H(0) + ∫(γ/α)∫h_t² converging under refinement. Its limits are stated now: it assumes writing proportional to mass and the symmetric kernel, it does not identify the physical reservoir, it provides no spatial causality, and it does not finish the momentum accounting for the reservoir and the fixed centre.

### What stage 4 is not

One planar ring around a fixed centre, an instantaneous Gaussian kernel, a writing rate set from a label, and no field energy budget for the grid model. The hashes of formation.py, longrun.py, rut4.py, rut3.py and rut1.py are recorded at launch, so the committed code can be verified to be the code that ran. A verdict classifies **this frozen equation under these conditions** and is not a verdict on gravitational memory in general. `rut4.py` is the science driver; `rut4_checks.py` is the suite job, rerunning a deterministic 5 T0 prefix against the committed series and recomputing the analysis from them.

## The vector extension is behind a narrower dependency

Retained as a proposal, not made the first formation implementation. Its sideways term does no direct work, since v·[v × ∇×A] = 0, but the force law also carries −∂A/∂t, which cannot be dropped while the field grows, and its source law, sign and evolution are unspecified. Its action, source coupling and conservation analysis proceed separately; a full orbital simulation waits for a declared evolution law, so that it stays an independent physical proposal rather than an adjustable repair applied to whichever scalar run needs help.

## What would make this informative either way

A regime with useful inward support and an angular-momentum-change time long compared with the build time would be a genuine self-reinforcing track. A regime where drag always dominates would identify which ingredient needs revision — writing, retention, source normalization or torque — without waiting for the rest of the programme. Neither outcome is a rotation curve, and a track that supports one orbit is not a galaxy.

## Files

`rut1.py` runs stages 0–2 and writes `rut1-results.json`, registered as its own suite job. It shares no state with PM-2A or PM-3 and writes no file under another experiment's directory.
