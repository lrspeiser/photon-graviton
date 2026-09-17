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

**Acceptance criteria.** Convergence is a gate: every reported run is repeated at half the timestep and half the grid spacing, and a conclusion that moves under refinement is reported as unresolved rather than as a result. The writing rate is set from the mature-ring prediction as a *label* only; the support actually achieved is **measured from the evolved field** and reported as such. No body's own footprint is excluded from its own force — the self-force is finite because the Gaussian's gradient vanishes at zero separation, and it is the drag mechanism under test. A body leaving the grid is flagged and its run reported as terminated, never silently truncated.

## The vector extension is behind a narrower dependency

Retained as a proposal, not made the first formation implementation. Its sideways term does no direct work, since v·[v × ∇×A] = 0, but the force law also carries −∂A/∂t, which cannot be dropped while the field grows, and its source law, sign and evolution are unspecified. Its action, source coupling and conservation analysis proceed separately; a full orbital simulation waits for a declared evolution law, so that it stays an independent physical proposal rather than an adjustable repair applied to whichever scalar run needs help.

## What would make this informative either way

A regime with useful inward support and an angular-momentum-change time long compared with the build time would be a genuine self-reinforcing track. A regime where drag always dominates would identify which ingredient needs revision — writing, retention, source normalization or torque — without waiting for the rest of the programme. Neither outcome is a rotation curve, and a track that supports one orbit is not a galaxy.

## Files

`rut1.py` runs stages 0–2 and writes `rut1-results.json`, registered as its own suite job. It shares no state with PM-2A or PM-3 and writes no file under another experiment's directory.
