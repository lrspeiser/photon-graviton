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

## Stage 3: release the body and let it write its own trajectory

Declared here and run only after stage 2's map exists. Start from Φ_mem = 0, evolve trajectory and field together, with no circular path imposed and no target speed supplied. Three outcomes are distinguished and **not conflated**:

- **Support** — the written field supplies additional inward force.
- **Stability** — disturbances stay bounded.
- **Settling** — disturbances decay by transferring energy elsewhere.

A stable oscillating orbit is not a failure merely for not circularizing, and a circularizing orbit is not successful support if it loses most of its angular momentum and spirals in. The same declared parameters carry over from stage 2, and an already-built-track control runs beside the empty-field formation so a formation problem can be told from a mature-field problem.

## The vector extension is behind a narrower dependency

Retained as a proposal, not made the first formation implementation. Its sideways term does no direct work, since v·[v × ∇×A] = 0, but the force law also carries −∂A/∂t, which cannot be dropped while the field grows, and its source law, sign and evolution are unspecified. Its action, source coupling and conservation analysis proceed separately; a full orbital simulation waits for a declared evolution law, so that it stays an independent physical proposal rather than an adjustable repair applied to whichever scalar run needs help.

## What would make this informative either way

A regime with useful inward support and an angular-momentum-change time long compared with the build time would be a genuine self-reinforcing track. A regime where drag always dominates would identify which ingredient needs revision — writing, retention, source normalization or torque — without waiting for the rest of the programme. Neither outcome is a rotation curve, and a track that supports one orbit is not a galaxy.

## Files

`rut1.py` runs stages 0–2 and writes `rut1-results.json`, registered as its own suite job. It shares no state with PM-2A or PM-3 and writes no file under another experiment's directory.
