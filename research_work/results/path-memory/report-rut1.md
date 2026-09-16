# RUT-1: the track's support grows with memory, its drag does not — and survival buys width, not time

[protocol-rut1.md](protocol-rut1.md), declared in 390bc37 before this run. Stages 0–2 pass. Stage 3, releasing the body to write its own trajectory, is declared and not run.

**The headline.** A body on a prescribed circular orbit does write an attractive track, and the track's inward support accumulates in proportion to the memory time. But the tangential drag it feels comes from the passage it has just made, so it **saturates after one revolution and is independent of the memory time entirely**. Their ratio is

    drag / support = C / (N_build · w/R),     C = 0.59 – 0.78,

with N_build = τ/T_orbit. Long memory and wide tracks both reduce the drag relative to the support. But the question that decides survival is not the ratio — it is whether the orbit outlives the build, and there the memory time **cancels**:

    L-change time / build time = (w/R) / (2π·C·f)     ⇒     the orbit survives only if w/R > 2π·C·f ≈ 4.4·f,

where f is the extra inward support as a fraction of Newtonian. **Buying 10% support needs a "track" about 44% of the orbital radius.** That is an annulus comparable to the orbit, not a rut — and at that width it is hard to tell from a smooth axisymmetric mass, which puts it back inside PM-3's steady-state theorem. Buying 1% needs only 4%, which *is* a genuine track. Neither the writing rate nor the memory time changes this; they scale support and drag together.

## Stage 0: the normalization check, run first

At fixed writing rate, retention and physical width, a family of self-written circular tracks is **inverse-square**, not inverse-radius:

| source-orbit radius | on-orbit inward acceleration | ×R | ×R² |
|---|---|---|---|
| 1 | 0.0200226 | 0.0200226 | 0.0200226 |
| 2 | 0.0049915 | 0.0099829 | 0.0199659 |
| 4 | 0.0012470 | 0.0049879 | 0.0199518 |
| 8 | 0.0003117 | 0.0024935 | 0.0199483 |

a·R² varies by 0.37% across a factor of eight in radius while a·R varies by 187%, converging on the predicted q·τ·w/(2√(2π)R²) = 0.0199471. This reproduces the owner's three-row table exactly. It does not invalidate the rut's attraction or stability; it means **"can create a supporting track" and "can explain flat rotation curves" are separate claims**, and holding D fixed independently at every radius and presenting the resulting 1/R force as a prediction is not permitted.

## Stage 1: the kernel, and moving the probe rather than the source

The closed form Φ_ring(r) = −A·exp[−(r−R)²/(2w²)]·I₀e(rR/w²) matches the declared angular integral to **1.3×10⁻¹⁴** over five cases. The scaled Bessel function is not a nicety: I₀(rR/w²) = I₀(100) overflows at the illustrative parameters.

**The mature-ring control reproduces the declared numbers** at GM = 1, R = 1, w = 0.1, D = 0.2: extra inward acceleration **10.03%** of Newtonian, supported speed **+4.89%**, κ² = **21.1500** against 1 without the track. This is motion in a fixed, already-built field, held still on purpose — a result to check, not evidence of formation.

**The source ring stays put while the probe moves.** Evaluating the same source at seven probe radii differs from the on-ring expression evaluated at the probe radius by up to a factor of 5 — PM-1's withdrawn row, re-tested rather than assumed retired. Off the plane the vertical acceleration is restoring at every height tested. The trough bottom sits **5.0×10⁻³ inside** the source path, which is why a body on its own track feels an inward force at all; a symmetric trough centred on R would give exactly zero there.

**A second, negligible-mass body feels the same track** — inward outside the trough bottom, outward inside it. The field is a function of position alone and nothing in the implementation keys on the writer, so this is a shared field and not a private trajectory constraint.

## Stage 2: the complete finite-memory force on a prescribed orbit

For a writer on a circle, the self-force separates exactly into one-dimensional integrals over its own past, with φ = Ω(t−s). Both integrands are 2π-periodic apart from the decay factor, so the entire history collapses to **one revolution times a geometric sum**, (1−ρ^N)/(1−ρ) with ρ = e^{−T/τ}. That reproduces the declared build law A(N) = A_∞[1−e^{−N/N_build}] identically, and it makes the integral exact rather than sampled.

That mattered. The radial integrand is **even** about each passage and a uniform grid handles it (converged to 2×10⁻¹¹ at the coarsest resolution tried); the tangential integrand is **odd** there, and uniform sampling got it wrong by 6×10⁻³ and converged only at second order. The tangential part is the drag — so a naive quadrature is least accurate on exactly the quantity the experiment is about.

| w/R | N_build | drag/support | C = ratio·N_build·(w/R) | L/build at f = 0.1 | at f = 0.01 |
|---|---|---|---|---|---|
| 0.05 | 3 | 5.163 | 0.775 | 0.103 | 1.028 |
| 0.05 | 100 | 0.156 | 0.781 | 0.102 | 1.019 |
| 0.1 | 10 | 0.763 | 0.763 | 0.209 | 2.087 |
| 0.2 | 10 | 0.361 | 0.722 | 0.441 | 4.408 |
| 0.4 | 10 | 0.156 | 0.622 | 1.023 | 10.229 |
| 0.4 | 100 | 0.0156 | 0.623 | 1.022 | 10.224 |

**Why the memory time cancels.** The support is ∝ N_build; the drag is not, because it comes from the most recent passage. So the angular-momentum-change time at fixed support fraction is ∝ N_build — and so is the build time. Longer memory buys a better instantaneous ratio and exactly as much extra time to build, leaving survival unchanged.

**The uniform-ring control has no tangential force at all**, identically: ∫₀^{2π} sin φ·E(φ) dφ = 0 exactly, because that integrand is a total derivative. Using the averaged ring to argue the writing process is drag-free would be circular — the averaging is what removed the drag. That is the owner's warning, made exact.

## What this does and does not establish

**Shown.** The footprint model produces a genuinely attractive, radially and vertically confining track that a second body also feels. Its support accumulates with memory while its drag saturates, so drag/support falls as 1/(N_build·w/R). And the survival question has a clean answer independent of writing rate and memory time: **w/R > 4.4·f**, which is comfortable at 1% support and untenable at 10%, where the "track" would have to be almost half the orbital radius.

**Not shown.** Stage 3 has not run: no body has yet been released to write its own trajectory from an empty field, so **support, stability and settling remain three separate unestablished claims**. Stages 0–2 are kinematic feasibility tests — the writing coefficient's dependence on source mass and the field's energy accounting are unspecified, and computing work and torque is informative but assigning a field-energy expression afterwards would not establish conservation. The instantaneous-Gaussian kernel is used throughout; a causal version, in which a disturbance written at X(s) cannot reach x before |x−X(s)|/c, is a separate obligation, since propagation speed and retention time are distinct properties. The vector extension is not implemented. And nothing here is a rotation curve: a track that supports one orbit is not a galaxy.

## Reproduce

```bash
python research_work/results/path-memory/rut1.py
```
