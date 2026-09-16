# RUT-1: a lone writer faces a real support–torque trade-off — and several writers sharing one track largely escape it

[protocol-rut1.md](protocol-rut1.md), declared in 390bc37; the corrections and the collective control below declared in 0cb2069, both before the runs they govern. Stages 0–2 and 2C pass. Stage 3, releasing the body to write its own trajectory, is declared and not run.

**The headline.** A body on a prescribed circular orbit writes a genuinely attractive track whose inward support accumulates in proportion to the memory time, while the mature tangential drag approaches a **retention-independent limit**, (q/ΩR)·[1 − I₀e(R²/w²)]. Their mature ratio is

    drag / support = C(w/R) / N_build,     C_long(b) = b³/(2π)·[1 − I₀e(α)]/[I₀e(α) − I₁e(α)],  α = 1/b²,

with C → √(2/π) = 0.797885 for a narrow track. The cost of holding a prescribed orbit while the track builds is substantial: at 10% mature support and τ/T = 10, the driver must supply **9.81, 4.79, 2.27 and 0.98 times the body's own angular momentum** at w/R = 0.05, 0.1, 0.2 and 0.4, by which time the support has reached only 63.2% of mature.

**But a lone writer is not the only arrangement the model permits.** With N evenly spaced writers sharing one track at **fixed total** writing rate, the radial contributions add while the tangential ones cancel: sixteen writers keep **99.93%** of the inward support at **2.4%** of the per-body drag, and thirty-two at 0.05%. A single-writer drag bound therefore cannot be applied to a collective source.

**The defensible statement is narrower than a survival theorem:** for an isolated, continuously writing Gaussian source, substantial support comes with substantial formation torque, and raising the retention alone does not remove that trade-off in the long-memory regime.

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

The closed form Φ_ring(r) = −A·exp[−(r−R)²/(2w²)]·I₀e(rR/w²) matches the declared angular integral to **1.3×10⁻¹⁴** over five cases. The scaled Bessel function is the right implementation over the scanned range — I₀(1000) is inf in double while I₀e(1000) is finite — but my stated reason was wrong: I₀(100) = 1.07×10⁴² does **not** overflow.

**The mature-ring control reproduces the declared numbers** at GM = 1, R = 1, w = 0.1, D = 0.2: extra inward acceleration **10.03%** of Newtonian, supported speed **+4.89%**, κ² = **21.1500** against 1 without the track. This is motion in a fixed, already-built field, held still on purpose — a result to check, not evidence of formation.

**The source ring stays put while the probe moves.** Evaluating the same source at seven probe radii differs from the on-ring expression evaluated at the probe radius by up to a factor of 5 — PM-1's withdrawn row, re-tested rather than assumed retired. Off the plane the vertical acceleration is restoring at every height tested. The trough bottom sits **5.0×10⁻³ inside** the source path, which is why a body on its own track feels an inward force at all; a symmetric trough centred on R would give exactly zero there.

**A second, negligible-mass body feels the same track** — inward outside the trough bottom, outward inside it. The field is a function of position alone and nothing in the implementation keys on the writer, so this is a shared field and not a private trajectory constraint.

## Stage 2: the complete finite-memory force on a prescribed orbit

For a writer on a circle, the self-force separates exactly into one-dimensional integrals over its own past, with φ = Ω(t−s). Both integrands are 2π-periodic apart from the decay factor, so the entire history collapses to **one revolution times a geometric sum**, (1−ρ^N)/(1−ρ) with ρ = e^{−T/τ}. That reproduces the declared build law A(N) = A_∞[1−e^{−N/N_build}] identically, and it makes the integral exact rather than sampled.

That mattered. The radial integrand is **even** about each passage and a uniform grid handles it (converged to 2×10⁻¹¹ at the coarsest resolution tried); the tangential integrand is **odd** there, and uniform sampling got it wrong by 6×10⁻³ and converged only at second order. The tangential part is the drag — so a naive quadrature is least accurate on exactly the quantity the experiment is about.

| w/R | N_build | drag/support | C measured | C_long closed form |
|---|---|---|---|---|
| 0.05 | 100 | 0.156 | 0.781 | 0.781230 |
| 0.1 | 10 | 0.763 | 0.763 | 0.763125 |
| 0.2 | 10 | 0.361 | 0.722 | 0.722625 |
| 0.4 | 10 | 0.156 | 0.622 | 0.622683 |

**These rows are the MATURE field** — the infinite-past, finite-retention solution — not a build from zero. The finite-age force is a separate calculation, through the exact decomposition a(t) = a(T)(1−ρⁿ)/(1−ρ) + ρⁿ·a(u) for t = nT + u, applied separately to each component.

**Correction: the drag does not saturate after one revolution.** At *integer* revolutions both components carry the same geometric factor, so the drag follows the same build law as the support; within an orbit it does not, and it overshoots its mature value:

| elapsed | drag / mature | support / mature |
|---|---|---|
| 0.5 orbit | **1.03956** | 0.04984 |
| 1 orbit | 0.09516 | 0.09516 |
| 1.5 orbits | **1.03580** | 0.14026 |
| 10 orbits | 0.63212 | 0.63212 |

What is true is that the *mature* drag approaches a retention-independent limit as τ grows — asymptotically, not after one revolution. The formation has a within-orbit structure the mature map cannot display.

**An exact drag identity replaces the cancellation-prone route.** With E(u) = exp[−R²(1−cos Ωu)/w²] and H(t) = ∫₀ᵗ e^{−u/τ}E du, integrating the tangential term by parts gives drag(t) = (q/ΩR)·[1 − e^{−t/τ}E(t) − H(t)/τ], tending to (q/ΩR)[1 − I₀e(α)]. Both routes are computed and agree to **2.5×10⁻¹³**. H's tail uses the same exact periodic decomposition — truncating its revolution sum is what makes a naive check of the identity disagree, and that, not the shipped force, was the cause of a 3.6×10⁻³ discrepancy I first saw.

## The timescale crossing is not a survival theorem

The earlier Boolean tested L_time > τ. That is a chosen timescale diagnostic that **no freely moving orbit was subjected to**, and it is renamed `mature_drag_timescale_exceeds_one_retention_time`. An orbit may migrate substantially without being destroyed, and may lose an unacceptable fraction of its angular momentum while passing the inequality.

The coefficient is also a function of width, not a constant, so the crossing is solved from the force integrals rather than from a representative C. At τ/T = 10:

| mature support fraction | critical w/R |
|---|---|
| 1% of Newtonian | **0.04907** |
| 10% of Newtonian | **0.39343** |

The earlier single coefficient "4.4·f", which gave 0.044 and 0.44, came from collapsing C(w/R) to a median and is withdrawn.

**The accumulated cost is the sharper statement.** Holding the prescribed orbit to t = τ at 10% mature support requires driver angular momentum of 9.81, 4.79, 2.27 and 0.98 times the body's own at w/R = 0.05, 0.1, 0.2, 0.4 — substantial even for the wide case that passes the Boolean. It measures the effort needed to *prevent* the orbit from changing, not a loss a released body would suffer: once it moves, the prescribed-circle history no longer describes it.

**The PM-3 inference is withdrawn.** A broad track is not excluded by PM-3. Width does not change the storage equation, and Φ_steady = −τ·(steady writing source) is nonzero whether the source is narrow or broad. That a broad collective field might admit a static description raises a question about telling mechanisms apart observationally; it is not a proof that the stored field vanishes, and "not a narrow rut" describes geometry rather than grounds for discarding a response.

## Stage 2C: several writers sharing one track

Nothing in the model requires each body to take its support from its own private trail. With N evenly spaced writers on the same prescribed orbit at **fixed total** writing rate:

| writers | support relative to one | drag relative to one |
|---|---|---|
| 1 | 1.000000 | 1.000000 |
| 4 | 0.999315 | 0.218801 |
| 8 | 0.999292 | 0.088597 |
| 16 | 0.999290 | **0.023559** |
| 32 | 0.999291 | **0.000506** |

Sixteen writers retain more than 99.9% of the inward support at about a forty-second of the drag. The radial contributions from writers ahead, behind and around largely add; the tangential ones increasingly cancel.

**This is not the numerical-subdivision loophole**, and the two controls are kept distinct. Dividing *one* writer into coincident copies whose rates sum to the original leaves the field unchanged to **2.2×10⁻¹⁶** — a relabelling must change nothing, which is exactly what withdrew PM-1's per-ring saturation. Placing *distinct* writers at different physical positions changes the source distribution, so the field may legitimately change. That is the second case.

**Nor is it proof that the arrangement forms.** Evenly spaced writers held on their orbits are a deliberately favourable symmetry. Nothing here shows that they build the track from zero, keep their spacing, tolerate phase disturbances, survive differential motion, or obey a completed matter–field energy law.

**The uniform-ring control has no tangential force at all**, identically: ∫₀^{2π} sin φ·E(φ) dφ = 0 exactly, because that integrand is a total derivative. Using the averaged ring to argue the writing process is drag-free would be circular — the averaging is what removed the drag. That is the owner's warning, made exact.

## What this does and does not establish

**Shown.** The footprint model produces a genuinely attractive, radially and vertically confining track that a second body also feels. For a lone continuously writing source the mature drag/support ratio is C(w/R)/N_build with C in closed form, the accumulated formation torque is substantial at useful support fractions, and raising the retention alone does not remove that trade-off. And several physically distinct writers sharing one track at fixed total writing rate keep essentially all the support while cancelling most of the drag, so a single-writer bound does not transfer to a collective source.

**Not shown.** Stage 3 has not run: no body has yet been released to write its own trajectory from an empty field, so **support, stability and settling remain three separate unestablished claims**. Stages 0–2 are kinematic feasibility tests — the writing coefficient's dependence on source mass and the field's energy accounting are unspecified, and computing work and torque is informative but assigning a field-energy expression afterwards would not establish conservation. The instantaneous-Gaussian kernel is used throughout; a causal version, in which a disturbance written at X(s) cannot reach x before |x−X(s)|/c, is a separate obligation, since propagation speed and retention time are distinct properties. The vector extension is not implemented. And nothing here is a rotation curve: a track that supports one orbit is not a galaxy.

## Reproduce

```bash
python research_work/results/path-memory/rut1.py
```
