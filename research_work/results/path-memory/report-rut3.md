# RUT-1 stage 3, corrected: the two-stage response wins every matched pair, and my ejections were numerical

[protocol-rut1.md](protocol-rut1.md); stage 3 declared in 390bc37, its numerical method in 868fe92, and the corrections and the two-stage candidate in 454155d — each before the run it governs. This supersedes [report-rut1-stage3.md](report-rut1-stage3.md), whose ejection claims did not survive the close-passage check.

## Two corrections to my own stage-3 report

**The ejections were numerical.** My `bounded` flag meant only that no coordinate crossed the domain edge, the timestep was fixed while the local dynamical time collapses inward, and the no-memory control was a *circular* orbit — which velocity-Verlet integrates almost exactly and which therefore cannot detect a mishandled close passage. On a bound eccentric Kepler orbit the same scheme conserves angular momentum to **1.2×10⁻¹³** while the energy is badly wrong. At the radii my plunges actually reached, fixed h = 0.01 gives a relative energy error of 4.3×10⁻⁹ at pericentre 0.11, **4.4×10⁻²** at 0.063, and at 0.03 throws a bound orbit to r = 247 and reports it unbound. The adaptive step cuts the worst of those to **6.7×10⁻⁴**, an improvement of 2.7×10⁴. The inward migration was real, because it is driven at radii where the integration is accurate. **The ejections and termination times were not.**

> **Correction 12, after the owner's review of 78b172a: those figures are endpoint errors.** Each was read after the integration finished, and a symplectic step returns close to its starting energy at the same orbital phase while carrying much larger error through pericentre. At pericentre 0.063 with fixed h = 0.002 the final error is 4.5×10⁻¹⁴ but the **maximum is 6.05×10⁻²**, so "converging away by h = 0.002" (454155d) was wrong. The gate now records the maximum energy error, the maximum position error against the analytic Kepler orbit, and the minimum radius. At the edge of the resolved domain, pericentre 0.25, the run stepping gives a **maximum** energy error of 9.4×10⁻⁵ and position error of 6.7×10⁻⁴; the two-stage runs never go below 0.25, so their results stand.

Three labels are now distinguished and never interchanged — *completed*, *left the domain*, *entered an unresolved central region* (r < 0.25).

> **Correction 13: the energy test overclaimed.** The quantity carried is a body's instantaneous *specific orbital energy* ε = |v|²/2 − GM/r − C, and in an evolving field it is not conserved: along the declared equations dε/dt = −∂C/∂t at the body. A positive value at one instant is therefore not a permanent-escape test, and a negative one is not a boundness guarantee. The fields are renamed accordingly. The same identity gives a stronger check, ε(t) − ε(0) + ∫∂C/∂t ds = 0, which stage 4 carries.

**"The obstacle is continued writing, not formation or field shape" was too strong**, and the split primed controls now say something sharper:

| control | 1 writer | 16 writers |
|---|---|---|
| frozen field | completed, L = **1.0000** | completed, L = **1.0000** |
| decaying, no new writing | completed, L = **1.0000** | completed, L = **1.0000** |
| decaying, with writing | **unresolved centre at 1.4 orbits**, L = 0.015 | completed, L = 0.873 |

Neither the field's presence nor its decay does anything. What destroys a lone writer is **the new uneven structure it adds behind itself**, and sixteen writers largely defuse it. That is a property of the *response process*, not of maintaining a field.

## The matched comparison: only the maturation changes

Same bodies, same total writing rate, same retention, same width, same initial positions and velocities, 2% phase and speed jitter, 20 orbits. τ_form = 3 periods.

| case | model | status | orbits | measured support | L retained | Δr | v_r rms | per-body L spread |
|---|---|---|---|---|---|---|---|---|
| 16, w/R 0.1 | one-stage | unresolved centre | 10.0 | −0.035 | 0.712 | −0.341 | 0.255 | 0.192 |
| 16, w/R 0.1 | **two-stage** | **completed** | 20.0 | **+0.021** | **0.996** | −0.016 | 0.086 | 0.050 |
| 32, w/R 0.2 | one-stage | completed | 20.0 | −0.040 | 0.837 | −0.261 | 0.139 | 0.286 |
| 32, w/R 0.2 | **two-stage** | completed | 20.0 | **+0.038** | **0.975** | −0.057 | 0.091 | 0.061 |
| 32, w/R 0.2, seed 2 | one-stage | completed | 20.0 | −0.017 | 0.864 | −0.203 | 0.154 | 0.273 |
| 32, w/R 0.2, seed 2 | **two-stage** | completed | 20.0 | **+0.037** | **0.994** | −0.029 | 0.059 | 0.035 |

> **Correction 14: the first two rows do not cover the same interval.** The one-stage sixteen-writer run reached the unresolved centre at 10.0 T0 and the two-stage run ran 20 T0, so their *completion* outcomes are directly comparable but their final-window support, migration, dispersion and spread are not measurements over identical intervals. Stage 4 records every model at common times.

The two-stage response wins every pair on every axis: it completes where the baseline does not, it retains 97.5–99.6% of angular momentum against 71–86%, it migrates by 1.6–5.7% against 20–34%, its radial-velocity dispersion is roughly half, its per-body angular-momentum spread is four to five times smaller — **and its measured support is positive where the baseline's is negative.** The baseline's negative support is the signature of bodies that have fallen inside their own trough; the two-stage rings have not.

Nothing about the equilibrium field was weakened to achieve this. The two models share the steady state C = τ_keep·S exactly, and the writing rate, retention and width are identical across each pair. The one-stage model is exactly τ_form = 0.

## Gates

| gate | requirement | measured |
|---|---|---|
| Close passage at the resolved edge (r = 0.25), **maximum** along the orbit | energy 10⁻³, position 5×10⁻³ | **9.4×10⁻⁵**, 6.7×10⁻⁴ (corrected from an endpoint test) |
| Transfer function, driving the shipped update with sinusoids | 2×10⁻³ | amplitude **4.0×10⁻⁶**, phase 1.3×10⁻⁶ |
| Two-state update against an independent ODE | 10⁻¹⁰ relative | **1.9×10⁻¹³** |
| Zero-frequency response unchanged | exact | 0.0 |
| Convergence **on a disturbed run** | 5×10⁻³ | **1.6×10⁻⁵** |
| No-memory control | Kepler exactly | L to 14 figures |

The transfer-function gate replaces an algebraic evaluation of H(ω) with sinusoids driven through the actual update; it realizes the owner's quadrature-ratio crossover — maturation attenuates fast variation but strengthens the phase-lagged response to variation slower than 2π√(τ_form τ_keep) = 34.4 T0, measured 0.08546, 0.68846, 0.99994 and 1.2554 at drive periods 5, 20, 34.41 and 100 T0 against the same values analytically. The convergence gate now refines the *disturbed* case that carries the claim — halved timestep and halved grid spacing — not the symmetric ring, which never approaches the centre and so could not have tested this. The ODE gate drives the shipped code path rather than a re-implementation of it, and it earned its place: **it caught a sign error in the two-stage update** (the exact step needs τ_k τ_f/(τ_f − τ_k), not (τ_k − τ_f)), which had made every two-stage run wrong until it was fixed.

## What this does and does not establish

**Shown.** The owner's two-stage candidate reproduces inside the repository, on matched initial conditions, with the baseline preserved unchanged. Letting a disturbance take time to become force-producing converts a migrating, heating, negatively-supported ring into one that holds its angular momentum and shows positive inward support. The earlier ejection claims are withdrawn as numerical. What destroys a lone writer is isolated to the fresh uneven structure it adds, not the field's presence or decay.

**Not shown.** These are 20-orbit runs, where a constant-source two-stage field has reached **80.72%** of mature (Correction 14: the 86% first reported here is the one-stage figure); the owner's exploratory 100-orbit runs found larger support (9.6–14.1%) and substantial radial motion, and that horizon has **not** been run here. The two-stage rings still migrate by a few percent and still heat relative to the no-memory control, so this is better finite-horizon survival, **not** demonstrated settling into steady collective orbits. Three configurations and two seeds is not a robustness map: differing radii, unequal masses and rates, vertical offsets, finite-width and overlapping annuli remain unrun. There is still no energy or momentum budget for the field, and a formation delay is **not** spatial causality. One planar ring around a fixed centre is not a galaxy.

The standing conclusion is the owner's: the tested immediate-response Gaussian develops severe migration and perturbed disruption; the two-stage excitation-to-force model substantially improves finite-horizon survival while retaining positive inward support; and that candidate still exhibits heating and secular migration and owes its physical accounting.

## Reproduce

```bash
python research_work/results/path-memory/rut3.py
```
