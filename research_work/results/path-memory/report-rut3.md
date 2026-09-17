# RUT-1 stage 3, corrected: the two-stage response wins every matched pair, and my ejections were numerical

[protocol-rut1.md](protocol-rut1.md); stage 3 declared in 390bc37, its numerical method in 868fe92, and the corrections and the two-stage candidate in 454155d — each before the run it governs. This supersedes [report-rut1-stage3.md](report-rut1-stage3.md), whose ejection claims did not survive the close-passage check.

## Two corrections to my own stage-3 report

**The ejections were numerical.** My `bounded` flag meant only that no coordinate crossed the domain edge, the timestep was fixed while the local dynamical time collapses inward, and the no-memory control was a *circular* orbit — which velocity-Verlet integrates almost exactly and which therefore cannot detect a mishandled close passage. On a bound eccentric Kepler orbit the same scheme conserves angular momentum to **1.2×10⁻¹³** while the energy is badly wrong. At the radii my plunges actually reached, fixed h = 0.01 gives a relative energy error of 4.3×10⁻⁹ at pericentre 0.11, **4.4×10⁻²** at 0.063, and at 0.03 throws a bound orbit to r = 247 and reports it unbound. The adaptive step cuts the worst of those to **6.7×10⁻⁴**, an improvement of 2.7×10⁴. The inward migration was real, because it is driven at radii where the integration is accurate. **The ejections and termination times were not.**

Three labels are now distinguished and never interchanged — *completed*, *left the domain*, *entered an unresolved central region* (r < 0.25) — and unbinding is asserted only from a total energy that includes the memory potential, which is why the scalar field is carried and not only its gradient.

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

The two-stage response wins every pair on every axis: it completes where the baseline does not, it retains 97.5–99.6% of angular momentum against 71–86%, it migrates by 1.6–5.7% against 20–34%, its radial-velocity dispersion is roughly half, its per-body angular-momentum spread is four to five times smaller — **and its measured support is positive where the baseline's is negative.** The baseline's negative support is the signature of bodies that have fallen inside their own trough; the two-stage rings have not.

Nothing about the equilibrium field was weakened to achieve this. The two models share the steady state C = τ_keep·S exactly, and the writing rate, retention and width are identical across each pair. The one-stage model is exactly τ_form = 0.

## Gates

| gate | requirement | measured |
|---|---|---|
| Close passage, adaptive step | 10⁻³, and ≥10³ better than the coarsest fixed step | **6.7×10⁻⁴**, 2.7×10⁴ better |
| Two-state update against an independent ODE | 10⁻¹⁰ relative | **1.9×10⁻¹³** |
| Zero-frequency response unchanged | exact | 0.0 |
| Convergence **on a disturbed run** | 5×10⁻³ | **1.6×10⁻⁵** |
| No-memory control | Kepler exactly | L to 14 figures |

The convergence gate now refines the *disturbed* case that carries the claim — halved timestep and halved grid spacing — not the symmetric ring, which never approaches the centre and so could not have tested this. The ODE gate drives the shipped code path rather than a re-implementation of it, and it earned its place: **it caught a sign error in the two-stage update** (the exact step needs τ_k τ_f/(τ_f − τ_k), not (τ_k − τ_f)), which had made every two-stage run wrong until it was fixed.

## What this does and does not establish

**Shown.** The owner's two-stage candidate reproduces inside the repository, on matched initial conditions, with the baseline preserved unchanged. Letting a disturbance take time to become force-producing converts a migrating, heating, negatively-supported ring into one that holds its angular momentum and shows positive inward support. The earlier ejection claims are withdrawn as numerical. What destroys a lone writer is isolated to the fresh uneven structure it adds, not the field's presence or decay.

**Not shown.** These are 20-orbit runs, where the field has reached 86% of mature; the owner's exploratory 100-orbit runs found larger support (9.6–14.1%) and substantial radial motion, and that horizon has **not** been run here. The two-stage rings still migrate by a few percent and still heat relative to the no-memory control, so this is better finite-horizon survival, **not** demonstrated settling into steady collective orbits. Three configurations and two seeds is not a robustness map: differing radii, unequal masses and rates, vertical offsets, finite-width and overlapping annuli remain unrun. There is still no energy or momentum budget for the field, and a formation delay is **not** spatial causality. One planar ring around a fixed centre is not a galaxy.

The standing conclusion is the owner's: the tested immediate-response Gaussian develops severe migration and perturbed disruption; the two-stage excitation-to-force model substantially improves finite-horizon survival while retaining positive inward support; and that candidate still exhibits heating and secular migration and owes its physical accounting.

## Reproduce

```bash
python research_work/results/path-memory/rut3.py
```
