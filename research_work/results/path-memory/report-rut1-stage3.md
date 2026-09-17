# RUT-1 stage 3: the collective helps enormously, and still does not survive a 2% perturbation

> **Superseded by [report-rut3.md](report-rut3.md).** Its ejection claims did not survive a close-passage check: the timestep was fixed while the local dynamical time collapses inward, and the no-memory control was a circular orbit, which cannot detect that. The inward migration reported below stands; the ejections, the `bounded` flags and the termination times do not. Its attribution of the trouble to "continued writing" is also narrowed there by split primed controls.

[protocol-rut1.md](protocol-rut1.md), stage 3 declared in 390bc37 with its numerical method in 868fe92, both before this run. Bodies and field are evolved together from Φ_mem = 0. No prescribed trajectory, no velocity resets, no target speed, and no writing rate adjusted to reach a chosen orbit. Support, stability, settling and formation cost are reported separately.

**The milestone is not reached.** The owner's next milestone is a *freely evolving, perturbed collective system* that builds its field from zero, gains measurable inward support, and stays bounded. The perfectly symmetric sixteen-writer ring does survive twenty orbits keeping 80% of its angular momentum. A **2% phase-and-speed perturbation breaks it up within thirteen orbits**, while the identical disturbed start with the memory switched off stays bounded and keeps its angular momentum to ten figures. The failure is attributable to the field, not to the disturbed initial orbits.

## What the runs show

All from an empty field, launched at the Newtonian circular speed, w/R = 0.1, τ/T = 10, and a writing rate **labelled** 10% mature support — the support reported is what the evolved field actually supplied.

| run | orbits | bounded | r: min / max / end | L retained | measured support, late |
|---|---|---|---|---|---|
| 1 writer | 5.5 | **no** | 0.063 / 2.475 / escaped | **0.075** | −0.015 |
| 4 writers | 20 | yes | 0.243 / 1.000 / 0.244 | 0.482 | −0.059 |
| 16 writers | 20 | yes | **0.889** / 1.000 / 0.899 | **0.802** | −0.249 |
| 16, jitter 0.02 | 12.9 | **no** | — | 0.544 | — |
| 16, jitter 0.06 | 10.7 | **no** | — | 0.565 | — |
| 16, jitter 0.02, **no memory** | 20 | yes | — | 1.0000 | — |
| 16, jitter 0.06, **no memory** | 20 | yes | — | 1.0000 | — |

The lone writer is destroyed almost exactly on the schedule stage 2's prescribed-orbit map predicted: at these parameters the mature drag gives an angular-momentum-change time of about 2.1 orbits, and the free body has lost 92.5% of its angular momentum by the time it leaves the domain at 5.5 orbits. **The prescribed-orbit calculation predicted the free outcome correctly**, which is the main reason to trust the map.

Collective cancellation then does most of what stage 2C suggested it would: sixteen writers sharing the same total writing rate turn a body destroyed in five orbits into a ring that has moved 10% inward after twenty.

## The problem is continued writing, not formation

The already-built-track control settles this. A lone body released into the *mature* ring field, launched in equilibrium with it — at the supported speed, not the Newtonian one, or the control would merely measure a body dropped into a deeper well — still collapses: r falls from 1.000 to 0.110 and it keeps 1.5% of its angular momentum.

So the difficulty is neither building the field nor the field's shape. Stage 1's mature-ring demonstration was stable because the field was held fixed **and the body was not writing**. Once a body writes continuously, it pays the drag whether or not the track already exists.

## Support, stability, settling and formation cost, kept apart

**Support.** The measured late-time inward memory acceleration is *negative* in every free run — the field pushes outward. That is not a sign error: a body that has migrated inward sits on the **inner flank of its own trough**, where the restoring force points outward. For the sixteen-writer ring that outward force is what arrests the contraction at r ≈ 0.90. The late-time state is therefore best described as *settled inside its own track*, not as a track supplying extra support. Calling it "10% support" would be reading the label rather than the measurement.

**Stability.** Bounded for 4 and 16 unperturbed writers over twenty orbits; unbounded for one writer and for every perturbed ring.

**Settling.** The sixteen-writer ring's radial oscillation decays (late/early amplitude 0.67); the four-writer ring damps harder (0.09) but only while spiralling to r = 0.24, which is migration, not settling. A circularizing run that loses three-quarters of its radius is not a success.

**Formation cost.** No driver acts in stage 3, so there is no external work: the whole cost appears as angular momentum and radius lost to the bodies' own field — 92.5%, 51.8% and 19.8% of L for 1, 4 and 16 writers.

## Numerical gates

| gate | requirement | measured |
|---|---|---|
| Grid against the stage-2 semi-analytic force | 2×10⁻³ | **5.5×10⁻⁴** inward, 3.5×10⁻⁶ backward |
| Deposition weight, h halved | field unchanged to 10⁻³ | **5.2×10⁻⁴** |
| Convergence, timestep halved | 5×10⁻³ | **5.3×10⁻⁵** |
| Convergence, grid spacing halved | 5×10⁻³ | **4.2×10⁻⁶** |
| No-memory control | Kepler exactly | radius 2.8×10⁻¹⁰, L to 10 figures |

The first is the cross-validation that licenses everything else: two entirely independent routes to the same force — a one-dimensional quadrature over the body's own past, and a two-dimensional grid carrying ∇Φ — agree to five parts in ten thousand. The second is the owner's requirement that halving the timestep must not double the writing rate; the deposit weight τ(1−e^{−h/τ}) → h makes the writing a rate rather than a per-step quantity, and that is tested rather than assumed.

## What this does and does not establish

**Shown.** A freely evolving body does build a field and does respond to it; the grid implementation reproduces an independent semi-analytic force; the lone writer's destruction matches the prescribed-orbit prediction quantitatively; collective writing at fixed total rate improves angular-momentum retention from 7.5% to 80%; and the obstacle is continued writing rather than formation or field shape. Runs that terminate early are reported as terminated, with their radius excursions, never silently truncated.

**Not shown, and the honest headline.** No perturbed collective configuration survived. Two percent is a small disturbance, and the memory-off control at the same disturbance is bounded with L = 1.0000, so this is the field breaking the ring. The declared failure is specific: **this Gaussian, continuously renewed, equal-writer planar ring at w/R = 0.1, τ/T = 10 and a 10% support label does not retain a perturbed collective configuration.** It is not a statement that collective gravitational memory cannot work — goal 2's sequence (differing radii, unequal rates, finite-width and overlapping annuli, several realizations) is exactly the space in which a viable region might still be found, and none of it has been run.

Also outstanding: one planar ring around a fixed centre, an instantaneous Gaussian kernel with no causal propagation, and no energy or momentum accounting for the field. Goals 3 and 4 remain untouched.

## Reproduce

```bash
python research_work/results/path-memory/rut1.py
```
