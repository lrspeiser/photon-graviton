# RUT-1 stage 4R: the two-stage response has an energy-accounted reciprocal form, and it corroborates the grid runs

[protocol-rut1.md](protocol-rut1.md), stage 4R declared in 8e2d681 before this run. This is the owner's proposed completion, not a change to the model under test: it asks whether the same two-stage response can be written so that the source coupling and the force come from one interaction term, and whether its energy then closes.

**All three checks pass.**

## The construction

With writing proportional to mass, q_i = α m_i, and ρ = Σ m_i δ(x − X_i), the two-stage equation is τ_form C_tt + γ C_t + C/τ_keep = α K∗ρ with γ = 1 + τ_form/τ_keep. The Gaussian kernel is a symmetric convolution square, W∗W = K, so writing C = W∗h gives

    τ_form h_tt + γ h_t + h/τ_keep = α W∗ρ,        matter acceleration = −∇Φ_N + ∇(W∗h),

and deposition (W∗ρ) and force readout (∇W∗h) are now one interaction. For a time-independent external potential,

    H = Σ m_i[|v_i|²/2 + Φ_N(X_i) − C(X_i)] + τ_form/(2α)∫h_t² + 1/(2ατ_keep)∫h²,        dH/dt = −(γ/α)∫h_t² ≤ 0.

The cross terms cancel identically — the matter energy changes by −Σ m_i ∂C/∂t(X_i), and the field energy by exactly that plus the damping — so the dissipated power has a **derived rate and a named destination**, the field's own damping, rather than being a variable introduced afterwards to absorb whatever went missing. The cancellation was checked in a finite-mode representation to 10⁻¹⁰ before declaration.

The implementation is spectral on a periodic box of side 5 with 64² modes, so deposition and readout use the same finite mode set and are exactly reciprocal in the truncated system; each mode's damped oscillator is advanced exactly for a source held at the step midpoint, and matter by the same velocity-Verlet step as the grid model. The two decay rates come out at exactly 1/τ_keep and 1/τ_form, the two stages of the cascade.

## A. It is the same response

Along prescribed circular trajectories for the 32-writer, w/R = 0.2 configuration, the spectral h-field and the grid two-stage field that every formation run uses are driven identically for 3 T0 and compared:

| where | C | ∇C |
|---|---|---|
| at the bodies | 1.2×10⁻⁵ | 6.6×10⁻⁵ |
| on probe rings at r = 0.9 and 1.1 | 7.4×10⁻⁶ | 4.3×10⁻⁵ |

Two independent spatial representations of one equation agree to a few parts in 10⁵, which bounds the grid representation error of the field in the formation runs.

## B. Its energy closes, at second order

For freely moving bodies over 3 T0, the balance H(t) − H(0) + ∫(γ/α)∫h_t² dt, relative to the bodies' kinetic energy:

| step | worst balance | dissipated |
|---|---|---|
| 0.02 | 5.33×10⁻⁶ | 0.3202 |
| 0.01 | 1.34×10⁻⁶ | 0.3202 |
| 0.005 | 3.35×10⁻⁷ | 0.3202 |

Halving ratios 3.977 and 3.996: the residual is the time-stepping error and nothing else, and the dissipated energy is converged. This is goal 4's first small coupled-system benchmark whose accounting closes.

## C. It corroborates the grid formation runs

From the same initial state and the same fixed step, free trajectories under the spectral field and under the grid field:

| time | max position difference | max velocity difference |
|---|---|---|
| 1 T0 | 4.8×10⁻⁷ | 4.7×10⁻⁷ |
| 2 T0 | 3.5×10⁻⁶ | 3.3×10⁻⁶ |
| 3 T0 | 9.2×10⁻⁶ | 7.7×10⁻⁶ |

An entirely different discretization of the same equation reproduces the trajectories the grid model computes. The difference grows with time, as it must between any two discretizations of a many-body system, so over the hundred-period horizon of stage 4 the comparison that matters is statistical rather than trajectory by trajectory.

## What this does and does not establish

**Shown.** The two-stage response admits a reciprocal form in which one interaction term both writes the field and exerts the force; its matter-plus-field energy closes against a derived dissipation rate to 3×10⁻⁷ and converges at second order; and it independently reproduces the grid field and the grid trajectories that the formation campaign relies on.

**Not shown.** The construction assumes writing proportional to mass and a symmetric Gaussian kernel — the source law of goal 3 remains a hypothesis. It does not identify the **physical reservoir** that receives the dissipated power; it names the channel, not the medium. It provides no **spatial causality** — the kernel is still instantaneous, and a formation delay is not a propagation speed. And it does not finish the **momentum accounting** for that reservoir and for the fixed central mass, which is external here. The reciprocal form is a completion of the effective model's energy, not evidence that such a field exists.

## Reproduce

```bash
python research_work/results/path-memory/rut4r.py
```
