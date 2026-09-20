# RUT-1 stage 12: formation from an empty field, with the reciprocal energy and angular-momentum budget

Declared before any stage 12 run, and pushed with the code it names: `reciprocal_budget.py`, `rut12_tasks.py`,
`rut12.py`, `rut12_checks.py`. Stage 4R's reciprocal field, stage 7's population builder and the owner's
sampler are used **by import, unchanged**.

## Why this stage exists

This is the fourth and last of the pieces of work the owner set in the review of 3a80fec:

> **Test formation and physical accounting.** Evolve from an empty field toward a verified supported
> population, carrying the reciprocal energy and angular-momentum budget through the transition. *Whether the
> acceptable state can form naturally and persist without an unaccounted source or sink.*

Half of that accounting does not exist yet. Stage 4R built the reciprocal form C = W∗h, with W∗W = K, precisely
so that deposition and force come from one interaction term, and derived an exact energy ledger. Its own
docstring then says what it did not do: *"it does not finish the momentum accounting for that reservoir and the
fixed centre."* Without it, "carrying the angular-momentum budget through the transition" cannot be done at
all, only asserted. So this stage derives it first.

## The angular-momentum ledger

The external potential −GM/r is axisymmetric and exerts no torque, so all of the matter's torque comes from the
field. With acceleration +∇C and ρ = Σᵢ mᵢ δ(x − Xᵢ),

    dL_m/dt = ∫ ρ (x × ∇C)_z = ∫ ρ ∂_θ (W∗h).

W is radially symmetric, so ∂_θ commutes with W∗, and convolution by a symmetric kernel is self-adjoint:

    ∫ ρ ∂_θ (W∗h) = ∫ (W∗ρ) ∂_θ h = (1/α) ∫ (τ_f h_tt + γ h_t + h/τ_keep) ∂_θ h.

The last term is ½∫∂_θ(h²) = 0, and ∫h_tt ∂_θ h = d/dt ∫h_t ∂_θ h because ∫h_t ∂_θ h_t vanishes for the same
reason. So with

    **J = ∫ h_t ∂_θ h,    L_field = −(τ_f/α) J,    torque = (γ/τ_f) L_field**

the budget closes in the same shape as the energy one:

    **d(L_m + L_field)/dt = −(γ/τ_f) L_field.**

Unlike the dissipation, which is a square and cannot change sign, the torque is of either sign: a memory that
decays gives angular momentum back as readily as it takes it. That is a derived transport rate with a stated
destination, not a term fitted to absorb a discrepancy — the same standard stage 4R set for the energy.

Two things make it approximate in the implementation, and both are measured rather than hoped for. The box is
periodic, so ∂_θ means what it should only where the field has died away before the edge; and the square
lattice of modes is not rotation invariant, so the truncation can exert a torque of its own. **L_field is a
physical quantity and must therefore not depend on the box or the mode count** — that is a gate below, not an
aspiration.

## The construction

**The start is an equilibrium that knows nothing about the field.** An annulus with the same f(E, L) parameters
as the verified supported population but at α = 0: an exact equilibrium of the bare point mass, with no field
and no support. At t = 0 writing is switched on at the target's α, and nothing else changes. That is the
formation event — a single declared switch, with no prescribed path and no target state built into the start.

**The mature control is primed from the smooth density.** "The verified population with the field it writes"
means the self-consistent axisymmetric field, not a snapshot of point sources; it is set from the Hankel
transform of the population's own σ(r).

**The field is the reciprocal one**, spectral on a periodic box, because the ledgers are exact only for the
form in which deposition and readout come from one interaction term. The grid field every earlier formation run
used deposits K and reads ∇(K∗ρ), which is why stage 4R had to be written at all.

**What is measured is structure, not the mean radius.** The mass-weighted mean radius mixes the circular radius
with the epicyclic amplitude, so a hotter population sits at larger ⟨r⟩ for no interesting reason. Each body's
angular momentum, its circular radius in the field the bodies are actually in, and the epicyclic amplitude left
over are measured separately.

## The prediction, fixed before any run

While the mean field stays axisymmetric each body's angular momentum is conserved — the ledger above says by
how little it is not — and the field grows over τ_keep = 10 orbital periods, slow compared with one, so the
radial action is an adiabatic invariant. Then the L-distribution comes through formation intact, which fixes
the circular radii in the final potential; and the epicyclic amplitude follows √(κ_bare/κ_final).

Computed by `task_prediction` on the declared population at 4,096 bodies, before any formation run:

| | predicted to form | the target | ratio |
|---|---|---|---|
| ⟨L⟩ | 1.010915 | 0.998445 | 1.0125 |
| σ_L | 0.060665 | 0.059857 | 1.0135 |
| ⟨r_c⟩ | 0.941395 | 0.934921 | 1.0069 |
| σ_rc | 0.031397 | 0.031349 | **1.0015** |
| rms epicycle | 0.144323 | 0.091231 | **1.5820** |

So: **formation should reproduce the supported population's circular-radius structure to better than 1% and
leave it 58% hotter radially.** Note the ⟨L⟩ offset of 1.2%: the α = 0 and α > 0 distribution functions do not
share an L-distribution, because the phase-space measure depends on the potential. That is a property of the
construction, it is not removable, and it is why ⟨r_c⟩ is predicted 0.7% above the target's rather than equal
to it.

If the measured epicyclic amplitude comes out at the target's instead, the adiabatic argument is wrong and
something restores the target; if it comes out at the prediction, what forms is a supported annulus with the
right angular-momentum structure and more radial motion than the verified one. Neither outcome is selected in
advance, and neither is a gate.

## What was run before this declaration, and what it showed

Scratch prototypes on the configurations declared below. Sample sizes are given; nothing here is a result.

1. **The ledger, on the declared reference** (one configuration, three steps). Energy closes at 3.26×10⁻⁶,
   8.15×10⁻⁷, 2.03×10⁻⁷ and angular at 2.51×10⁻⁸, 6.30×10⁻⁹, 1.58×10⁻⁹ for steps 0.02, 0.01, 0.005 — both
   second order, ratios 4.01/4.00 and 3.98/3.98, against a total angular momentum of 1.0166.
2. **L_field is a property of the field** (two configurations). 1.874977×10⁻⁵ at box 5 with 64 modes and
   1.874977×10⁻⁵ at box 7 with 96 — identical to seven digits.
3. **The half-box control** (one configuration). `ifft2` returns the field on [0, L), so the signed coordinate a
   moment needs is that wrapped into [−L/2, L/2). Displacing it by half a box makes the angular closure
   5.16×10⁻⁴ instead of 6.30×10⁻⁹ — 82,000 times worse — and moves L_field by a factor 2.3. **The energy ledger
   does not notice at all**, because it is computed by Parseval and never touches a position. I made this
   mistake first and found it only because L_field varied with the box; it is declared as a control below.
4. **A formation run** (one realization, 256 bodies, 40 periods, 1,689 s). Support climbs from 0 to 0.1111
   against the target's 0.10257, the annulus contracting from ⟨r⟩ 1.094 to 1.017 and its spread from 0.2165 to
   0.130, while both ledgers stay closed at 3×10⁻⁷ throughout — through a transition in which the field goes
   from nothing to its full depth. The support was still rising slowly at the horizon, which is why a longer
   run is declared as a variant below.
5. **A reference configuration that failed, and why.** Priming a mature field from point sources puts every
   body at the bottom of its own Gaussian well; the field is then as lumpy as the sample, it relaxes violently,
   and the energy balance reaches 4.6 while the structure drifts by a factor 20. What fails there is the
   integrator, not the ledger. Two things follow, both declared below. The reference for part A is the
   formation configuration itself — an empty field. And the mature control is primed from the population's own
   **smooth** axisymmetric density, through the Hankel transform of σ(r), which is what "the field it writes"
   means; that field reproduces the one the population was built with to 1.2% at its worst and 0.46% rms over
   0.5 < r < 1.6, and gate A3 requires it.
6. **The step of the derivation the box can spoil.** ∫h ∂_θ h = ½∫∂_θ(h²) vanishes on a rotationally closed
   domain, and a periodic square is not one. On a localised field — eight sources on a ring, decayed long
   before the edge — it vanishes to 2.7×10⁻¹⁷, machine precision. On a field of random noise filling the whole
   box it does not vanish at all. That is the approximation the whole ledger rests on, it is checked directly
   in the suite job, and A2 gates it independently by requiring L_field to be the same in different boxes.

## Three statuses, kept separate (the owner's ruling)

| status | question | consequence of failure |
|---|---|---|
| reproduction | does the run reproduce its frozen archive? | the job exits non-zero |
| numerical verification | do both ledgers close and converge; is L_field a property of the field rather than of the box; do the runs complete and survive refinement? | blocks every dependent conclusion; the job exits non-zero |
| scientific outcome | what forms, against the prediction and against the verified population | archived as it fell; **never** affects the exit status |

No population is claimed stable, and whether what forms counts as "the acceptable state" is not a gate.

## Part A. The ledger

The reference is the formation configuration at a short horizon: the declared population's α = 0 equilibrium,
256 bodies, an empty field, two orbital periods.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| A1 both ledgers close | at the reference step the energy balance is under 10⁻⁵ and the angular balance under 10⁻⁷, and each falls by a factor in [3.5, 4.5] for each halving of the step | — (a convergence condition; A2 carries the comparison) |
| A2 L_field belongs to the field | across the declared boxes and mode counts, L_field agrees to 10⁻⁶ relative and the angular closure does not degrade | the position grid displaced by half a box must be seen to break the closure by more than 10⁻⁵ **and** move L_field by more than 10% |
| A3 the spectral field is the population's | primed from the population's own σ(r), the spectral field at the declared box and mode count reproduces the field the population was built with to 2% of its peak over 0.5 < r < 1.6 | — (without this neither the mature control nor the comparison means anything) |


The order is required only of pairs whose finer residual is still above 10⁻⁹: a convergence ratio measures
round-off rather than the leading error once the residual has reached it, and which pairs were judged is
recorded in the archive.

A2's control is the mistake I made: `ifft2` returns the field on [0, L), and a moment integral needs that
coordinate wrapped into [−L/2, L/2). The energy ledger cannot detect it, which is exactly why it is a control.

## Part G. Formation

The declared population's α = 0 equilibrium, 1,024 bodies, an empty field, writing switched on at t = 0, forty
orbital periods, four realizations; with two controls and five refinements on the first realization.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| G1 the run is sound | every declared run reaches its horizon, and at every recorded instant the energy balance is under 10⁻⁵ and the angular balance under 10⁻⁶ — through the transition, not only at its ends | — (a completeness and closure condition) |
| G2 the start is an equilibrium | with writing frozen, the α = 0 annulus must hold its structure to 2% over the horizon | — (this gate *is* the control that the formation seen in the other runs is caused by the writing) |
| G3 the target is a fixed point | started at the verified population with the field it writes, the state must hold its structure to 10% over the horizon | — (this gate *is* the control that the comparison below is against something the integrator can hold) |
| G4 refinement | halving the step, refining the modes, enlarging the box and doubling the bodies each move the settled structure by under 5% | — (a refinement condition) |

"Settled" means averaged over the records from period 30 to the horizon, declared before any run.

## The declared reading

For the four realizations, the settled structure — ⟨L⟩, σ_L, ⟨r_c⟩, σ_rc, rms epicyclic amplitude, support —
with its spread across realizations, set beside the prediction fixed above and beside the verified population.
The reading states which of the two the formed state matches in each quantity, and says plainly that a state
matching the prediction rather than the target is a supported annulus with the right angular-momentum structure
and more radial motion than the verified one. Whether that counts as the acceptable state forming naturally is
the owner's judgement, not a gate.

The energy dissipated and the angular momentum passed through the field over the transition are reported as
quantities, with their signs, since they are what "without an unaccounted source or sink" refers to.

## What this stage does not establish

The reservoir. The ledgers say where the energy and the angular momentum go — into the field, and out of the
field at a derived rate — but nothing here identifies what physically receives them, and a decaying memory is
still a postulate rather than a mechanism. Spatial causality: the field responds instantaneously in space.
The fixed centre: the central mass is held, so it absorbs linear momentum silently, and the linear-momentum
budget is not attempted. Anything beyond the declared horizon, which is forty periods and one run at eighty —
stage 11's m = 1 mode has an e-folding time of 181 periods, far longer than either, so this stage cannot see
it and does not speak to it. And one population: B24 alone is run, chosen because stage 11 showed its only
mode is that slow one.

## Declared thresholds (read by the driver; the code holds none of its own)

```json
{
    "R": {
        "population": "B24",
        "tau_form_periods": 3.0
    },
    "A": {
        "population": "B24",
        "start": "bare",
        "bodies": 256,
        "seed": 1,
        "orbits": 2.0,
        "samples": 8,
        "steps": [0.02, 0.01, 0.005],
        "reference_step": 0.01,
        "configurations": [
            [5.0, 64],
            [7.0, 96],
            [10.0, 128]
        ],
        "energy_closure": 1e-05,
        "angular_closure": 1e-07,
        "order_range": [3.5, 4.5],
        "field_invariance": 1e-06,
        "control_offset": 0.5,
        "control_closure_min": 1e-05,
        "control_field_min": 0.1,
        "order_min_residual": 1e-09,
        "representation_agreement": 0.02
    },
    "P": {
        "population": "B24",
        "bodies": 4096,
        "seed": 0
    },
    "G": {
        "population": "B24",
        "bodies": 1024,
        "step": 0.01,
        "modes": 64,
        "box": 5.0,
        "horizon": 40.0,
        "samples": 20,
        "seeds": [0, 1, 2, 3],
        "start": "bare",
        "controls": ["mature", "no_writing"],
        "variants": ["half_step", "fine_modes", "big_box", "more_bodies", "long_horizon"],
        "bodies_by_variant": {
            "more_bodies": 2048
        },
        "modes_by_variant": {
            "fine_modes": 96
        },
        "box_by_variant": {
            "big_box": 7.0
        },
        "horizon_by_variant": {
            "long_horizon": 80.0
        },
        "energy_closure": 1e-05,
        "angular_closure": 1e-06,
        "refinement_agreement": 0.05,
        "settled_from": 30.0,
        "no_writing_drift": 0.02,
        "mature_drift": 0.1
    },
    "F": {
        "ring_points": 64,
        "sample_radii": 120,
        "chunk": 512,
        "compare_r": [0.5, 1.6]
    },
    "D": {
        "structure_agreement": 0.05,
        "epicycle_agreement": 0.15
    },
    "C": {
        "replay_prefix": 2
    }
}
```
