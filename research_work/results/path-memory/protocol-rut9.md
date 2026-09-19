# RUT-1 stage 9: the corrected response calculation

Declared before any stage 9 run. Its own protocol file. Stages 4 to 8 stay frozen: their archives keep
reproducing, stage 6's H6 stays failed, **stage 8's L3 and L3b stay failed in stage 8's archive and stage 8's
checks job keeps exiting non-zero**, and no file of theirs is edited — including `warm_modes.py`, whose
omission this stage repairs in a new library beside it.

## Why this stage exists

The owner's ruling on 3a80fec. Stage 8 has two separate unfinished issues, and repairing the first will not
by itself explain the second:

1. **the ring-limit verification is not passed.** A third verification was approved, as a numerical repair and
   not a threshold negotiation: *computed annulus root − exact ring root = physical finite-width difference +
   numerical representation error*, and the two must be separated. Each narrow annulus is to be resolved on
   its own — kernel spacing, radial coverage, retained basis rank, a shifted placement of the nodes, the
   orbit quadrature and the population's own discretization — with an explicit error budget comfortably
   under the ring-limit tolerance, demonstrated on the annuli and not borrowed from `A_cold` or B16. The
   amended 0.5% stays the physical target, the frequency difference is reported separately, the fitted order
   is a diagnostic, and stage 8's under-resolved three-node construction becomes a negative control. The same
   checks are extended to B13, B24 and the near-threshold branch, and stage 8's forecasts are recalculated
   **beside** the originals, never over them;
2. **B13 grows 22% more slowly than predicted, at 5.4 standard errors, for an unknown reason.** Of the four
   diagnostics the owner asked for, this stage runs the three that need no new simulation of the live model:
   the boundary-derivative audit of the sharp truncation, an eigenvalue sensitivity, and a time-domain check
   that imposes **B13's own eigenfunction at its predicted complex rate** — stage 8's L2 tested one kernel
   column in each of two populations, and no mode. The fourth, the seeded full-state experiment, is stage 10
   and has its own protocol. No outcome is selected in advance.

## The correction

The population that is built, drawn and simulated is exactly zero outside declared bounds:

    f0(E, L) = g(E, L) Θ(E_max(L) − E) Θ(L − L_lo) Θ(L_hi − L),
    g = A exp[−(L−L0)²/2dL²] exp[−(E−E_c(L))²/2dE²],   E_max = E_c(L) + 4 dE,   L_lo,hi = L0 ∓ 4 dL.

`warm_modes.py` differentiates g and omits what the bounds contribute. With dE_max/dL = dE_c/dL = Ω_c(L),

    ∂f0/∂E     = g_E Θ… − g δ(E − E_max) …
    ∂f0/∂L|_E  = g_L Θ… + g Ω_c δ(E − E_max) … + g [δ(L − L_lo) − δ(L − L_hi)] …

so F_l = ν_l ∂f0/∂E + m ∂f0/∂L|_E gains three families of boundary rows, each a one-dimensional integral:

| family | F_l on it | orbits | measure |
|---|---|---|---|
| energy edge | (m Ω_c − ν_l) g(E_max(L), L) | (L, E_max(L)) | (2π)² dL / Ω_r |
| lower L edge | + m g(E, L_lo) | (L_lo, E) | (2π)² dE / Ω_r |
| upper L edge | − m g(E, L_hi) | (L_hi, E) | (2π)² dE / Ω_r |

The lower energy bound E = E_c(L) is not a truncation: it is J_r = 0, which no perturbation carries a body
across, and g_E vanishes there. The cutoff is **not** replaced by a taper: a taper appears below only as the
independent check of these rows. `warm_response.py` carries them; on stage 8's quadrature and without them
it is stage 8's calculation, which gate B1 verifies to the last bit of the orbit arrays.

## What was run before this declaration, and what it showed

A scratch prototype, not archived. It is stated because the campaign will reproduce these numbers, and because
every number below was obtained **on the configuration it is declared for** — the last of them, the refined
grid of the narrowest annulus, took three hours, and the declaration waited for it:

* **without the boundary rows the new library is stage 8's**: orbit arrays identical, |T − T_stage8| ≤
  1.1×10⁻¹⁶, and the four archived roots recovered;
* **the boundary rows are small for the populations and are not the cause of B13's shortfall**: they move
  the m = 2 growth rate by +0.004% (`A_cold`), +0.005% (`B_cold`), **+0.006% (B13)** and +0.015% (B16); the
  root moves by 1 to 4 parts in 10⁴ of |s|, almost all of it in the frequency and almost all from the energy
  edge. On the narrow ring-limit annuli they matter more: 0.08% to 0.10% of the growth rate, 9×10⁻⁴ to 10⁻³ of
  |s|, ten times the convergence standard — they belong in the calculation;
* **they are the limit of a resolved taper**: with every bound replaced by a quintic smoothstep 0.08, 0.04 and
  0.02 widths wide, on dedicated 24-node panels, and with ∂f/∂E and ∂f/∂L taken as central differences of the
  function so that none of the algebra above enters, the tapered response minus the sharp interior tends to
  the boundary rows linearly (1.9%, 0.91%, 0.44% of their norm for `A_cold`; 1.1%, 0.56%, 0.28% for B13) and
  the Richardson limit agrees with them to 2×10⁻⁵ (`A_cold`) and 7×10⁻⁵ (B13); each family alone to 10⁻³ or
  better; without the rows the same limit misses by 100% of them;
* **the repaired node rule** — kernel centres every 0.05 = w/4 from 0.3 = 1.5 w inside the smallest pericentre
  to 0.3 outside the largest apocentre — gives the narrow annuli 13 to 18 nodes where stage 8's gave 3 to 5.
  **On each of the four annuli, separately:** every kernel, rank, quadrature, angle and harmonic variant moves the
  root by under 3×10⁻⁷ of |s|; the annulus's own grid is the largest term, 3.4×10⁻⁶, 1.1×10⁻⁶, 1.8×10⁻⁶ and
  1.9×10⁻⁶; the sums are 3.4×10⁻⁶, 1.2×10⁻⁶, 1.8×10⁻⁶ and 2.2×10⁻⁶, against a budget of 5×10⁻⁴; between the nodes
  the mode equation holds to 10⁻⁹ or better. Stage 8's rule, put through the same kernel refinements, moves its
  root by 1.6×10⁻⁴, 6.0×10⁻³, 4.1×10⁻³ and 7.5×10⁻³, and between its nodes the mode equation is violated by
  7×10⁻⁶, 2.2×10⁻⁴, 2.7×10⁻³ and 3.3×10⁻³. **With the boundary rows and the new rule the growth rate differs
  from the exact-history ring by +23.45%, +4.79%, +1.13% and +0.278%**, the frequency by −6.3×10⁻³, −6.8×10⁻⁴,
  −1.4×10⁻⁴ and −3.4×10⁻⁵; the observed orders between successive annuli are 2.15, 2.05 and 2.02. Without
  the boundary rows the same sequence reads +23.37%, +4.70%, +1.03%, +0.18%: it heads for −0.10% and not for
  zero, because the rows' share of the response does not shrink with the width. That, and the three-node
  basis, are what bent the sequence L3 fitted an order of 2.5 to;
* **B13 under the same refinements**: every kernel, quadrature, angle and harmonic variant moves its root by
  under 2×10⁻⁹ of |s|; **the population's own grid is the largest term, 2.1×10⁻⁵** (stage 7's `all_together`
  refinement); the sum is 2.1×10⁻⁵. Recalculated root 0.03216180 − 0.08786641i against stage 8's forecast
  0.03215979 − 0.08785351i: the growth rate moves by +0.006%. A deliberately coarse rule — four nodes 0.3
  apart — moves it by 2.1%, which is the control;
* **the mode equation between the nodes** separates an adequate basis from an inadequate one as sharply as
  refinement does: on the second annulus the regenerated field differs from the imposed one by 10⁻¹⁰ between
  the nodes under the new rule and by 2.2×10⁻⁴ under stage 8's, whose root is off by 2.5×10⁻⁴;
* **sensitivity**: B13's mode is carried entirely by the l = −1 radial harmonic, 2Ω_θ − Ω_r, far from any
  resonance (|ν − i s| ≈ 1.1); +1% of the whole response moves the growth rate by +0.95% and a first-order
  formula reproduces an actual re-solve to 0.4%. **A uniform change of the response of −23% would be needed to
  lower the rate by 22%.** The same formula with the memory transfer function held fixed in dT/ds is off by
  111%, which is the control;
* **the root-free sweep** of B24's m = 2 strip under the new representation returns winding number zero in
  all eleven rectangles, and the same sweep finds B_dE0.018's root at 0.010586 − 0.104756i (stage 8: 0.0106);
* **B13's own eigenfunction imposed at its predicted rate**: six batches of 200,000 bodies at ε = 0.002 implied
  a shift of the growth rate of −0.63 ± 0.19% and of the frequency of (−0.8 ± 0.6)×10⁻⁴, every batch negative.
  On ONE draw of 50,000 bodies the implied shift is −0.180%, −0.251%, −0.536% and −1.672% at ε = 0.0005,
  0.001, 0.002 and 0.004: **exactly quadratic in ε, −0.095% per (10⁻³)²** — the offset was the finite amplitude
  of my own test and not the response, so the declared ε is 0.0005, where it is −0.02%. Halving the step moves
  the implied shift by 6×10⁻⁶ of the rate. Lengthening the switch-on from 8 to 12 e-foldings is not a
  refinement of the same kind — the same bodies reach t = 0 in other orbital phases, which is another sample of
  the response — and it moved the result by 0.5%, the size of a 50,000-body batch's noise; the residual
  transient of an eight e-folding switch-on is e⁻⁸ = 3×10⁻⁴ of the mode. A batch of 200,000 takes 69 minutes;
* **the other populations of N3 and N4**: `A_cold`, `B_cold`, B16 and family B at 0.017 and 0.018 each move by a
  sum of 1.9×10⁻⁵ to 2.6×10⁻⁵, all of it the population's own grid, with the mode equation holding to 6×10⁻¹¹
  between the nodes; the coarse rule is caught on each (0.5% to 5.7%; at 0.018 it loses the root altogether).
  Family B at 0.019, the root-free population nearest the branch, returns winding number zero in every
  rectangle of its m = 2 strip under all four declared configurations, and so do its m = 1, 3 and 4 strips
  and all four strips of `B_warm`, the warmest, under the base representation; the control's root at
  dE = 0.018, 0.010586 − 0.104756i, is found under every configuration. An m = 4 strip takes 35 minutes and
  the quadrature refinement twice as long: part N4 is about forty core-hours;
* **seen before this declaration, and belonging to stage 10, not to this stage**: a prototype of the seeded
  experiment. A four-fold symmetric sample of B13 keeps |C₂|/C₀ at 10⁻¹⁶ for twenty periods and does not heat,
  and seeded with the predicted eigenmode at 10⁻⁶ it grows as a clean exponential over five e-foldings: at
  0.0329, 0.0309 and 0.0316 in three sampling realizations at 1,024 bodies, pattern speeds +0.0454, +0.0459
  and +0.0452, against a predicted 0.0322 and +0.0439. A prototype, not a measurement. It bears on none of
  this stage's gates; it is written here because I have seen it.

## Three statuses, kept separate (the owner's ruling)

| status | question | consequence of failure |
|---|---|---|
| reproduction | does the run reproduce its frozen archive? | the job exits non-zero |
| numerical verification | do the boundary rows, the representation, the ring limit, the sensitivity formula and the time-domain check pass, and reject their negative controls? | blocks every dependent conclusion; the job exits non-zero |
| scientific outcome | what the corrected calculation says about stage 8's forecasts and about B13's shortfall | archived as it falls; never affects the exit status |

**This stage's status is its own.** Whatever it finds, stage 8's archive, report and checks job are not
edited, and L3 and L3b are failed there.

## The declared representation

Base: Gauss-Legendre nodes in L (192) and in u = √[2(E − E_c)] (64) over the declared four-width support; 128
midpoint nodes in the orbit's eccentric-type angle; radial harmonics l = −8 to 8; **kernel centres every 0.05
from 0.3 inside the smallest pericentre to 0.3 outside the largest apocentre**, clipped to [0.35, 2.9]; T(s)
posed in the orthonormal basis of the kernel matrix's subspace above 10⁻¹⁰ of its largest eigenvalue; **the
three families of boundary rows on the same L and u nodes.** Roots are located and counted as in stage 8.

Nine refinement variants, each changing one thing: kernel spacing 0.025; kernel margin 0.5; every centre
shifted by half a spacing; rank cutoff 10⁻¹²; rank cutoff 10⁻⁸; orbit quadrature 288 × 96; 256 angle nodes;
harmonics to ±12; and the population's own discretization refined — for a ring-limit annulus the radial grid
doubled with 240 × 144 nodes in (L, v_r), for a population stage 7's `all_together` refinement (559 radial
nodes on [0.2, 3.3], 320 × 192). A population of family A is *defined* by its mean support, so on a refined
grid its coupling is solved again; one of family B is defined by the family's coupling, which is held.

## Part B. The boundary rows

| gate | requirement | negative control the gate must reject |
|---|---|---|
| B1 reproduction | for `A_cold`, `B_cold`, B13 and B16, on stage 8's quadrature and stage 8's nodes and WITHOUT the boundary rows: the orbit arrays equal stage 8's exactly, T(s) at the archived root agrees to 10⁻¹³, and the archived root is recovered to 10⁻⁹ | the same comparison WITH the boundary rows must be seen to differ, by more than 10⁻⁹ in T |
| B2 the rows are the limit of a resolved taper | for `A_cold` and B13, at the population's stage 8 root, in the orthonormal basis: the tapered response (above) minus the sharp interior, extrapolated to zero taper width, equals the boundary rows to 1% of their norm — for all three bounds together and for each alone — and the unextrapolated difference falls at every halving | the same limit against the interior alone must miss by more than 50% of the rows' norm |

The shift each family of rows makes in each root is a reading, reported with its first-order estimate.

## Part N. The numerical representation

| gate | requirement | negative control the gate must reject |
|---|---|---|
| N1 every annulus resolved on its own | for each of the four annuli of stage 8's L3b, (dL, dE) = (0.04, 0.004), (0.02, 0.001), (0.01, 0.00025), (0.005, 0.0000625), each solved for a 9.99% mean support: exactly one root in the declared rectangle with the winding number agreeing; **each of the nine variants moves it by under 10⁻⁴ of \|s\|; the nine movements sum to under 5×10⁻⁴, a tenth of the ring-limit tolerance;** and between the nodes the regenerated field differs from the imposed one by under 10⁻⁶ of its largest value | **stage 8's node rule** (spacing 0.1, margin 0.1), put through the same three kernel refinements, must be caught moving its root by more than 10⁻⁴ on each of the three narrowest annuli, where it has three or four nodes |
| N2 the ring limit | N1 passed; the growth rate's relative difference from the exact-history ring of stage 5, at the annulus's own mean radius and total writing rate, falls at every step and **the last is under 0.5%** | the ring at R = 1 instead of the annulus's radius must miss the last annulus by more than 2% |
| N3 every population root resolved | for `A_cold`, `B_cold`, B13, B16 and family B at dE = 0.017 and 0.018, in the family rectangle: the same requirements as N1 | a deliberately coarse rule (spacing 0.3, margin 0.1), with its spacing halved, must be caught moving the root by more than 10⁻⁴ on every population |
| N4 the root-free conclusions survive refinement | for family B at dE = 0.019, `B_mid`, B24, `A_warm` and `B_warm`: stage 8's declared strips, m = 1 to 4, re-swept with the base representation in both partitions and, in the first partition, under three refinements — basis (spacing 0.025, margin 0.5, shifted), quadrature (288 × 96 orbits, 256 angle nodes) and harmonics (±12): the winding number is zero and nothing is located in every rectangle | the same sweeps of the m = 2 strip of family B at dE = 0.018, whose root is the slowest stage 8 resolved, must find exactly one root under every configuration |

In N2 the frequency difference is reported beside the growth difference and carries no bound, and the
observed order between successive annuli is a diagnostic with no bound: I have no analytical expectation for
it, and a bound without one is how L3 failed. **N4 is a statement about the declared domain only**: m = 1 to
4, growth rates of at least 0.006, the declared frequencies. It is not a stability result, and a damped mode
is invisible to it.

## Part F. Stage 8's forecasts, recalculated beside the originals

Family B at dE = 0.010 to 0.024 under the corrected representation, with and without the boundary rows; the
difference of every recalculated root from the archived forecast of aaf0bfe; and the zero-growth crossing
extrapolated from the resolved roots with Re s < 0.03 by polynomials of degree 1, 2 and 3, for the archived
forecasts and for the recalculated ones. **The forecasts of record remain those pushed in aaf0bfe before
anything was simulated.** These are readings. The crossing is an extrapolation of the resolved unstable
branch in every case; no stability boundary is measured by it.

## Part S. What the root is sensitive to

For `A_cold`, `B_cold`, B13 and B16: δs = −(uᴴ δT v)/(uᴴ T′ v) for +1% of the whole response, of each radial
harmonic, of each family of rows, and of each memory time in H(s).

| gate | requirement | negative control the gate must reject |
|---|---|---|
| S1 the formula | for +1% of the whole response it reproduces an actual re-solve to 5% of the movement | the same formula with H(s) held fixed in T′ must miss by more than 50% |

The uniform change of the response that would account for B13's shortfall is a reading, and it sets the
negative control of part E.

## Part E. B13's own eigenfunction, in the time domain

Bodies of B13 drawn with the owner's sampler move in the frozen potential plus ε e^{γt} Re[δC_m(r) e^{i(mθ+βt)}],
where δC_m is **the predicted eigenfunction**, scaled to a largest value of 0.1, and s = γ + iβ **the predicted
root**, with ε = 0.0005, switched on eight e-foldings before t = 0, with +ε and −ε on identical bodies. At t = 0 the
kernel-weighted second harmonic of the difference gives −M(s)c at every node. Projected on the left
eigenvector it is the root the sampled, sharply truncated population would have, to first order:
δs = −(uᴴ δT v)/(uᴴ T′ v) with δT v = α H Qᵀ[M_measured c − M c]. Eight independent batches of 200,000 bodies. A
separate group of three runs on ONE draw of 50,000 bodies — as declared, with half the step, with twice the
amplitude — measures what the step and the amplitude do, free of sampling noise.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| E1 | the mean implied shift of the growth rate is under 2% of the rate, with a standard error under 1%; in the refinement group, halving the step and doubling ε each move the implied shift by under 0.2% of the rate | the prediction scaled by 0.77 — the uniform change that would account for the shortfall (part S) — must imply a shift of more than 10% |

E1 says whether the *linear response of the simulated population to this mode* is what the calculation says it
is. It says nothing about the live model's nonlinear or finite-number behaviour, which is stage 10's subject.

## What this stage cannot establish

Whether the live model's B13 grows at the predicted rate: that needs the seeded experiment. Anything below a
growth rate of 0.006, or outside m = 1 to 4 and the declared frequencies. Any damped mode. Stability of any
population. A cause for B13's shortfall: this stage can exclude candidates — the boundary terms, the basis, the
quadrature, the population's grid, the response to the mode itself — and cannot by itself name one. The source
law is still a label with a declared mass, no energy or angular-momentum budget is carried, and no
observational comparison is made.

## Declared thresholds (read by the driver; the code holds none of its own)

```json
{
  "R": {
        "base": {"n_L": 192, "n_u": 64, "n_eta": 128, "l_max": 8, "node_spacing": 0.05, "node_margin": 0.3, "node_shift": 0.0, "node_clip": [0.35, 2.9], "kernel_rank_cutoff": 1e-10, "edges": true},
        "variants": {"kernel_spacing": {"node_spacing": 0.025}, "kernel_margin": {"node_margin": 0.5}, "kernel_shift": {"node_shift": 0.5}, "rank_more": {"kernel_rank_cutoff": 1e-12}, "rank_fewer": {"kernel_rank_cutoff": 1e-08}, "orbit_quadrature": {"n_L": 288, "n_u": 96}, "angle_nodes": {"n_eta": 256}, "harmonics": {"l_max": 12}, "population_grid": {"refined": true}},
        "stage8_rule": {"node_spacing": 0.1, "node_margin": 0.1},
        "coarse_rule": {"node_spacing": 0.3, "node_margin": 0.1},
        "control_refinement": {"kernel_spacing": {"node_spacing": 0.5}, "kernel_margin": {"node_margin": 0.3}, "kernel_shift": {"node_shift": 0.5}},
        "contour": {"left_panel": 0.015, "first_panel": 0.02, "growth": 1.5, "right_panel": 0.25},
        "singular_gap": 1000.0,
        "relative_residual": 1e-08,
        "polish_movement": 0.0001,
        "winding_phase_step": 0.39269908169872414
  },
  "B": {
        "populations": ["A_cold", "B_cold", "B13", "B16"],
        "reproduction_T": 1e-13,
        "reproduction_root": 1e-09,
        "reproduction_control_min": 1e-09,
        "taper_populations": ["A_cold", "B13"],
        "taper_selections": [["energy", "lower_L", "upper_L"], ["energy"], ["lower_L"], ["upper_L"]],
        "taper_eps": [0.08, 0.04, 0.02],
        "taper_panel": 24,
        "taper_relative_step": 0.001,
        "taper_agreement": 0.01,
        "taper_control_min": 0.5
  },
  "N": {
        "ring_limit": [[0.04, 0.004, 600], [0.02, 0.001, 1200], [0.01, 0.00025, 2400], [0.005, 6.25e-05, 4800]],
        "ring_limit_rect": [0.006, 0.1, -0.2, 0.2],
        "ring_limit_support": 0.0999,
        "refined_annulus": {"n_r_factor": 2, "n_L": 240, "n_vr": 144},
        "each_variant": 0.0001,
        "sum_of_variants": 0.0005,
        "off_node_residual": 1e-06,
        "control_min": 0.0001,
        "control_annuli": [1, 2, 3],
        "ring_limit_last": 0.005,
        "ring_control_min": 0.02,
        "populations": ["A_cold", "B_cold", "B13", "B16", "B_dE0.017", "B_dE0.018"],
        "family_rect": [0.006, 0.2, -0.5, 0.3],
        "root_free": ["B_dE0.019", "B_mid", "B24", "A_warm", "B_warm"],
        "root_free_control": "B_dE0.018",
        "root_free_configs": {"basis": {"node_spacing": 0.025, "node_margin": 0.5, "node_shift": 0.5}, "quadrature": {"n_L": 288, "n_u": 96, "n_eta": 256}, "harmonics": {"l_max": 12}},
        "m_values": [1, 2, 3, 4],
        "strip_re": [0.006, 0.5],
        "strip_im_per_m": 2.0,
        "strip_im_offset": 0.5,
        "strip_im_retrograde": 1.0,
        "rectangle_height": 0.5
  },
  "F": {
        "family_B_dE": [0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024],
        "threshold_fit_below": 0.03,
        "fit_degrees": [1, 2, 3]
  },
  "S": {
        "populations": ["A_cold", "B_cold", "B13", "B16"],
        "relative_change": 0.01,
        "first_order_agreement": 0.05,
        "control_min": 0.5
  },
  "E": {
        "population": "B13",
        "bodies": 200000,
        "batches": 8,
        "refinements": ["half_step", "double_eps"],
        "refinement_bodies": 50000,
        "replay_bodies": 2000,
        "efoldings": 8.0,
        "eps": 0.0005,
        "step": 0.01,
        "peak_field": 0.1,
        "growth_agreement": 0.02,
        "standard_error_max": 0.01,
        "refinement_agreement": 0.002,
        "control_scale": 0.77,
        "control_min": 0.1
  }
}
```
