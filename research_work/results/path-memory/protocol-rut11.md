# RUT-1 stage 11: the quiet region — the response below stage 8's floor, wider mode coverage, and the answer to a controlled disturbance

Declared before any stage 11 run, and pushed with the code it names: `resonant_response.py`, `rut11_tasks.py`,
`rut11.py`, `rut11_checks.py`. Stage 9's response, stage 10's quiet sample, stage 7's simulator and population
builder and the owner's sampler are used **by import, unchanged**.

## Why this stage exists

The owner's review of 3a80fec set four pieces of work in order. Stage 9 did the first — the repaired response
calculation, verified on every population and on the ring limit. Stage 10 did the second — a deliberately
seeded warm mode grows at the predicted rate. This is the third:

> **Establish the quiet region's scope.** Improve resonant integration below the current growth cutoff; expand
> relevant mode coverage; measure response to controlled disturbances. *Whether the warmer population is
> genuinely stable, weakly unstable or simply untested on longer timescales.*

and it is the last clause of the milestone: *identifies a supported quiet region with a stated timescale and
mode-domain bound*.

Three things stand in the way of such a statement, all of them named in the review.

**The floor.** Stage 8's solver stops resolving growth below Re s = 0.006, an e-folding time of 26.5 reference
periods, so "no root" has never meant more than "no fast root". That floor is not a tolerance anyone chose. The
Gauss-Legendre rule of stages 8 and 9 replaces the continuum of orbital frequencies by 12,288 of them, so the
discretised response has a **pole at s = −iν for each one**, all sitting on the imaginary axis a mean spacing of
about 1.2×10⁻⁴ apart. Where Re s is not comfortably above that spacing the sum is not a quadrature of anything.

**The mode domain.** The search covered m = 1 to 4 and a finite strip of complex frequencies. Nothing was said
about m = 5 and above.

**The horizon.** B24's forty periods exclude nothing slower than the run itself, and nothing was ever done *to*
a quiet population to see what it does when pushed.

## The construction

**The resonant integral, done in closed form.** On each cell of a refined (L, u) grid the frequency ν and the
numerator are interpolated **linearly from the vertices**, and the integral of the result against 1/(ν − i s) is
written down exactly. For a triangle with vertex frequencies n₁ ≤ n₂ ≤ n₃ each barycentric coordinate needs only

    L0(x) = log(1+x)/x,   L1(x) = [x − log(1+x)]/x²,   L2(x) = [log(1+x) − x + x²/2]/x³,

taken by their series where the closed forms cancel. What this evaluates is the exact integral of a
piecewise-linear model of the true integrand, so **the pole is integrated rather than sampled** and the error is
second order in the cell size *uniformly in Re s*. Two refinements are carried and Richardson-extrapolated; the
difference between the pair and its coarser member is computed on every contour and is the error the rule
declares about itself.

Only the harmonics that can resonate in the rectangle being searched are treated this way — the ones whose
frequency range comes within a declared margin of −Im s. Every other harmonic keeps stage 9's rows, whose poles
are then at least that margin from the contour, and the margin actually achieved is recorded for each
rectangle. The orbits, the kernel basis and the orbit transform are stage 9's: the refined grid is reached by
barycentric interpolation from the Gauss-Legendre nodes, where those quantities are known to spectral accuracy,
so nothing is recomputed and nothing new is approximated.

This buys growth rates, not damped modes. For Re s > 0 the integral above **is** the Laplace transform, with no
continuation. A damped mode lives across the cut on Re s = 0 and is not attempted here; that is stated again
below among the things this stage cannot establish.

**The operator bound.** T(s) = 1 + α H(s) M(s) in the orthonormal basis. Wherever the spectral norm of
α H(s) M(s) stays below one, T is invertible and **no mode exists there at any resolution** — a statement about
the operator, not about a search. It is evaluated on a declared grid over the whole strip, from the new floor to
Re s = 0.5, and the grid is refined as its own control.

**The controlled disturbance.** A quiet sample as in stage 10 — `quarter` bodies from the owner's sampler and
their images under rotation, four-fold for an even m and two-fold for an odd one, exact in floating point, so
the sample carries no harmonic of the one being disturbed. Into its field goes

    δC = ε Re[ k_m(R, r) e^{imθ} ],   δE = δC / τ_keep,

the memory kernel centred on the population's own mean radius: a declared shape with no free parameter but its
amplitude, and the excitation the settled relation C = τ_keep E requires. The bodies are **not** moved — there
is no mode to prepare them in. What is measured is what the live system does with the kick.

## What was run before this declaration, and what it showed

Scratch prototypes, on the exact configurations declared below. Sample sizes are given; nothing here is a
result of this stage.

1. **The rule against a known integral** (`reference_integral`, one integral, five growth rates). Richardson
   pair (192×64, 384×128): errors 9×10⁻¹⁰ to 5×10⁻⁸ at growth rates from 3×10⁻² down to 10⁻⁶. Stage 9's
   Gauss-Legendre rule on the same integral: 1.2×10⁻¹⁰ at 3×10⁻², then 7×10⁻⁴ at 6×10⁻³, **6.1×10⁻² at 10⁻³,
   192% at 10⁻⁴ and 596% at 10⁻⁶**. Each doubling of the refined grid divides the rule's error by 4.00.
2. **A second toy with a band edge** (one integral, five positions, three growth rates) — where ν has an
   interior extremum the rule's convergence degrades: at 10⁻⁵ the Richardson pair reaches only 1×10⁻³ against
   4×10⁻⁴ away from the edge. At 10⁻³ and above it is 10⁻⁶ or better. This is why the floor below is 5×10⁻⁴ and
   not lower, and why the rule's declared error is monitored on every contour.
3. **The rule on B24, m = 2** (one population, eleven values of s). Against the Richardson pair of (384×128,
   768×256): the pair (192×64, 384×128) differs by 2.9×10⁻⁸ at Re s = 0.05, 4.4×10⁻⁷ at 10⁻³, 3.0×10⁻⁵ at
   10⁻⁴; stage 9's rule differs by 4.9×10⁻⁷ at 0.05, 1.1×10⁻⁴ at 10⁻³, **6.6×10⁻³ at 10⁻⁴, and 30% at
   10⁻⁴ − 0.3i**, where the resonance sits mid-band. One evaluation of the pair costs 0.36 s against 0.06 s for
   stage 9's rule. Dropping stage 9's row pruning changes the response by 2×10⁻¹⁵.
4. **The same root by both rules** (B13, m = 2, stage 9's family rectangle, one rectangle). Stage 9:
   0.03216180 − 0.08786641i. The analytic rule: 0.03216180 − 0.08786641i, a relative difference of
   **2.7×10⁻⁸**.
5. **A root below stage 8's floor** (B_dE0.019, m = 2, box Re s ∈ [5×10⁻⁴, 6×10⁻³], Im s ∈ [−0.5, 0], one
   rectangle). The analytic rule locates **s = 0.0052280 − 0.1083226i**, winding number 1 in agreement, an
   e-folding time of 30.4 periods. The finer pair moves it by 6.1×10⁻⁸. The rule's declared error on the floor
   edge is 2.8×10⁻⁴ for the pair and 1.1×10⁻³ for its coarse member — a ratio of 4.0, the second-order
   behaviour the extrapolation assumes. **Stage 9's rule on the same box returns a winding number of 3 while
   locating one root, and puts it at 0.0057893, 10.7% away.** This population is one of the five stage 9
   reported no root in, having searched only above 0.006.
6. **The contour may be coarser** (same box, four panel settings). With the analytic rule the response has no
   pole beside the contour, and the located root is identical to six digits at left panels of 0.015, 0.05, 0.1
   and 0.25, with the winding number 1 throughout and the singular gap falling from 1.3×10⁵ to 2.0×10⁴ — far
   above the 10³ stage 8 requires. 0.1 is declared, and stage 9's 0.015 is a declared control.
7. **The operator norm** (B24, twelve azimuthal numbers, a grid over the whole strip above 0.006). The maximum
   of ‖αH(s)M(s)‖ is 0.248 at m = 1, 0.546 at m = 2, then **0.0095 at m = 3 and falling monotonically to
   8×10⁻⁴ at m = 12**. Carried down to the new floor with the analytic rule the maximum rises: **3.07 at m = 2
   and 0.308 at m = 3**, both at Re s = 5×10⁻⁴, Im s = 0 — a nearly stationary pattern. Halving the grid
   spacing in Im s did not move either. So m = 2 must be searched and cannot be bounded; m ≥ 3 can be.
8. **A controlled disturbance** (B24, m = 2, 1,024 bodies, one realization, forty periods) — the run this
   declaration's horizon and window are taken from. The disturbance decays. |C₂|/C₀ falls from 1.002×10⁻⁵ to 1.744×10⁻⁶ over sixteen periods — a factor of 5.7 — with a fitted rate of −0.018 over periods 5 to 16 and −0.027 over 10 to 16, a root-mean-square residual in the log amplitude of about 0.1, and a bump near period 10: this is phase mixing and a decaying mode together, not one clean exponential, which is why the window below starts at 10 and why the plain ratio of end to start is reported beside the fitted rate. The pattern barely turns (a speed of +0.003). The four-fold sample holds: |C₁|/C₀ and |C₃|/C₀ stay at 2×10⁻¹⁶ for the whole run. **m = 4 is not folded away** and rises from 3.2×10⁻⁷ to 8.9×10⁻⁵ as the 64 independent bodies write their own shot noise into the field — expected, reported as a diagnostic, and the reason the disturbed harmonic is the only one read.

## Three statuses, kept separate (the owner's ruling)

| status | question | consequence of failure |
|---|---|---|
| reproduction | does the run reproduce its frozen archive? | the job exits non-zero |
| numerical verification | does the analytic rule do what it claims — reproduce stage 9 where stage 9 is valid, hold its declared error on every contour it is used on, count what it locates, and would the apparatus see growth if there were any? | blocks every dependent conclusion; the job exits non-zero |
| scientific outcome | which populations hold no mode above the new floor, where the family's boundary is bracketed, and what a disturbance does to a quiet population | archived as it falls; **never** affects the exit status |

No population may be called stable. Which populations turn out to be quiet is an outcome, not a gate: the gates
below ask only whether the search was sound.

## Part A. The rule

| gate | requirement | negative control the gate must reject |
|---|---|---|
| A1 on a known integral | at every declared growth rate the Richardson pair is within 10⁻⁶ of a reference good to fourteen digits | stage 9's quadrature of the same integral must be seen to miss it by more than 5% at the rates at or below 10⁻³ |
| A2 against stage 9 | on every population with a root stage 9 resolved, the analytic rule's root differs from stage 9's by less than 10⁻⁵ relative, and moves less than 10⁻⁵ under the finer pair | stage 9's own quadrature, applied in the low box where its poles lie, must be seen to disagree with the analytic rule by more than 1% |
| A3 the declared error | on every contour used in parts L and T, the rule's estimate of its own error stays under 2×10⁻³ | — (a budget, checked on every contour, not a comparison) |

## Part L. The low box

Re s from 5×10⁻⁴ — an e-folding time of **318.3 reference periods** — to stage 8's floor of 6×10⁻³, over the
same strip in Im s stage 8 declared, m = 1 to 4, in two partitions offset by half a rectangle.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| L1 the search is sound | in every rectangle, of every population searched, in both partitions: the winding number equals the number of roots located, the rule's declared error is under 2×10⁻³, and the nearest remaining Gauss-Legendre pole is at least 0.2 from the rectangle | stage 9's quadrature on the declared control rectangles must be seen to break that agreement or misplace the root by more than 1% |
| L2 the search can see | the root the prototype found in B_dE0.019 is located in both partitions, with agreeing counts, and the two agree to 10⁻³ relative | — (this gate *is* the detectability control the readings rest on) |

## Part M. The mode domain

| gate | requirement | negative control the gate must reject |
|---|---|---|
| M1 the grid is adequate | for every population and azimuthal number in the declared set, refining the grid changes the maximum of ‖αH(s)M(s)‖ by less than 0.05 | B13 at m = 2, which holds a root, must be seen to reach 1 |

Where that maximum stays below 0.9 the reading may say no mode exists at that m at any resolution. Where it
does not, the reading says only what was searched.

## Part T. The family's boundary, bracketed

Family B members at the declared widths, each searched for its fastest m = 2 root from the new floor upward.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| T1 each rung | counts agree, the declared error is under 2×10⁻³, and every located root moves less than 10⁻³ relative under the finer pair | — (L1's control covers the rule; this is a refinement condition) |

The reading is a **bracket**: the warmest member with a resolved root, and the coolest with none above the
floor. It is not an extrapolation, and it says nothing about growth below the floor. Stage 8's extrapolated
crossings (linear 0.02044, quadratic 0.019927, cubic 0.01991) are reported beside it, never replaced by it.

## Part D. The controlled disturbance

| gate | requirement | negative control the gate must reject |
|---|---|---|
| D1 the disturbance is as declared | at t = 0 the ring's \|C_m\|/C₀ is the declared target to 2%, and the source the folded sample would write carries \|S_m\|/S₀ under 10⁻¹⁰ **at the harmonic being disturbed**, so what is watched can only be what was put there | an ordinary draw of the same size must be seen to carry that source above 10⁻³ |
| D2 the apparatus would see growth | B13, given the same disturbance at m = 2, grows, and its fitted rate over the declared window agrees with stage 9's recalculated root to 15% or three standard errors | — (this gate *is* the control every quiet reading rests on) |
| D3 completeness | every declared run reaches its horizon | — |

## The declared reading

For each population: the azimuthal numbers searched and the floor searched to, expressed as an e-folding time;
the azimuthal numbers bounded by the operator norm; the roots found; and what the disturbance did — the fitted
rate over the declared window with its spread over realizations, and the plain ratio of the amplitude at the
end of the horizon to the amplitude at the start.

A population with no root found and a disturbance that does not grow is described as **"no mode above an
e-folding rate of 5×10⁻⁴ in the searched domain, and a declared disturbance that did not grow over forty
periods"**. It is not described as stable.

## What this stage cannot establish

Damped modes: the rule is used only for Re s > 0, and a decaying mode needs a continuation across the cut that
is not attempted. Anything outside the declared strip in Im s. Azimuthal numbers above the declared set, except
through the operator norm, which is itself only as good as its grid. Growth slower than 5×10⁻⁴, which the
band-edge prototype shows is where this rule's convergence begins to degrade. The sharply truncated
distribution function is still the model; nothing here changes it. And the disturbance is one declared shape at
one declared amplitude: another shape may excite something this one does not.

## Declared thresholds (read by the driver; the code holds none of its own)

```json
{
    "R": {
        "base": {
            "n_L": 192,
            "n_u": 64,
            "n_eta": 128,
            "l_max": 8,
            "node_spacing": 0.05,
            "node_margin": 0.3,
            "node_shift": 0.0,
            "node_clip": [0.35, 2.9],
            "kernel_rank_cutoff": 1e-10,
            "edges": true
        },
        "refines": [1, 2],
        "control_refines": [2, 4],
        "pattern": "alternate",
        "resonant_margin": 0.3,
        "pole_margin_min": 0.2,
        "gamma_min": 0.0005,
        "gamma_stage8": 0.006,
        "contour": {
            "left_panel": 0.1
        },
        "spread_max": 0.002
    },
    "A": {
        "toy": {
            "a": 0.3,
            "b": 2.5,
            "c": 1.2,
            "omega": 0.42,
            "dL": 0.06,
            "dE": 0.024,
            "reach": 4.0,
            "n_L": 192,
            "n_u": 64,
            "growth_rates": [
                0.03,
                0.006,
                0.001,
                0.0001,
                1e-06
            ],
            "accuracy": 1e-06,
            "control_below": 0.001,
            "control_min": 0.05
        },
        "roots": ["A_cold", "B_cold", "B13", "B16", "B_dE0.017", "B_dE0.018"],
        "rect": [0.006, 0.2, -0.5, 0.3],
        "stage9_agreement": 1e-05,
        "refinement_agreement": 1e-05,
        "control_min": 0.01
    },
    "L": {
        "quiet": ["B_mid", "B24", "A_warm", "B_warm"],
        "control_population": "B_dE0.019",
        "control_rectangles": [8],
        "control_agreement": 0.001,
        "m_values": [1, 2, 3, 4],
        "strip_im_per_m": 2.0,
        "strip_im_offset": 0.5,
        "strip_im_retrograde": 1.0,
        "rectangle_height": 0.5
    },
    "M": {
        "m_values": [1, 2, 3, 4, 5, 6, 8, 10, 12, 16],
        "refined_m": [3, 4, 5, 6],
        "bounded_from": 3,
        "norm_bound": 0.9,
        "refine_rise": 0.05,
        "re_max": 0.5,
        "n_re": 6,
        "im_spacing": 0.25,
        "refine_factor": 2,
        "analytic_at_floor": true,
        "control_population": "B13",
        "control_m": 2
    },
    "T": {
        "dE_values": [0.019, 0.0193, 0.0195, 0.0197, 0.0198, 0.0199, 0.0201, 0.0203],
        "re_max": 0.05,
        "im_lo": -0.5,
        "im_hi": 0.0,
        "refinement_agreement": 0.001,
        "stage8_extrapolations": {
            "linear": 0.02044,
            "quadratic": 0.019927,
            "cubic": 0.01991
        }
    },
    "D": {
        "quiet": ["B_mid", "B24", "A_warm", "B_warm"],
        "control_population": "B13",
        "control_predicted_rate": 0.0321618,
        "control_agreement": 0.15,
        "m_primary": 2,
        "odd_population": "B24",
        "odd_m": [1, 3],
        "odd_realizations": 2,
        "quarter": 256,
        "realizations": 4,
        "target": 1e-05,
        "horizon": 40,
        "window": [10, 40],
        "windows_reported": [
            [5, 40],
            [20, 40],
            [30, 40]
        ],
        "record_every": 0.25,
        "h_max": 0.01,
        "eta": 0.01,
        "grids_per_w": {
            "fine_grid": 8,
            "coarse_grid": 1.25
        },
        "variant_population": "B24",
        "variants": ["half_step", "fine_grid", "coarse_grid", "no_excitation", "ordinary_draw"],
        "variant_target_small": 1e-06,
        "amplitude_agreement": 0.02,
        "other_harmonic_max": 1e-10,
        "control_harmonic_min": 0.001
    },
    "C": {
        "replay_prefix": 2
    }
}
```
