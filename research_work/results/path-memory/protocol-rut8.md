# RUT-1 stage 8: the linear modes of the constructed populations, predicted before they are measured

Declared before any stage 8 run. Its own protocol file. Stages 4 to 7 stay frozen: their archives keep
reproducing, stage 6's H6 stays failed, and no file of theirs is edited — including `spectrum.py`, in which
this stage found a defect that stage 7 never exercised (below).

## Why this stage exists

Stage 7 found that, with the physics held fixed, an annulus at dE = 0.010 roughly doubles its radial spread
within thirty periods at every body count, while annuli at dE = 0.020 and 0.030 are quiet. It could not say
why, where the threshold lies, or whether the quiet populations hold slower modes: "only the linear mode
calculation around these populations can", and that calculation — declared in stage 6, item 6 of the
owner's handoff — had not been run. This stage runs it.

Stage 7 also recorded two readings of mine that were badly designed: a mode-fit window sized on the wrong
system, and a discreteness exponent that cannot discriminate before saturation. The repair adopted here is
to **predict first**. The linear calculation is archived and pushed before any new simulation is run, the
measurement window of each new run is sized from the prediction, and what counts as agreement is fixed now.

## The calculation

A perturbation of azimuthal number m and time dependence e^{st} in the memory field, delta C_m(r), is a
potential perturbation -delta C_m for the bodies. In action-angle variables of the unperturbed potential
the linearised collisionless response is, for each radial harmonic l,

    delta f_l(J) = F_l(J) psi_l(J) / (nu_l(J) - i s),    nu_l = l Omega_r + m Omega_theta,
    F_l = nu_l df0/dE + m df0/dL|_E,                     df0/dL|_E includes dE_circ/dL = L / r_c^2,

and the field equation is local in time in the Laplace variable: delta C = alpha H(s) K * delta Sigma, with
H(s) = tau_keep / [(1 + s tau_keep)(1 + s tau_form)] and the footprint's azimuthal harmonic
k_m(r, r') = exp[-(r-r')^2/2w^2] I_m^e(r r'/w^2). Expanding delta C_m in the kernel's own harmonics centred on
radial nodes r_b, and collocating at the same nodes,

    T(s) c = 0,     T(s) = P + alpha H(s) M(s),     P_ab = k_m(r_a, r_b),
    M_ab(s) = (2 pi)^2 sum_l int dJ_r dL  G_al(J) G_bl(J) F_l(J) / (nu_l(J) - i s),

with G_al the orbit transform of k_m(r_a, .). T(s) is analytic for Re s > 0, and its roots there are the
growing modes. They are located by the contour method of stage 7 and counted independently by the argument
principle. **The method cannot see a damped mode.** A population with no root above the declared floor has
"no unstable mode faster than the floor in the declared strip". It is not thereby stable.

## What was run before this declaration, and what it showed

A scratch prototype, not archived, run while stage 7's suite was running. It is stated here because the
campaign will reproduce these numbers:

* the orbit library returns each population's mass and four density moments to 1.4x10^-5;
* **`A_cold`, m = 2: one root, s = 0.050876 - 0.052004i** (e-folding 3.13 T0, pattern speed +0.0260), against
  stage 7's exploratory measurement of 0.0480 +/- 0.0071 and +0.0258 +/- 0.0009 at 1,024 bodies. It is
  unchanged to six figures under four-fold refinement of the orbit quadrature, from 2 to 10 radial
  harmonics, and from node spacing 0.1 to 0.075. No root with Re s >= 0.01 at m = 1 or m = 3;
* **`B_cold`, m = 2: s = 0.041315 - 0.077451i** (3.85 T0, +0.0387) against 0.0370 +/- 0.0055 and
  +0.0448 +/- 0.0056;
* `B_mid` and `A_warm`, m = 2: winding number zero everywhere in Re s >= 0.01;
* **through family B the m = 2 rate falls smoothly with temperature**: 0.0413, 0.0386, 0.0355, 0.0322, 0.0285,
  0.0245, 0.0202, 0.0156, 0.0106 at dE = 0.010 to 0.018, and no root above 0.004 at 0.020. The threshold is
  near dE = 0.020, so stage 7's `B_mid` sits at marginal stability;
* **the narrow cold limit approaches the exact-history ring theory quadratically**: against
  `ring_modes.Ring` evaluated at the annulus's own radius and total writing rate, the m = 2 rate differs by
  23.4%, 4.7%, 1.03% and 0.18% at radial widths 0.028, 0.013, 0.0066 and 0.0033 (observed orders 2.2, 2.2,
  2.5). At the narrowest the annulus gives 0.032103 - 0.004171i and the ring 0.032045 - 0.004143i: two
  formulations that share no code agree to 0.2%. The fourth annulus takes an hour and is not in the gate;
* **the response matrix agrees with a time-domain experiment**: test particles under an imposed growing
  perturbation, +eps and -eps on identical bodies, give -M_ab(s) at s = 0.05 - 0.05i, real and imaginary
  parts at every node, to 2.1% of its largest element with 100,000 bodies and to 0.58% with 400,000 — so
  most of the 2.1% was sampling noise, and the gate is sized accordingly;
* a defect: `spectrum.located()` takes the determinant on the contour with the 2x2 formula whatever the
  size of T. Stage 7 only ever gave it 2x2 matrices, so stage 7 is unaffected, and `spectrum.py` is pinned by
  stage 7's archive and is not edited. For an n x n problem the winding number was that of the wrong function:
  it read 0 around `B_cold`'s root. This stage carries its own contour routine, and the defect is its
  negative control (L5).

**The new populations of part M have been predicted in the prototype and have not been simulated.** Nothing
has been measured for them.

## Three statuses, kept separate (the owner's ruling)

| status | question | consequence of failure |
|---|---|---|
| reproduction | does the run reproduce its frozen archive? | the job exits non-zero |
| numerical verification | do the orbit library, the response matrix, the limits and the root search pass their checks, and reject their negative controls? | blocks every dependent conclusion; the job exits non-zero |
| scientific outcome | what the linear theory predicts, and whether the simulations agree | archived as it falls; never affects the exit status |

The thresholds below are machine-readable and the driver reads them from this file.

## Part L. The calculation, verified

Declared discretization: Gauss-Legendre nodes in L over the declared four-width support (192) and in
u = sqrt[2(E - E_circ)] (64); 128 midpoint nodes in the orbit's eccentric-type angle, with the radial angle
and the azimuthal lag from cosine series; radial harmonics l = -8 to 8; kernel nodes every 0.1 from 0.1
inside the smallest pericentre to 0.1 outside the largest apocentre, clipped to [0.35, 2.9]. Overlapping
Gaussians are nearly dependent (the kernel matrix P has a condition number of 2x10^5 at this spacing and
5x10^8 at 0.075), so T(s) is posed in the orthonormal basis of P's subspace above 10^-10 of its largest
eigenvalue: the same roots, a determinant of order one, and a size set by the kernel's bandwidth. In the
prototype the `A_cold` root is then identical to seven figures at node spacings 0.1, 0.075, 0.05 and 0.025.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| L1 the orbit library | through the orbits, each population's mass and its moments R, R^2, 1/R and a Gaussian window agree with the solver's density to 5x10^-5; for the most nearly circular orbit of every L, Omega_r and Omega_theta agree with the epicyclic and circular frequencies of the same potential to 10^-3 | the same comparison with the memory field removed from the potential must fail |
| L2 the response matrix, in the time domain | 200,000 bodies drawn with the owner's sampler move in the frozen potential plus an imposed perturbation eps e^{gamma t} psi(r) cos(m theta + beta t), switched on six e-foldings before t = 0, once with +eps and once with -eps on identical bodies so that shot noise cancels. At t = 0 the kernel-weighted m-th harmonic of the difference must equal -M_ab(s) eps, column by column, to 3% of the largest element — for `A_cold` at s = 0.05 - 0.05i and `A_warm` at s = 0.03 - 0.10i | the prediction with only the l = 0 harmonic kept must fail |
| L3 the ring limit | three annuli at (dL, dE) = (0.04, 0.004), (0.02, 0.001), (0.01, 0.00025), each solved for a 9.99% mean support: the m = 2 root against the exact-history ring of stage 5 at the annulus's own mean radius and total writing rate. The relative difference in growth rate must fall at every step, the last must be under 2%, and the observed order in the radial width must lie between 1.5 and 2.5 | the ring evaluated at R = 1, stage 5's radius, instead of the annulus's must fail the 2% |
| L4 convergence | for `A_cold` and for family B at dE = 0.016: doubling the L and u nodes, doubling the angle nodes, raising the harmonics to 12, refining the kernel nodes to 0.075 and widening their range by 0.2 each move the m = 2 root by under 10^-4 relative | keeping only l = 0 must be caught moving it by more |
| L5 completeness | in every sub-rectangle of every declared strip: the winding number of det T, taken with the full n x n determinant, equals the number of located roots with multiplicity; a located root has a singular gap of at least 10^3, a polished residual relative to the contour's median |det| under 10^-8, and moves under 10^-4 on polishing; and a second partition sharing no cut returns the same roots to 10^-6 | the 2x2 determinant formula applied to the n x n problem must be caught returning the wrong winding number around `B_cold`'s root |

The declared strip is Re s from 0.006 to 0.5 and Im s from -(2 m + 0.5) to +1, for m = 1 to 4, for the five
stage 7 populations and the three new ones: pattern speeds from -1/m to 2 + 1/2m, which covers every prograde
pattern faster than any orbit in these populations (the largest orbital frequency is 1.63) and a margin of
retrograde ones. It is cut into sub-rectangles 0.5 high; the second partition is offset by half a rectangle.
The contour's left edge has Gauss-Legendre panels of 0.015 beside the axis, where the discretised response
has its poles, and the horizontal edges are graded from 0.02 at the left corner.

**The floor of 0.006 (an e-folding of 26.5 T0) is a limit of this discretization, measured in the
prototype:** against a quadrature of 384 x 96 orbits, the declared 192 x 64 reproduces M(s) to 5x10^-4 at
Re s = 0.01, to 3x10^-3 at 0.006, to 1.4x10^-2 at 0.004 and to 7x10^-2 at 0.002. Below the floor the
resonant denominators are not resolved and nothing is claimed.

## Part P. Predictions, pushed before anything new is simulated

1. Every root in the declared strip for the five stage 7 populations. For `A_cold` and `B_cold` this is a
   **retrospective** comparison with stage 7's exploratory growth-phase fit, and is labelled so: I have seen
   both numbers.
2. Family B — mass, alpha, footprint and response times fixed — at dE = 0.010 to 0.024 in steps of 0.001: the
   m = 2 root, and the threshold dE* from a quadratic least-squares fit of Re s against dE over the roots with
   Re s < 0.03, with its spread under dropping each point in turn.
3. **Three new populations, never simulated: family B at dE = 0.013, 0.016 and 0.024.** For each: the m = 2
   growth rate and pattern speed, or the statement that the strip holds no root.

The predictions are archived and committed, and the commit is pushed, before part M is run.

## Part M. Measurement of the three new populations

The three new populations pass through stage 7's construction and sampling gates unchanged (V1, V2, V3, V5,
and V6 with four realizations of the declared draw and of each of its three negative controls), using stage
7's code by import. They are then evolved with stage 7's
simulator and records: live and frozen from byte-identical bodies, 256 bodies with 8 realizations and 1,024
with 4. The horizon is sized on the prediction — six periods for the shot-noise imprint plus five predicted
e-foldings, rounded up to a multiple of ten: **40 T0 at dE = 0.013, 50 T0 at 0.016**, and 40 T0 at 0.024
where no growth is predicted.

**The growth-phase rule, declared here.** It is the rule stage 7 used as exploratory; its three numbers were
chosen there after seeing the cold annuli, and they are fixed now before these populations exist. For each
population and body count: the geometric mean over realizations of |C_2|/C_0 on the ring at the initial mean
radius; a fit of log amplitude and of phase, per realization, from 6 T0 to the first time that mean reaches
half its maximum over the run. If that window is shorter than 2 T0 there is no growth phase.

**Declared readings.**

* *Confirmed*, for a population predicted unstable: a growth phase exists at both body counts, and at 1,024
  bodies the ensemble-mean rate lies within three standard errors of the prediction or within 15% of it, and
  the pattern speed within three standard errors or within 0.01.
* *Contradicted*: no growth phase at 1,024 bodies; or a rate differing from the prediction by more than three
  standard errors **and** more than 30%.
* *Unresolved*: anything else.
* For the population predicted to hold no root: *confirmed* if there is no growth phase at either body count;
  *contradicted* if a growth phase exists at 1,024 bodies with a rate above 0.006 by three standard errors;
  otherwise *unresolved*.

D_Q is reported for the new populations exactly as in stage 7, with no reading attached to it here.

## What this stage cannot establish

The floor. Below a growth rate of 0.006, an e-folding of 26.5 T0, this discretization does not resolve the
resonances, and the cold
ring's spectrum shows that growth times of hundreds of periods exist in this model. Damped modes are
invisible to a search of the right half-plane, so marginal populations are described only from the unstable
side. The calculation is for one footprint, one pair of response times, one angular-momentum width and one
radius; the threshold it locates is a threshold in dE for that family. The source law is still a label with
a declared mass, no energy or angular-momentum budget is carried, and no observational comparison is made.

## Declared thresholds (read by the driver; the code holds none of its own)

```json
{
  "L": {"n_L": 192,
        "n_u": 64,
        "n_eta": 128,
        "l_max": 8,
        "node_spacing": 0.1,
        "node_margin": 0.1,
        "node_clip": [0.35, 2.9],
        "kernel_rank_cutoff": 1e-10,
        "m_values": [1, 2, 3, 4],
        "strip_re": [0.006, 0.5],
        "strip_im_per_m": 2.0,
        "strip_im_offset": 0.5,
        "strip_im_retrograde": 1.0,
        "rectangle_height": 0.5,
        "contour": {"left_panel": 0.015, "first_panel": 0.02, "growth": 1.5, "right_panel": 0.25},
        "moments": 5e-05,
        "epicyclic": 0.001,
        "time_domain_bodies": 200000,
        "time_domain_eps": 0.002,
        "time_domain_efoldings": 6.0,
        "time_domain_step": 0.01,
        "time_domain_relative": 0.03,
        "time_domain_cases": {"A_cold": [0.05, -0.05], "A_warm": [0.03, -0.1]},
        "ring_limit": [[0.04, 0.004, 600], [0.02, 0.001, 1200], [0.01, 0.00025, 2400]],
        "ring_limit_rect": [0.006, 0.1, -0.2, 0.2],
        "ring_limit_support": 0.0999,
        "ring_limit_last": 0.02,
        "ring_limit_order": [1.5, 2.5],
        "convergence": 0.0001,
        "convergence_populations": ["A_cold", "B16"],
        "singular_gap": 1000.0,
        "relative_residual": 1e-08,
        "polish_movement": 0.0001,
        "partition_agreement": 1e-06,
        "winding_phase_step": 0.39269908169872414},
  "P": {"family_B_dE": [0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024],
        "family_rect": [0.006, 0.2, -0.5, 0.3],
        "threshold_fit_below": 0.03,
        "new_populations": {"B13": 0.013, "B16": 0.016, "B24": 0.024}},
  "M": {"bodies": [256, 1024],
        "realizations": [8, 4],
        "horizon": {"B13": 40.0, "B16": 50.0, "B24": 40.0},
        "growth_start": 6.0,
        "growth_end_fraction": 0.5,
        "growth_shortest": 2.0,
        "confirm_sigma": 3.0,
        "confirm_relative": 0.15,
        "confirm_pattern": 0.01,
        "contradict_relative": 0.3,
        "stable_rate_floor": 0.006}
}
```
