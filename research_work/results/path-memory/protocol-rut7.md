# RUT-1 stage 7: the spectrum certified to the axis, the population drawn faithfully, and live compared against frozen with a measured uncertainty

Declared before any stage 7 run. Its own protocol file. Stages 4, 5 and 6 stay frozen: their archives keep
reproducing, stage 6's failed H6 stays failed, and what they got wrong is corrected in their reports, not
in their numbers.

## Why this stage exists

The owner's review of 72aef44 found that **the stage 6 sampler does not draw the population the solver
constructs**: in the variables (R, L, v_r) the phase-space measure is

    dM = f(E, L) R dR dtheta dv_r dv_theta = f(E, L) dR dtheta dv_r dL,        since dv_theta = dL / R,

with no factor of R left over, and `sample()` multiplies its node weights by `2 pi R`. Confirmed before
this declaration to five digits. The same review ruled that a run's status has three separate meanings,
asked for the distribution's support and discretization accuracy to be established, for the temperature
comparison to hold the physics fixed, for dispersion to mean what is left after organised motion is
removed, and for H6 to be replaced by a construction check and a separate physical experiment.

While this stage was being prepared the owner published a corrected sampler of their own, with its checks
and a six-item handoff (`research_work/annulus_sampling`, protocol 2f1d5ed, code bde7a7c). **This stage is
the execution of that handoff**: their sampler is the sampler of record at every source-sampling call
site, the historical `equilibrium.py` and `rut6.py` stay byte-identical, and H6 stays failed.

I also had stage 6 audited and re-verified the consequential findings myself. They are recorded in
`report-rut6.md`; the ones this stage has to repair are:

* **The certified spectrum stopped at Re s = 0.001.** Below that edge the two-stage ring has unstable roots
  at every m from 1 to 8 — the fastest a static lopsided m = 1 mode at e-folding 285 T0, on a branch
  continuation cannot reach — and Beyn's own moments had located them before the rectangle filter threw
  them away. Stage 5's continuation returns a duplicated root at every m >= 1.
* **H2 compared Newton's fixed points, not what the contour solver located.** Unpolished, the m = 2 root
  moves 5.8e-6 against a declared 1e-6. The rank rule silently dropped one of six eigenvalues in the
  solver's own transcendental test. The region test and the deduplication ran before the polish.
* **Both archived annulus runs carry a growing, inertially stationary m = 2 mode** (e-folding 8.2-9.1 T0
  warm, 5.3-5.7 T0 cold) that `_drift`'s azimuthal averages could not see; the cold run's "+14%
  dispersion" is that mode, its residual dispersion changing by -0.6%.
* **"10% support" was read at the argmax grid node beside the field's gradient zero**: a sawtooth in alpha
  with several roots, moving 10-35% with the grid. The mass-weighted support was 13.4% and 20.3%.
* **`circular_energy()` inverts a map that is not monotonic** in these potentials, so the wide annuli are
  not the declared distribution; and with five-width windows the most extended populated orbit reaches
  r = 2.9 to 3.0, past the simulator's limit.
* **The frozen-field stationarity statistic cannot tell an equilibrium from a kicked sample**, because
  window averages in a static axisymmetric field are near-invariants of any start.
* **Two declared items were changed without being named**: H3's tolerance (coded 1e-5, declared 1e-6) and
  H6's observable (support replaced by radial spread). H4's amplitude-scaling check could not fail.

## Two safeguards, adopted from here on

**The thresholds below are machine-readable, and the driver reads them from this file.** A tolerance
cannot be coded differently from how it was declared, because the code holds no tolerance of its own.

**Every numerical verification gate carries a negative control** — a deliberately broken input that the
gate must reject. A gate whose negative control passes is itself recorded as failed. Several stage 6 checks
could not have failed; a check that cannot fail verifies nothing.

Any departure from this protocol is written into the archive's `deviations` list and into the report,
named as a departure, whether or not it changes an outcome.

## Three statuses, kept separate (the owner's ruling)

| status | question | consequence of failure |
|---|---|---|
| reproduction | does the run reproduce its frozen archive? | blocks reproducibility claims; the job exits non-zero |
| numerical verification | do identities, samplers, integrators and solvers pass their correctness checks? | blocks every dependent physical conclusion; the job exits non-zero |
| scientific outcome | does the model achieve the physical objective? | a legitimate result, archived as it falls; never affects the exit status |

## Part E. The cold-ring spectrum, certified down to the axis

The search strip is **Re s in [-0.008, 1.5], |Im s| <= 6**, for every m from 0 to 16, for the instantaneous,
one-stage and two-stage rings and the no-attraction control at the stage 5 configuration (32 writers,
w/R = 0.2, 10% label). Its left edge lies to the LEFT of the imaginary axis, halfway to the nearest poles
of B_m(s) at Re s = -1/tau_keep = -0.0159, so marginal and slowly growing roots are *inside* the contour
and are located and classified, not excluded. There is then no lower cutoff on the growth rate other than
the resolution floor: a located root is called unstable if Re s > 1e-7 and marginal if |Re s| <= 1e-7.

The strip is cut into sub-rectangles in Im s so that each holds few roots, with panel lengths of 0.02 along
the left edge (node spacing about 0.001, against a distance of 0.008 to the nearest root or pole). In each
sub-rectangle, from one set of contour evaluations:

* the eigenvalues are located by Beyn's method with **shifted and scaled moments** and the rank chosen at
  the **largest gap** in the singular spectrum, which is archived;
* the number of zeros of det T is counted independently by the **argument principle**, with the contour
  refined until no phase increment exceeds pi/8 and the smallest |det| on the contour reported. The
  winding number counts every zero whatever its residue, which is what the moment method cannot promise;
* located roots are **polished first and tested afterwards**: region, deduplication and residual are all
  evaluated on the polished value, the polish step is bounded by the distance to the contour and to the
  nearest other root, and the distance each root moved is archived.

**A declared rule for roots near an edge, so that applying it is not a deviation:** if a located root, or a
minimum of |det| along the contour, puts a root within 0.001 of any edge of a sub-rectangle, that
sub-rectangle is re-run with the edge moved 0.002 further from the root, and both runs are archived.

Continuation is kept as a cross-check only, with a pairwise-distinctness assertion and the lost partner
recovered by deflation; its disagreement with the certified list is reported and is not a gate.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| E1 known spectra | exp(-s) = 1/2: three roots, **six with multiplicity**, raw error < 1e-8. Zero-lag quartic: raw error < 1e-8 against its companion roots. Kepler control: its four roots 0, 0, +/- i Omega are **located** and none is classified unstable | a solver run with the stage 6 rank rule (1e-8 of the largest singular value, unscaled moments) on the transcendental case must be caught returning five |
| E2 raw robustness | the UNPOLISHED located roots move < 1e-6 under doubled quadrature and under a shifted contour, for one-stage m = 2, 5 and two-stage m = 1, 2, 3 | the stage 6 quadrature (12 panels of 16 nodes on the stage 6 rectangle) must be caught moving the two-stage m = 2 root by more than 1e-6 |
| E3 completeness | in every sub-rectangle, for every rung and m: Beyn's count with multiplicity equals the winding number; no moment saturation; singular gap >= 1e3; every root's polished residual < 1e-10 with polish movement < 1e-4; and a second, shifted partition of the strip returns the identical root set to 1e-8 | a synthetic T(s) holding a root whose residue is 1e-9 of its neighbours' must be reported as a count mismatch, not silently dropped |
| E4 edge handling | a root placed 1e-6 inside a sub-rectangle's edge is found | the stage 6 order of operations (filter, then polish) must be caught losing it |

## Part S. The population, defined, and drawn faithfully

    f(E, L) = A exp[-(L-L0)^2/2dL^2] exp[-(E-E_circ(L))^2/2dE^2]    for |L-L0| <= 4 dL and
                                                                       0 <= E - E_circ(L) <= 4 dE,
            = 0                                                       otherwise,

with `E_circ(L)` the **global minimum over radius of Phi(r) + L^2/2r^2** — single-valued for every L,
a lower bound on the energy of every orbit of that L, and a function of L alone, which is all Jeans'
theorem needs. The truncation is part of the definition, not an integration window: a Gaussian factor in
energy is not by itself a finite-mass guarantee. `A` is fixed by a declared total mass of 1, so the writing
coupling `alpha` is a physical constant of a family and not something the distribution's width rescales.

The descriptor of a population is its **mass-weighted support** — the mean and the spread, over the
population, of the inward memory acceleration over the Newtonian — solved for with a bracketed root
finder. The value at one radius is reported and is not used.

Two families, L0 = 1, dL = 0.06, w = 0.2, tau_keep = 10 T0, tau_form = 3 T0:

* **A, matched support:** dE = 0.010 and 0.030, each with alpha solved so the mass-weighted mean support
  is 10%. Asks whether supported populations exist across temperatures.
* **B, fixed physics:** alpha calibrated once, at dE = 0.020, to a 10% mean support, then held fixed with
  the mass, the footprint and the response times while dE = 0.010, 0.020, 0.030. The support that results
  is reported, not forced. Temperature is the only thing that changes in this family.

**Bodies are drawn with the owner's corrected node sampler** (`annulus-phase-space-v2`: nodal masses
4 pi dR dL w_vr f, radial trapezoid weights, the density routine's own equal L weights, both signs of v_r,
uniform angle), built on the redefined distribution with its declared four-width support. It is the
sampler of record for every run in this stage. An **independent continuous sampler** is its cross-check:
rejection in (R, L, u) with v_r = a + u(b - a) and target density f(E, L)(b - a) — no factor of R, (b - a)
the Jacobian of u — which places no body on a quadrature node. Two samplers built differently, agreeing
with the same quadrature and with each other, is a stronger statement than either alone.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| V1 the measure | the owner's `verify_distribution` on every population: radial marginal, total mass and the first two radius moments agree with the density routine to 1e-12, and the source is left unmutated | the stage 6 weights, with their extra factor of R, must be rejected, with a mean-radius shift equal to Var(R)/Mean(R) to 1e-12 |
| V2 drawn moments | 200,000 bodies per population **from each of the two samplers**: means of R, L and E within 4 standard errors of quadrature and of each other; variances of R, L, E and v_r within 3% | a sample re-weighted by R, as stage 6 drew it, must fail |
| V3 the support fits | every orbit the declared support allows is bound, lies inside the solver's radial domain and inside 0.9 of the simulator's box, and E_circ's minimising radius is continuous and increasing across the L window | the stage 6 five-width support must fail for the warm annulus |
| V4 discretization | refining the radial grid, the L and v_r quadratures and the radial domain, separately and together, moves alpha, the mean support, mean radius, radial width and peak field by < 1e-4 relative, and the density and field profiles by < 2e-3 of their peaks | the stage 6 label (support at the argmax node, unbracketed secant) must be caught moving alpha by more than 1e-4 under the same refinement |
| V5 the drawn source writes the intended field | the azimuthally averaged field of N drawn bodies approaches C0 through N = 256, 1024, 4096, 16384, sixteen realizations each: fitted exponent in [-0.7, -0.3] and rms error below 0.5% of the peak field at 16384 | the R-re-weighted sample must show an error floor above that bound |
| V6 the draw is stationary | 16,384 bodies in the FROZEN field for 30 T0, four realizations per population: each of mean radius, radial spread, rms v_r and mean v_t at t = 0 agrees with its own time average over 5-30 T0 within 4 sampling standard errors, and the instantaneous mean radius fluctuates about its mean by less than 3 standard errors rms | three, each of which must fail: every body kicked outward by 0.15; every velocity scaled by 1.05; the R-re-weighted sample |
| V7 restart | a live run saved mid-way with its full state — bodies, E and C — and resumed reproduces the uninterrupted run bit for bit | a restart from the bodies alone, with the field re-primed, must differ |
| V8 the ensemble is what it says | every declared run completes its horizon with status recorded, and each live run and its frozen partner start from byte-identical bodies | — (an identity) |

V6 replaces the construction half of stage 6's H6. It compares the ensemble at the initial instant with
its own long-time average, which a non-stationary draw cannot satisfy, where stage 6 compared two window
averages, which any draw satisfies.

## Part R. The refinement stage 5 quoted, archived

Stage 5 reported a grid and timestep refinement of the cold ring's m = 2 growth from a scratch run with no
archive and no code, and called it convergence when it was scatter. Here it is declared and archived: the
unseeded two-stage cold ring, 30 T0, m = 2 growth rate fitted over 10-30 T0, at (spacing, step) = (w/5,
0.01), (w/10, 0.01) and (w/5, 0.005).

| gate | requirement | negative control |
|---|---|---|
| R1 simulator against linear theory | each of the three growth rates within 2% of the predicted 0.018552, and each frequency within 0.5% of 2.096535 | the same measurement on a run launched with no primed field must fail |

## Part X. The experiment — a scientific outcome, archived however it falls

Identical bodies are evolved for 30 T0 with the field **live** and with it **frozen**. Populations: the
two of family A and the three of family B. Body counts 64, 256 and 1024 at fixed total mass and coupling,
with 8, 8 and 4 realizations. Recorded every quarter period: the decomposed kinematics, the support felt
by the bodies and the axisymmetric support on a ring at the initial mean radius, and the **complex**
azimuthal coefficients, m = 0 to 8, of the source S, the excitation E and the force-producing field C on
that ring — phases, not only powers.

**Dispersion means what is left after organised motion is removed**: a mean radial flow, a breathing term
linear in radius and the m = 1 and m = 2 azimuthal harmonics are fitted out by least squares, with the
variance a p-parameter fit absorbs from noise subtracted, as in stage 5. rms v_r is reported under that
name and no other.

What is reported for each quantity Q is the signed difference of the two runs' changes at common times,

    D_Q = [Q_live(late) - Q_live(early)] - [Q_frozen(late) - Q_frozen(early)],

early being the mean over the first 2 T0 and late the mean over the last 5 T0, normalised by the initial
value of Q — a declared physical scale, never another run's drift — with its standard error over
realizations. Q runs over the mean radius, the radial spread, the residual radial and azimuthal
dispersions, rms v_r, the mean and spread of the angular momenta, the support on the bodies and on the
ring, the m = 1 and m = 2 streaming amplitudes and the breathing coefficient.

**Declared readings, fixed now.**

* **A collective mode.** For m = 1 to 4 the growth rate of |C_m| is fitted over 15-30 T0, late enough that
  the shot-noise structure of a freshly live field — which builds over tau_keep and would mimic growth at
  about 0.002 — has saturated. A mode of the underlying distribution grows at a rate that does not depend
  on the body count, from an amplitude that falls as N^-1/2. It is read as **present** if the
  ensemble-mean rate exceeds 0.006 by more than three standard errors at all three body counts and the
  three rates agree within three standard errors; as **absent** if the rate plus three standard errors
  lies below 0.006 at 1024 bodies; and as **unresolved** otherwise.
* **Quiet over 30 T0**, for a quantity, only if |mean D| plus two standard errors fits inside the bound at
  1024 bodies: 1% for the mean radius, 10% for the radial spread, each residual dispersion and the
  support. Failing to detect a drift is not a demonstration of a small one; a bound the uncertainty does
  not fit inside is reported as unresolved.
* **Discreteness against collective heating.** The change in residual radial dispersion is fitted as N^p
  across the three body counts: p + 2 sigma < -0.5 is read as discreteness-dominated, p - 2 sigma > -0.5
  as collective, anything else as unresolved.
* **Temperature**, in family B only, where nothing else changes: the trend of the m = 2 growth rate and of
  D with dE at 1024 bodies, with its uncertainty.

## What this stage cannot establish

Thirty periods. A population quiet for thirty periods may hold a mode with a growth time of hundreds, and
the only thing that can exclude one is the linear mode calculation around these populations — the
transport of a perturbation along the unperturbed orbits, coupled to the memory response — which remains
declared from stage 6 and is **not run here**. Part E certifies the cold ring's strip and says nothing
about a warm annulus's spectrum. The source law is still a label and a declared mass, the dissipation
reservoir is unidentified, the kernel is instantaneous in space, and the reciprocal energy and
angular-momentum accounting has not been carried through any of these runs.

## Declared thresholds (read by the driver; the code holds none of its own)

```json
{
  "E": {"strip_re": [-0.008, 1.5], "strip_im": 6.0, "m_max": 16, "unstable_floor": 1e-7,
        "raw_known_spectra": 1e-8, "raw_robustness": 1e-6, "polished_residual": 1e-10,
        "polish_movement": 1e-4, "singular_gap": 1e3, "partition_agreement": 1e-8,
        "winding_phase_step": 0.39269908169872414, "edge_offset": 1e-6},
  "S": {"reach": 4.0, "mean_support": 0.10, "measure_identity": 1e-12, "draws": 200000,
        "mean_sigma": 4.0, "variance_relative": 0.03, "box_fraction": 0.9,
        "convergence_integrated": 1e-4, "convergence_profile": 2e-3,
        "field_exponent": [-0.7, -0.3], "field_error_16384": 0.005, "field_realizations": 16,
        "stationary_bodies": 16384, "stationary_realizations": 4, "stationary_horizon": 30.0,
        "stationary_sigma": 4.0, "stationary_fluctuation_sigma": 3.0,
        "kick": 0.15, "velocity_scale": 1.05},
  "R": {"growth_relative": 0.02, "frequency_relative": 0.005, "fit_window": [10.0, 30.0], "horizon": 30.0,
        "predicted_growth": 0.018552, "predicted_frequency": 2.096535},
  "X": {"horizon": 30.0, "record_every": 0.25, "bodies": [64, 256, 1024], "realizations": [8, 8, 4],
        "early_window": [0.0, 2.0], "late_window": [25.0, 30.0], "mode_fit_window": [15.0, 30.0],
        "mode_rate_floor": 0.006, "mode_sigma": 3.0, "quiet_sigma": 2.0,
        "quiet_bounds": {"mean_radius": 0.01, "radial_spread": 0.10, "residual_radial_dispersion": 0.10,
                         "residual_azimuthal_dispersion": 0.10, "support_on_bodies": 0.10},
        "discreteness_exponent": -0.5}
}
```
