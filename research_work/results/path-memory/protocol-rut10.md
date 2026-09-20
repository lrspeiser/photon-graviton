# RUT-1 stage 10: the seeded, full-state warm-mode experiment

Declared before any stage 10 run. Its own protocol file. Stages 4 to 9 stay frozen: their archives keep
reproducing, stage 6's H6 stays failed, stage 8's L3 and L3b stay failed in stage 8's archive, and no file of
theirs is edited.

## Why this stage exists

B13 — family B at dE = 0.013, never simulated before its forecast was pushed — grew an m = 2 pattern at the
forecast pattern speed and 22% more slowly than forecast, at 5.4 standard errors. The owner ruled that this is
the central scientific discrepancy of stage 8, and asked for four things. Stage 9 ran the three that need no
new simulation of the live model. Its verification passes, and the answer is no three times over: the boundary terms of the sharp
truncation move B13's forecast rate by +0.006%; no refinement of the representation moves its root by more than
2.1×10⁻⁵ of |s|; and B13's own eigenfunction, imposed at its predicted complex rate on 1.6 million bodies of the
simulated population, regenerates itself with an implied shift of the growth rate of −0.30 ± 0.24%, where a
uniform change of the response of −23% would be needed for the shortfall. The 22% is not an error of the linear
calculation. This stage is the fourth:

> the seeded test needs to initialize the warm population's perturbation, not merely reuse the cold-ring
> particle displacement … construct the corresponding perturbation of the distribution function and the
> excitation field. For a mode proportional to exp(st), delta_E = (s + 1/tau_keep) delta_C … Use multiple
> perturbation amplitudes, paired signs where useful, several sampling realizations and numerical refinement.
> The seed must be above the measured noise floor but below the nonlinear regime. A quieter sampling
> construction may be more efficient than increasing body count indiscriminately.

and the question it is to answer: **does the measured rate approach the predicted rate as the seed amplitude
and the numerical error decrease?** No outcome is selected in advance.

## The construction

**A quiet sample.** N/4 bodies are drawn with the owner's sampler, the sampler of record, and each is joined
by its images under rotation by 90°, 180° and 270°, which are exact in floating point. The population is
axisymmetric, so this is a sample of it; it has four-fold symmetry, so its m = 1, 2 and 3 harmonics vanish
identically. The simulator's field lives on a Cartesian grid that is symmetric under the same rotations, and
its equations are covariant, so **the live system keeps the symmetry**: an m = 2 pattern can then grow only
from what is put there, or from rounding. It is a quieter sampling of the azimuth and of nothing else: the
sample still holds only N/4 independent orbits, and that is what limits it (below).

**The seed is the predicted eigenmode in the full state**, from stage 9's corrected calculation — root s and
field δC_m(r) = Σ_b c_b k_m(r_b, r), scaled to a largest value of 0.1:

* *the bodies* carry δf because they have been **moved by the mode's own field**: the quiet sample is evolved in
  the frozen equilibrium potential plus ε e^{γt} Re[δC_m(r) e^{i(mθ+βt)}] from eight e-foldings before t = 0.
  That is the first-order response δf_l = F_l ψ_l/(ν_l − i s) of *this* sample of *this* sharply truncated
  population, boundary terms included because the bodies are the population, with nothing of the cold ring's
  displacement in it. It is the operation stage 9's gate E1 verified to regenerate the mode's field;
* *the force-producing field* is the constructed C₀ plus ε Re[δC_m(r) e^{imθ}], value and analytic gradient, on
  the simulator's grid;
* *the excitation* is C₀/τ_keep plus ε Re[(s + 1/τ_keep) δC_m(r) e^{imθ}]: for a mode growing as e^{st},
  dC/dt = E − C/τ_keep requires δE = (s + 1/τ_keep) δC, and the equilibrium relation δE = δC/τ_keep would be wrong
  by the factor 1 + sτ_keep = 3.0 − 5.5i. That mistake is kept as a negative control.

ε is set so that the seeded |C₂|/C₀ on the ring at the population's own mean radius is the declared target. The
seeded field does not depend on the realization and is **saved whole**; every run's initial bodies are **saved
whole**; the mode's complex coefficients, its nodes and its root are saved with every run; and the suite job
**replays the exact preparation** of one run — sample, preparation, field — and compares bodies and field bit
for bit before replaying a prefix of the run itself.

## What was run before this declaration, and what it showed

A scratch prototype, not archived. **I have therefore seen results of the kind this stage will produce, and
every threshold below was set with them in view.** They are stated so that the thresholds can be judged:

* **the quiet sample stays quiet, and does not heat.** B13, 1,024 bodies, unseeded, live for 20 T0: |C₂|/C₀ never
  above 7×10⁻¹⁷, |C₁|/C₀ never above 2.4×10⁻¹⁶; the m = 4 noise it does carry is 1 to 6×10⁻⁵; its radial
  dispersion wanders between 0.097 and 0.107 as the frozen runs of stage 8 do, with no trend. B16, 40 T0:
  |C₂|/C₀ never above 2.3×10⁻¹⁵. An ordinary draw of the same size reaches 2.3×10⁻⁵ in one period and 2.4×10⁻⁴
  in six (stage 8's archive), and the pattern saturates near 2×10⁻³: **an unseeded run has under two e-foldings
  between its own noise and saturation, and a seeded quiet run has seven;**
* **the field is seeded as declared**: a target of 10⁻⁶ gives |C₂|/C₀ = 1.000×10⁻⁶ at t = 0, and E₂/C₂ =
  0.048078 − 0.087867i against s + 1/τ_keep = 0.048077 − 0.087866i;
* **the source written by the prepared bodies is too noisy at this body count to be a gate**: S₂/E₂ = 2.19 − 2.23i,
  1.51 − 2.08i and 2.84 − 2.62i in three realizations against 1 + sτ_form = 1.61 − 1.66i — the phase is right and
  the modulus scatters by 40%, which is what 256 independent orbits give (stage 9 measured the same response to
  1% with 200,000). Unprepared bodies write S₂ = 0 to rounding. The first four periods of every seeded run are a
  transient for this reason, and the declared window starts after them;
* **B13 seeded at 10⁻⁶, 1,024 bodies, three sampling realizations: growth rates 0.03291, 0.03093 and 0.03162**
  (mean 0.0318 ± 0.0006) **against the forecast 0.03216, and pattern speeds +0.04542, +0.04588 and +0.04523**
  (0.04551 ± 0.0002) **against +0.04393.** The unseeded reading of stage 8 was 0.0252 ± 0.0013. The seeded rate
  agrees with the forecast; the pattern speed is 3.6% high in every realization
  — **and it drifts**: started at 3, 10 and 30 times the seed the same fit gives rates of 0.03182, 0.03184 and
  0.03187 (each ± 0.0006) and pattern speeds of +0.04551, +0.04530 and +0.04490 (± 0.0002 to 0.00005). The 10⁻⁵
  run shows the same pattern speed at the same TIME, not at the same amplitude, so it is a slow transient or a
  slow change of the background and not an amplitude effect; I do not know which. In the first realization the
  local rate over successive
  four-period windows is 0.0390 (the transient), then 0.0338, 0.0334, 0.0332, 0.0321, 0.0331, 0.0333 — a
  wander of ±3% in one realization — and then 0.0355 and 0.0302 once |C₂|/C₀ passes 4×10⁻⁴;
* **linear in the seed**: the 10⁻⁵ run's amplitude over the 10⁻⁶ run's keeps its initial value to 0.19% while
  the larger is under 2×10⁻⁴, and has fallen by 30% by the time it reaches 2.5×10⁻³: hence the ceiling. **The −10⁻⁶
  run is not independent information**: a sign flip of an m = 2 seed is a rotation by 90°, the quiet sample is
  symmetric under it, and the two runs agree to 2×10⁻¹¹ in amplitude and 3×10⁻¹¹ in phase. It is kept as a test
  that the simulator is rotation-covariant, which the quiet start relies on, and not as a paired measurement;
* **the grid does not matter and neither does the step**: on a w/8 grid the rate is 0.03292 against 0.03291
  and the pattern speed +0.04542 against +0.04542. Halving h_max alone changes nothing because it is not the
  binding limit — the step is min(h_max, η√(r³/GM)) and η = 0.01 gives about 0.005 — so the declared refinement
  halves both: with both halved the rate is 0.03291 against 0.03291 and the pattern speed +0.04542 against
  +0.04542. Steps four and eight times COARSER still give 0.03292 and 0.03293; a grid of w/2.5 gives 0.03276
  (−0.5%), one of w/1.25 gives 0.03004 (−8.7%) and one of w gives 0.0162. The measurement is far inside its
  numerical margins, the w/1.25 grid is the control that G2 must catch, and **the 3.6% in the pattern speed is
  not the simulator's step or grid**;
* **body count**: one realization at 4,096 bodies gave 0.03027 and +0.04539 — a steadier exponential, with local rates
  between 0.0297 and 0.0318 from the first window on and no visible transient, 5.9% below the forecast, and the
  same 3% in the pattern speed, which is therefore not a small-number effect either. One realization is not a
  mean: every sample is a slightly different population, and the spread between realizations, not the step, the
  grid or the seed, is what limits this measurement. Hence eight realizations at 1,024 bodies and four at 4,096;
* **the forecast is a sensitive function of the background**: with the distribution function held, scaling the
  equilibrium field C₀ by +0.1% moves B13's predicted rate by +0.67% and its pattern speed by −0.42%; +1% in
  dE moves them by +0.06% and +0.37%; +1% in dL by +0.47% and +0.16%. In two prototype runs the live field on
  the ring fell slowly — by 0.03% in the first eight periods and 0.25% in the last eight of 32 in one, 0.01% and
  0.09% in the other — as the sample's
  mean radius crept outward by 0.3%. That is worth about −1% in rate and +0.5% to +1% in pattern speed by the
  proxy of a uniform scaling: some of the offsets above, not all of the 3.6%. In stage 8's ordinary runs the
  same field stays within 0.01% of its primed value for ten periods;
* **B16, nearer the branch's end, cannot be read at 1,024 bodies**: two realizations gave 0.0216 and 0.0129 against
  a forecast of 0.0202, with local rates wandering between 0.005 and 0.027. A sample of 256 independent orbits is
  a slightly different population each time, and near zero growth that matters more than the mode. B16 is
  therefore a reading at 4,096 bodies only, and it carries a power requirement (below) so that a wide error bar
  cannot pass for agreement, which is how stage 8's B16 was "confirmed";
* **seen in stage 8's archive, exploratory**: B24, which holds no unstable root, carries an m = 2 shot-noise
  pattern of 1 to 3×10⁻⁴ at 1,024 bodies. B13's unseeded |C₂|/C₀ is 2.8×10⁻⁴ at 6 T0, where stage 8's fit window
  starts, and 9.6×10⁻⁴ where it ends. Subtracting B24's power from B13's at those two times alone gives a rate
  of about 0.035. That is two points and a proxy; the declared diagnostic below makes the full version of it once.

## Three statuses, kept separate (the owner's ruling)

| status | question | consequence of failure |
|---|---|---|
| reproduction | does the run reproduce its frozen archive, preparation included? | the job exits non-zero |
| numerical verification | is the sample quiet and the population's own; is the field seeded as declared; is the growth linear in the seed and insensitive to the step and the grid; did every run complete? | blocks every dependent conclusion; the job exits non-zero |
| scientific outcome | the seeded growth rate and pattern speed against the forecast of record | archived as it falls; never affects the exit status |

## Part Q. The quiet sample

| gate | requirement | negative control the gate must reject |
|---|---|---|
| Q1 the sample | for every quiet sample used: the mean radius, the mean angular momentum and the radial velocity variance agree with the owner's node distribution within 4 standard errors of the N/4 INDEPENDENT bodies; the source the bodies would write on the ring has \|S_m\|/S₀ under 10⁻¹² for m = 1, 2, 3 | an ordinary draw of the same size must be seen to carry an m = 2 source above 10⁻³ |
| Q2 the symmetry survives the live run | an unseeded quiet run of B13 over the full horizon completes with \|C_m\|/C₀ under 10⁻¹⁰ for m = 1, 2, 3 at every record | an unseeded ORDINARY draw must be seen to grow \|C₂\|/C₀ above 10⁻⁵ from its own noise within 8 T0 |

Whether the unseeded quiet population heats is a reading, with stage 7's decomposition, and carries no gate.

## Part P. The preparation

| gate | requirement | negative control the gate must reject |
|---|---|---|
| P1 the seeded field | at t = 0 the ring's \|C₂\|/C₀ is the declared target to 2%, and E₂/C₂ = s + 1/τ_keep to 10⁻³ | a field seeded with the equilibrium relation δE = δC/τ_keep must be seen to miss that ratio by more than 50% |

The source written by the prepared bodies is *not* a gate, for the reason in the prototype: at 1,024 bodies it
scatters by 40%. What verifies the bodies' δf is stage 9's E1. **A second wrong preparation is run and read, not
gated: the declared field seeded on bodies that were not prepared**, which is the mistake the owner warned of
in its simplest form. Both wrong preparations are followed for the full horizon, so that what each does to the
measured rate is on record.

## Part G. The growth of the seeded mode

B13: 1,024 bodies (eight sampling realizations) and 4,096 bodies (four), seeded at +10⁻⁶, 40 T0; for the first
realization at 1,024 also +10⁻⁵ and −10⁻⁶, the step halved, a field grid of w/8 in place of w/5, a grid of w/1.25
as that gate's control, and the two wrong preparations.

**The window, declared here.** For each run, from the first record at which |C₂|/C₀ exceeds three times its
value at t = 0 to the last record before it first exceeds 2×10⁻⁴; at least 8 T0 long. The growth rate is the
slope of log |C₂|/C₀ in it and the pattern speed is minus half the slope of the phase of C₂.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| G1 linear in the seed | read on the SAME records of two runs, so that one realization's wandering cancels: while the +10⁻⁵ run is under the ceiling, the ratio of its amplitude to the +10⁻⁶ run's keeps its initial value to 1%; and the −10⁻⁶ run has the +10⁻⁶ run's amplitude to 10⁻³ and a phase π away to 10⁻³ | the same ratio followed ABOVE the ceiling, to 2.5×10⁻³, must be seen to fall by more than 10% |
| G2 the step and the grid | halving the step — h_max and η together, since η is the binding limit — and refining the grid to w/8 each change the first realization's rate by under 0.5% | a grid of w/1.25, coarse enough to matter, must be seen to change it by more than 2% |
| G3 | every run completes (a completeness condition, not a numerical gate) | — |

## The declared reading

Against **the forecast of record — stage 8's, pushed in aaf0bfe before anything was simulated** (B13: 0.03216
and +0.04393; B16: 0.02023 and +0.04897). Stage 9's recalculated values are reported beside them. At the largest
body count of each population:

* *underpowered*, and no more is said: the standard error of the mean rate exceeds 5% of the forecast. A wide
  error bar is not an agreement;
* *reproduced quantitatively*: the mean rate is within 5% of the forecast or within three standard errors of
  it, **and** the mean pattern speed is within 5% or three standard errors;
* *contradicted*: the mean rate differs from the forecast by more than three standard errors **and** more than
  15%;
* *unresolved*: anything else.

**B16 is read at 4,096 bodies only (four realizations, 55 T0) and carries no gate**, for the reason in the
prototype. And, as a table with no reading attached: the rate and the pattern speed against seed amplitude,
step, grid, body count and preparation — the owner's question.

## Four declared diagnostics

No new runs, and no verdict; declared so that each is made once, as written.

0. **Where the window starts.** The same fit started at 3, 10 and 30 times the seed, for every seeded
   population and body count: in the prototype the rate did not move and the pattern speed drifted towards the
   forecast (below), and a reader should see which.

1. **Stage 8's own growth-phase rule applied to this stage's seeded quiet runs**: is that rule biased when
   there is no shot noise to start from?
2. **Stage 8's archived unseeded runs read again with a noise floor.** At each body count the ensemble-mean
   m = 2 power of B24 — the root-free population, whose m = 2 pattern is shot noise and nothing else — is
   subtracted from that of B13 and of B16, and the rate is refitted over stage 8's own window. B24's noise is a
   *proxy* for theirs: a population nearer instability amplifies its noise more. The result is a rough
   indication of how much of the 22% is the noise under the measurement, and nothing finer.

3. **The live background.** For each seeded run, the axisymmetric field on the ring averaged over the window,
   against its primed value; and what the predicted root's sensitivity to a uniformly scaled C₀ — a declared
   finite difference, ±10⁻³ — makes of it. The forecast is a sensitive function of the background (prototype,
   above), so a run whose field has drifted is not testing quite the forecast population. A uniform scaling is
   a *proxy* for whatever the live field actually does, and the reading is rough.

## What this stage cannot establish

Stability of anything. A mechanism for the unseeded reading beyond the rough diagnostic above. Anything about
the nonlinear saturation. Modes other than m = 2, or slower than the declared horizon can show. The rate of a
population as close to the end of the branch as B16, unless 4,096 bodies prove enough. The source law is still
a label with a declared mass, no energy or angular-momentum budget is carried, and no observational comparison
is made.

## Declared thresholds (read by the driver; the code holds none of its own)

```json
{
  "Q": {
        "moment_sigma": 4.0,
        "symmetric_harmonic_max": 1e-12,
        "control_harmonic_min": 0.001,
        "live_floor_max": 1e-10,
        "live_control_min": 1e-05
  },
  "P": {
        "peak_field": 0.1,
        "efoldings": 8.0,
        "step": 0.01,
        "amplitude_relative": 0.02,
        "excitation_relative": 0.001,
        "excitation_control_min": 0.5
  },
  "G": {
        "populations": ["B13"],
        "quarters": [256, 1024],
        "realizations": [8, 4],
        "targets": [1e-06, 1e-05],
        "horizon": {"B13": 40.0},
        "unseeded_ordinary_horizon": 8.0,
        "record_every": 0.25,
        "h_max": 0.01,
        "eta": 0.01,
        "refinements": ["half_step", "fine_grid"],
        "window_lo_factor": 3.0,
        "window_ceiling": 0.0002,
        "window_shortest": 8.0,
        "control_window": [0.001, 0.0025],
        "amplitude_agreement": 0.01,
        "sign_agreement": 0.001,
        "sign_phase_agreement": 0.001,
        "ceiling_control_min": 0.1,
        "refinement_agreement": 0.005,
        "noise_proxy": "B24",
        "reading_only": {"B16": {"quarter": 1024, "realizations": 4, "horizon": 55.0}},
        "quiet_windows": [5.0, 5.0],
        "window_lo_factors_reported": [3.0, 10.0, 30.0],
        "background_delta": 0.001,
        "grids_per_w": {"fine_grid": 8.0, "coarse_grid": 1.25},
        "refinement_control": "coarse_grid",
        "refinement_control_min": 0.02
  },
  "D": {
        "confirm_sigma": 3.0,
        "confirm_relative": 0.05,
        "confirm_pattern_relative": 0.05,
        "contradict_relative": 0.15,
        "power_relative": 0.05
  },
  "C": {
        "replay_prefix": 2.0
  }
}
```
