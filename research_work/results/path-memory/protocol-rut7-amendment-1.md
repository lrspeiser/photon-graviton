# RUT-1 stage 7, amendment 1: two negative controls that could not be met as worded, and the definitions the protocol left open

Declared before the stage 7 campaign is run, in its own file. `protocol-rut7.md` (e59ef6f) is not edited
and its thresholds are not changed: **this amendment introduces no number**. The driver still reads every
tolerance from the protocol's machine-readable block.

## What was run before this amendment, and what it showed

Building the driver meant running its pieces on the real populations. No part X run of a declared size has
been made, and nothing has been archived. What has been seen, because the seeds are declared and the
campaign will reproduce these numbers exactly:

* the five populations (alpha = 0.03434181, 0.02943086 for family A; 0.03156101 for family B, whose mean
  supports come out at 9.29%, 10.00% and 10.61%);
* gates E1, E2, E4 and the E3 control, and the spectrum at a handful of m; V1, V2, V3, V4, V5 and V7 on
  all five populations; V6 for one realization of the declared draw and of the R-re-weighted control, for
  `A_cold` and `A_warm`; and the four part R runs;
* two smoke campaigns at toy sizes (16 to 64 bodies for 4 T0, 1,024 bodies for the stationarity check),
  which exercise every code path and are never archived, and 64-body timing and restart tests of 4 T0.

Three things did not pass as first coded, and they are of two different kinds.

**One was my estimator, and it is fixed in the code, not here.** V4 failed on three of the five
populations because I read the "peak field" as the largest node value, which moved by 1.3e-4 to 2.0e-4
under radial refinement against a declared 1e-4 — while the field profile itself agreed to 1e-5 of its
peak everywhere. The largest node is not the peak: sampling a smooth maximum on a grid costs C'' dr^2/8,
about 3e-4 here. It is the mistake this stage records against stage 6's label, made again by me. With the
peak taken as the vertex of the parabola through the three nodes around it, the same refinement moves it
by 1.4e-7 and V4 passes on all five (worst integrated quantity 3.2e-5, worst profile 1.0e-5). The
definition is written down below so it cannot drift.

**Two are negative controls that cannot pass as worded for a cold annulus**, for a reason that is itself
worth recording: the stage 6 sampler's defect shifts the mean radius by Var(R)/Mean(R), which is 0.0087
for the warm annulus and only 0.0025 for the cold one.

| gate | the control as worded | warm annulus | cold annulus |
|---|---|---|---|
| V5 | the R-re-weighted sample "must show an error floor above that bound" (0.5% of the peak field at 16,384 bodies) | 0.96%: above | **0.31%: below** |
| V6 | the R-re-weighted sample "must fail" a per-realization test at 4 standard errors | z = 8.6: fails, as it must | **z = 3.1: passes** |

In both cases the requirement half of the gate passed, and the broken sample *was* distinguishable — the
controls were worded around the warm prototype they were sized on. Left alone, the protocol's own rule ("a
gate whose negative control passes is itself recorded as failed") would record V5 and V6 as failed for
every cold population because of how I worded a control, not because of anything the sampler did. That
would be an accurate record of a badly specified check and a useless one. So the wording is corrected
here, before the campaign, with the thresholds untouched, and both the original and the amended
evaluations are archived so that the difference stays visible.

## Amendments

**A1. V5's negative control goes through the same gate as the sample it imitates.** The R-re-weighted
sample is drawn at the same four body counts with the same sixteen realizations, and the same two
declared criteria are applied to it: the fitted exponent in [-0.7, -0.3] and the rms error at 16,384
bodies below 0.5% of the peak field. The control is rejected if it fails either. A biased sample has an
error floor, so its error stops falling with N and its exponent collapses toward zero; that is what a
"floor" means and it does not depend on how large the floor is. Archived for every population: the
control's exponent, its error at 16,384, and whether that error exceeds the bound — the original wording.

**A2. V6 is evaluated per realization as declared, and also on the ensemble of four.** For each of the
four statistics the mean over realizations of z, times the square root of the number of realizations,
must lie within the same 4 standard errors. This makes the requirement **stricter**, with no new
threshold. Each negative control is run for four realizations as well, and is rejected if any realization
fails or the ensemble fails. A re-weighted cold annulus relaxes by about three standard errors in every
realization, always in the same direction: a per-realization test cannot see that and an ensemble test
can.

**A3. V8 covers the declared runs.** A negative control may end early — a kicked population that leaves
the box has been caught — so its status is recorded and is not required to be `completed`.

**A4. D_Q for a quantity that starts at zero.** The protocol normalises every D_Q "by the initial value
of Q". The m = 1 and m = 2 streaming amplitudes and the breathing coefficient start at the noise level,
and dividing by that is dividing by another run's noise, which the protocol itself forbids. They are
normalised instead by the initial residual radial dispersion — the breathing coefficient first
multiplied by the initial radial spread, to make it a velocity. The other nine quantities are normalised
by their own initial values, as declared.

## Definitions the protocol left open, fixed here

* **E1, Kepler control.** "Located" means Beyn's count with multiplicity is four, it equals the winding
  number, and each root lies within the declared resolution floor (`unstable_floor`) of 0, 0, +/- i Omega.
* **E3's negative control** is satisfied if the small-residue root is located, or if the disagreement
  between Beyn's count and the winding number reports it. Either way it is not silently dropped.
* **E4** places the root `edge_offset` inside each of the four edges in turn and requires it found at
  all four. Its control is rejected if the stage 6 solver — filter, then polish — loses the root at any.
* **The edge rule** is applied to located roots inside a sub-rectangle, and to polished roots that fall
  just outside an outer edge of the strip; in both cases the edge moves 0.002 away from the root and the
  whole strip is swept again. Roots still within 0.001 of an edge after that are listed.
* **The R-re-weighted control**, wherever it appears (V2, V5, V6), is a correct draw of four times the
  bodies resampled without replacement with probability proportional to R — exactly what multiplying the
  node weights by 2 pi R did. In V2 it is a quarter of the 200,000 draws, and a body's energy is its
  kinetic energy plus a cubic spline of the solver's potential at its radius.
* **V3.** "Inside" is strict: an orbit that reaches the first or last radial node does not fit. "E_circ's
  minimising radius is continuous" is tested by halving the spacing in L: a continuous function's largest
  step halves, a switch of the global minimum between two wells keeps its size, and the rule is that the
  largest step on 321 points is under three quarters of the largest on 161. The control — the five-width
  support — is evaluated on the populations with dE = 0.030.
* **V4.** The refinements are: radial nodes 280 to 559 (the same nodes and one between each pair); L nodes
  160 to 320; v_r nodes 96 to 192; radial domain [0.3, 3.0] to [0.2, 3.3] at the same spacing; and all
  four together. For the fixed-physics populations alpha is imposed, so it is the mean support that is
  followed. **The peak field is the vertex of the parabola through the largest node and its neighbours**,
  not the largest node value: the latter carries a sampling error of order C'' dr^2/8, about 3e-4 of the
  peak, which is larger than the declared 1e-4 and is the same mistake this stage records against stage
  6's label. The control is stage 6's own solve at 280 and 559 nodes, dL = 0.06, dE = 0.010 and 0.030.
* **V5's error norm** is the rms over the solver's radial grid of the difference from C0, over the peak
  of C0, combined over realizations in quadrature.
* **V6.** The kick is an outward radial velocity of 0.15 added to every body. Sampling standard errors
  come from the quadrature moments: sqrt(Var/N) for the two means, and sqrt[(mu_4 - sigma^4)/(4 sigma^2 N)]
  for the radial spread and for rms v_r.
* **V7** uses 64 bodies of `B_mid`, saved at 2 T0 and compared at 4 T0.
* **Part R** sets the integrator's accuracy parameter equal to the declared step, so that halving the
  step halves the step actually taken.
* **Part X.** Seeds are the first 32 bits of the SHA-256 of a label naming the population, the body
  count and the realization. Rates are alpha times the mass over the body count. A frozen run's field
  coefficients never change, so only its source coefficients are archived; a live run also archives the
  coefficients of the inward pull. Series are rounded to eight significant figures and every gate and
  reading that uses them is computed from the rounded, archived files. The pattern speed of mode m is
  minus the phase rate of its coefficient over m. The discreteness exponent is fitted only if the change
  is resolved from zero at two standard errors, with one sign, at all three body counts; otherwise that
  reading is `unresolved`. The quiet reading has a third value, `exceeds the bound`, when |mean D| minus
  two standard errors lies outside it.

## What this amendment does not do

It does not change a threshold, a population, a body count, a horizon, a window or a declared reading,
and it does not touch the scientific outcome: part X has not been run.
