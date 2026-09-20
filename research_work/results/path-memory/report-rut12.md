# RUT-1 stage 12: formation from an empty field, and the angular-momentum ledger

Protocol `protocol-rut12.md` (48cb882), pushed with the code before any stage 12 run. Units of work
`rut12_tasks.py` and `reciprocal_budget.py`, driver `rut12.py`, archive `rut12-results.json`, series
`rut12-series/`, suite job `rut12_checks.py`. Stage 4R's reciprocal field, stage 7's population builder and
the owner's sampler are used by import, unchanged.

This is the fourth and last of the pieces of work the owner set in the review of 3a80fec.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | see `rut12_checks.py`: the identities of the ledger recomputed from scratch, the prediction recomputed, one formation run drawn and replayed, every gate and reading re-derived |
| **numerical verification** | **PASSED**, 7 gates of 7 |
| **scientific outcome** | archived as it fell; it never sets the exit status |

## What had to be derived first

Stage 4R built the reciprocal form C = W∗h, with W∗W = K, so that deposition and force come from one
interaction term, and derived an exact energy ledger. Its docstring then says what it did not do: *"it does
not finish the momentum accounting for that reservoir and the fixed centre."* Without that half, the owner's
instruction could not be carried out, only asserted.

The external potential is axisymmetric and exerts no torque, so all of the matter's torque comes from the
field; W is radially symmetric, so ∂_θ commutes with W∗, and a symmetric kernel is self-adjoint. That turns
the matter torque into a statement about the field alone, and with **J = ∫h_t ∂_θ h**:

**L_field = −(τ_f/α) J,  torque = (γ/τ_f) L_field,  d(L_m + L_field)/dt = −(γ/τ_f) L_field.**

Unlike the dissipation, which is a square and cannot change sign, the torque takes either sign: a memory that
decays gives angular momentum back as readily as it takes it.

## Part A. The ledger

| step | energy balance | angular balance |
|---|---|---|
| 0.02 | 1.26e-06 | -5.01e-09 |
| 0.01 | 3.11e-07 | -1.24e-09 |
| 0.005 | 7.76e-08 | -3.10e-10 |

Observed orders angular_balance:0.01->0.005 4.02, angular_balance:0.02->0.01 4.03, energy_balance:0.01->0.005 4.01, energy_balance:0.02->0.01 4.03 — both
second order. The angular pair at the finest step was excluded from the order requirement by the declared
floor, its residual having reached -3.1e-10, where a
ratio measures round-off rather than the leading error.

**L_field is a property of the field, not of the box it is computed in:**

| box, modes | L_field |
|---|---|
| box10.0:modes128 | 4.5861302251e-06 |
| box5.0:modes64 | 4.5861302255e-06 |
| box7.0:modes96 | 4.5861302251e-06 |

— identical to ten digits, a relative spread of 6.04e-11 against a bound of
1.0e-06. Displacing the position grid by half a box, the mistake this control
exists to catch, breaks the closure to 0.002 — 1.61e+06
times worse — and moves L_field by a factor 374. The energy ledger does not
notice it at all.

The spectral field at the declared resolution reproduces the field the population was built with to
0.0013 of its peak at worst and 4.93e-04 rms, against a bound of
0.02.

## Part G. Formation

An annulus in exact equilibrium with the bare point mass — the same f(E, L) parameters as the verified
population but α = 0 — an empty field, and writing switched on at t = 0.

With writing frozen the start holds to **8.16e-06** over the
horizon, so what the other runs do is caused by the writing and by nothing else. The verified population, in
the smooth field it writes, holds to 0.012: the thing being
compared against is one the integrator can hold.

| variant | movement of the settled structure |  |
|---|---|---|
| B24:bare:s0:big_box | 1.39e-09 | bounded at 0.05 |
| B24:bare:s0:fine_modes | 8.68e-15 | bounded at 0.05 |
| B24:bare:s0:half_step | 2.08e-05 | bounded at 0.05 |
| B24:bare:s0:long_horizon | 0.055 | reported, not bounded |
| B24:bare:s0:more_bodies | 0.028 | bounded at 0.05 |

Across every run the worst closure at any recorded instant was
7.65e-07 in energy and
7.35e-08 in angular
momentum — **through the transition, not only at its ends.**

## What formed

Settled over the records from period 30 to the horizon, 4
realizations, against the prediction fixed in the protocol before any run:

| quantity | formed | predicted | target | formed / predicted | formed / target |
|---|---|---|---|---|---|
| mean L | 1.0115 ± 6.2e-04 | 1.0109 | 0.99844 | 1.0005 | 1.0130 |
| sigma L | 0.060099 ± 4.3e-04 | 0.060665 | 0.059857 | 0.9907 | 1.0040 |
| mean rc | 0.98814 ± 0.002 | 0.9414 | 0.93492 | 1.0497 | 1.0569 |
| sigma rc | 0.041456 ± 4.5e-04 | 0.031397 | 0.031349 | 1.3204 | 1.3224 |
| rms epicycle | 0.15359 ± 0.001 | 0.14432 | 0.091231 | 1.0642 | 1.6836 |
| support | 0.10849 ± 4.3e-04 | — | 0.10257 | — | 1.0577 |

**Angular momentum comes through formation intact.** ⟨L⟩ lands at 1.0005 of the prediction and
σ_L at 0.9907 — the distribution the bodies started with, carried across a transition in which the field went
from nothing to its full depth. That is what the ledger above guarantees, seen rather than assumed.

**What forms is not the verified population: it is a hotter, wider relative of it.** Its radial motion is
1.6836 times the target's, the spread of its circular radii 1.3224 times, and the whole annulus sits
1.0569 times further out while carrying 1.0577 times the support. Adiabatic invariance predicted the heating
and got its size right to 6%: 1.0642 of the prediction against 1.6836 of the target. So the answer to
whether the acceptable state can form naturally is **a supported annulus forms from nothing, with the right
angular-momentum structure, and it is not the verified state** — whether that counts as the acceptable state
forming is the owner's judgement, and it is not a gate.

**Where the prediction failed, and why.** It put the circular radii at 0.9527 of what formed, because it read
the start's angular momenta in the TARGET's field. The formed state writes its own, and being hotter and wider
it writes a different one, so its bodies settle further out. The prediction had no self-consistent step; the
measurement shows exactly that gap, and σ_rc misses by the same cause, 0.7573 of what formed, while
σ_L is right to 0.9907. The adiabatic argument is right about what is conserved and silent about what the
conserved thing is embedded in.

## The budgets through the transition

| realization | energy dissipated | torque integral | L_field at the horizon | energy residual | angular residual |
|---|---|---|---|---|---|
| s0 | 0.0517 | 3.48e-05 | -1.38e-05 | -3.3e-07 | 1.1e-08 |
| s1 | 0.051 | 1.92e-04 | 3.86e-05 | -2.5e-07 | -3.2e-08 |
| s2 | 0.0522 | 2.37e-04 | 1.15e-04 | -4.3e-07 | -4.9e-08 |
| s3 | 0.0521 | 3.20e-05 | 7.52e-06 | -9.1e-08 | 7.3e-08 |

The energy dissipated over the transition is the same to 2% across realizations, and it is not small: about
0.0518 against an initial ⟨L⟩ of order one and a support that reaches 0.108. The angular momentum that
passed through the field is four orders smaller and **changes sign between realizations**, which is what a
nearly axisymmetric system should do — the field takes angular momentum from one part of the annulus and
returns it to another rather than draining it. Neither is a fitted quantity: both are what the derived rates
integrate to, and the residuals in the last two columns are what is left over.

## What I got wrong

Three, all of them in this stage’s own construction and all caught before the campaign, which is why they
appear above as controls rather than as corrections.

The position grid. `ifft2` returns the field on [0, L), so a moment integral needs that coordinate wrapped
into [−L/2, L/2). I had it displaced by half a box. L_field then varied threefold with the box and the
closure was destroyed — and **the energy ledger never noticed**, because it is computed by Parseval and never
touches a position. I found it only because L_field was not invariant, which is now gate A2.

The mature control. Priming it from instantaneous point sources puts every body at the bottom of its own
Gaussian well; the field is as lumpy as the sample, it relaxes violently, and the energy balance reached 4.6.
What failed there was the integrator, not the ledger. It is now primed from the population’s own smooth
σ(r), and gate A3 requires that field to be the one the population was built with.

The identity ∫h ∂_θ h = 0. It is ½∫∂_θ(h²), which vanishes on a rotationally closed domain — and a
periodic square is not one. My first test of it used random noise filling the box, where it does not vanish
at all; on a localised field it vanishes to 2.7×10⁻¹⁷. The suite job now tests it the right way, and A2
gates the same approximation independently through box invariance.

And one thing the prediction got wrong rather than the code: it had no self-consistent step, as the reading
above says.

## What this stage does not establish

The reservoir. The ledgers say where the energy and the angular momentum go — into the field, and out of it
at a derived rate — but nothing here identifies what physically receives them. A decaying memory remains a
postulate rather than a mechanism, and that is the honest limit of "without an unaccounted source or sink":
the accounting closes, the destination is still named rather than found.

Spatial causality: the field responds instantaneously in space. The fixed centre: the central mass is held,
so it absorbs linear momentum silently, and the linear-momentum budget is not attempted. The horizon: forty
periods, with one run at eighty that moved the settled structure by 5.5%, so the settled state is settled only
on the declared horizon — and stage 11’s m = 1 mode has an e-folding time of 181 periods, far beyond either,
so nothing here speaks to it. One population: B24 alone. And the formed state was not itself put through the
mode calculation, so whether the hotter annulus that forms is one this family would call supported is a
question this stage raises and does not answer.
