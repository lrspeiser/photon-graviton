# RUT-1 stage 6, part 1: both repairs hold, a supported population exists in exact equilibrium, and the declared stationary-control gate failed

[protocol-rut6.md](protocol-rut6.md), declared in b5e0401 before any of this ran. Five of six gates pass.
**H6 failed as declared and is recorded as failed**; it is discussed on its own below, without being
reworded into something that passes.

Part 1 is the two repairs the owner's review of ceb0c86 named as preconditions, plus the construction of
the state part 2 will test. **No configuration is claimed stable here.**

## R1. The spectrum is now located, not counted — and stage 5's spectra turn out to have been complete

The stage 5 contour count returned a false unstable root at every mode of the neutral control. It is
replaced by a contour-integral eigensolver of Beyn's kind ([beyn.py](beyn.py)), which needs no initial
guess and returns every eigenvalue inside a declared region.

| check | requirement | result |
|---|---|---|
| H1 `exp(-s) = 1/2`, roots `ln 2 + 2πik` | located to 1e-8, count exact | 3 of 3, **1.8×10⁻¹⁷** |
| H1 Kepler control, `s²(s² + Ω²)` | **zero** roots with Re s > 0 | **0** at m = 0, 2, 5 (the old counter said 1 at every mode) |
| H1 zero-lag quartic vs its exact companion roots | 1e-8 | **1.1×10⁻¹⁵** |
| H2 doubled quadrature, shifted contour | roots move < 1e-6, no count changes | **1.8×10⁻¹⁶**, none |
| H3 against continuation, 32 writers, m ≤ 8 | every root found by both | agreed on every mode, no saturation |

**The certified result, for the declared region Re s ∈ [0.001, 1.5], |Im s| ≤ 6, m = 0…8:**

| rung | unstable modes | fastest |
|---|---|---|
| no attraction | none | — |
| instantaneous | every m from 1 to 8 | 0.155 T0 |
| one-stage | every m from 1 to 8, two roots per mode | 4.451 T0 |
| **two-stage** | **m = 2 only, one root** | **8.579 T0** |

So the two-stage cold ring has exactly one unstable mode in that region, and continuation had found it.
Stage 5's spectra were complete; the only thing broken was the counter that was supposed to certify them.
The region is a declared limit, not the whole plane: it excludes Re s < 0.001, where marginal roots sit on
the contour and a root is resolved slowly, and it stops at m = 8.

One mistake of mine on the way: my first comparison against the exact quartic reported a missed root. The
solver was right — the root sat at Re s = 1.026, outside the region's upper edge, and my comparison had
filtered the exact list on the lower bound only. Both lists are now filtered by the same region.

## R2. The seeded test works once the whole state is perturbed — exactly as the owner diagnosed

`seed_mode` perturbed positions and velocities only. The repair prescribes the mode's trajectory through a
40-period spin-up and advances the **shipped** field update through it, so the excitation and
force-producing fields carry the mode's own history, then releases the bodies.

| seed amplitude | field history | measured growth | measured frequency |
|---|---|---|---|
| prediction | — | 0.018552 | 2.096535 |
| 10⁻⁴ | built in | 0.018555 (**0.014%**) | 2.096583 (**0.002%**) |
| 10⁻³ | built in | 0.018598 (0.25%) | 2.096577 (0.002%) |
| 10⁻³ | **none** | **0.00367** | **2.372** |

Initial amplitudes scale as **10.0000** for a tenfold seed. H4 passes at 0.25% against a 5% tolerance.
The last row is the control that settles the question: with positions and velocities perturbed but the
fields left unperturbed, the run measures a different rate and a different frequency altogether. **My
stage 5 explanation — that the seeded test was intrinsically ill-conditioned by a nearby branch — was
wrong; it had a setup defect, and the owner identified it.** G4's three-line status is now: declared
seeded form failed; unseeded verification passed; full-state seeded verification **passed**.

A second mistake of mine here, caught by a zero-amplitude baseline: my first version of this gate handed
the spun-up state to `run_ring`, which builds its own field, so the bodies were launched at the supported
speed into an *empty* field. That measured a formation transient — growth 0.040, frequency 1.77 — and
seeding made no difference to it, which is what gave it away.

## R3. A supported population in exact self-consistent equilibrium

A distribution function of the conserved quantities, `f(E, L) = exp[-(L-L0)²/2dL²]
exp[-(E-E_circ(L))²/2dE²]`, solved together with the field it writes for writing proportional to mass
([equilibrium.py](equilibrium.py)). Four annuli, each solved to the declared 10% support label:

| annulus | dL | dE | FWHM | σ_r | mass | support | consistency residual |
|---|---|---|---|---|---|---|---|
| narrow, cold | 0.06 | 0.010 | 0.087 | 0.094 | 0.0268 | 0.10000 | 8×10⁻¹³ |
| narrow, warm | 0.06 | 0.030 | 0.135 | 0.167 | 0.0686 | 0.10000 | 3×10⁻¹⁵ |
| wide, cold | 0.12 | 0.010 | 0.135 | 0.092 | 0.0643 | 0.10003 | 1×10⁻¹¹ |
| wide, warm | 0.12 | 0.030 | 0.174 | 0.162 | 0.1834 | **0.09645** | 3×10⁻¹³ |

H5 passes: every annulus satisfies the coupled equations to 10⁻¹¹ or better against a 10⁻⁶ requirement, at
finite mass with edge density below 10⁻³ of peak. The wide warm annulus is an exact equilibrium at 9.6%
support rather than the requested 10%: the outer solve for the writing rate stopped short, and the achieved
value is reported rather than the label.

**And it is an equilibrium of the simulator, not just of the equations.** Sampled to 64 bodies and evolved
with the field held fixed, the narrow warm annulus changes its mean radius by 0.02%, its radial spread by
0.06% and its dispersion by 0.08% over 20 periods.

Getting there took four corrections, each of which had produced a "converged" answer that was wrong:

* A uniform velocity grid cannot resolve a cold distribution. With dE = 0.01 the energy Gaussian is
  narrower than one velocity cell, so the density was quantisation noise and the fixed point was the
  discretisation's. Integrating over L and v_r instead puts both Gaussians on the integration variables,
  and `E = E_min + v_r²/2` turns `dE/|v_r|` into `dv_r`, removing the turning-point singularity.
* A finite-difference `dΦ/dr` left a 10⁻³ noise floor that a cold distribution amplifies tenfold through
  `E_circ(L)`. The kernel's derivative is analytic — only I0e and I1e are needed — and using it removed the
  floor.
* Rescaling the writing rate *inside* the fixed point mixed two couplings and stalled. They are now
  separated: an inner solve at fixed rate, and a one-dimensional outer solve for the rate.
* Plain under-relaxed iteration is not reliably contractive for wide annuli — it converged or stalled at
  10⁻³ depending on damping and starting point. Anderson acceleration fixed that and cut the narrow solves
  from thirty seconds to three.

## H6 failed as declared

The gate: evolving the constructed equilibrium with no added perturbation for 20 periods, its drift must be
within 3× that of a control built by the same sampling.

| run (narrow warm, 64 bodies, 20 T0) | mean radius | radial spread | radial dispersion | angular momentum |
|---|---|---|---|---|
| **field written** | 0.86530 → 0.86637 | 0.04726 → 0.04719 | 0.14326 → 0.14350 | 1.00623 → 1.00617 |
| field frozen | 0.86515 → 0.86498 | 0.04727 → 0.04724 | 0.14323 → 0.14334 | unchanged |
| field decaying, no writing | 0.86732 → 0.91620 | 0.04849 → 0.07458 | 0.14058 → 0.10713 | unchanged |
| no field at all (the declared control) | 1.12112 → 1.12273 | 0.27143 → 0.27210 | 0.16358 → 0.16361 | unchanged |

Against the frozen-field control the ratios are **6.5**, 2.2 and 2.2, so the mean-radius ratio exceeds 3
and **the gate fails**. Against the control the protocol actually declared — no field — the ratios are 0.9,
0.6 and **7.0**, so it fails that way too.

Two things about this gate are my errors in declaring it, and I am stating them rather than using them to
turn the failure into a pass. **The declared control is not a control**: with no field the sampled bodies
are not the constructed population at all — they settle at mean radius 1.12 with spread 0.27 against 0.865
and 0.047 — so its drift says nothing about this state. **And a ratio of two near-zero drifts is
ill-posed**: the written run's mean radius moved by 0.12% and the frozen control's by 0.02%, both inside
what 64 bodies can resolve, and their quotient is mostly sampling noise. A well-posed version would compare
absolute drifts against a threshold fixed in advance, with enough bodies to resolve it. That gate has not
been declared, so it has not been passed.

What the same table does show, as a measurement: the third row is the control my first version of this
gate ran by accident. `write=False` lets the field **decay** rather than freezing it, which removes the
support over a retention time and pushes the population outward by 6%. Stage 3 drew exactly this
distinction between a frozen field and a decaying one, and I reintroduced the confusion here.

## The first evidence on temperature — a measurement, not a stability result

The same two runs for the narrow **cold** annulus (σ_r = 0.094 against the warm annulus's 0.167):

| narrow cold, 20 T0 | radial spread | radial dispersion | angular momentum |
|---|---|---|---|
| field written | 0.03906 → 0.04448 (**+13.9%**) | 0.09071 → 0.10341 (**+14.0%**) | 0.99769 → 0.99570 |
| field frozen | 0.03905 → 0.03908 | 0.09077 → 0.09078 | unchanged |

Under the written field the colder annulus heats by 14% in twenty periods while the warmer one changes by
0.2%, and both are stationary when the field is frozen. That is the direction one expects if velocity
dispersion suppresses the collective instability, and it is the first thing in this programme that bears
on the owner's question.

**It is not a stability result and should not be read as one.** It is one seed of 64 bodies over 20
periods with no mode analysis. Heating by discreteness — two-body relaxation mediated by the written field
— is also stronger for a colder population, and nothing here separates it from a collective mode; that
needs the body-count scan and the mode analysis of part 2. A warm annulus that is quiet for 20 periods has
not been shown to be stable, and the cold ring of stage 5 would itself have looked quiet over its first
twenty periods at this amplitude.

## What this does and does not establish

**Shown.** The spectrum of the cold-ring ladder is now located by a method validated against three spectra
known exactly, and it confirms that the two-stage ring has a single unstable mode in the declared region.
The seeded verification works once the fields are perturbed along with the bodies, to 0.014%, and fails
without that. An orbital population can be constructed in exact equilibrium with the field it writes, at
finite mass, across narrow and wide and cold and warm, and that construction is stationary under the
simulator to a few parts in 10⁴.

**Not shown.** No mode analysis has been done around any constructed population, so **no population is
claimed stable or unstable.** The certified spectra belong to the cold-ring family of stage 5, in a
declared region and for m ≤ 8. The equilibria use the stage 4/5 writing label, which is still a label and
not a source law; the dissipation reservoir is still unidentified; the kernel is still instantaneous in
space. Part 2 — the complete mode analysis around these states, targeted nonlinear tests with full-state
perturbations, the formation test, the reciprocal accounting carried through growth and saturation, and
complex mode coefficients with full restart states — is declared and not run.

**A note on the exit status.** `rut6.py` exits on whether it reproduces its archive, not on whether every
declared gate passed; a failed gate is listed in `failed_gates`, here and in the changelog, and the archive
comparison still turns the suite red if any gate's outcome ever changes. The alternative, tying the exit
status to the gates, leaves two options when one fails honestly — withhold the result, or reword the gate
until it passes — and the second is the habit the owner's review told me to stop. This is a change from
how rut4 and rut5 behave, and it is the owner's to overrule.

## Reproduce

```bash
python research_work/results/path-memory/rut6.py    # about ten minutes; also the suite job
```
