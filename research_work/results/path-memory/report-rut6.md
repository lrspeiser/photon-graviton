# RUT-1 stage 6, part 1: both repairs hold, a supported population exists in exact equilibrium, and the declared stationary-control gate failed

[protocol-rut6.md](protocol-rut6.md), declared in b5e0401 before any of this ran.

> ## Corrections after the owner's review of 72aef44 and an audit of this stage
>
> **Status, in the owner's three separate senses.** *Reproduction:* the archive reproduces. *Numerical
> verification:* H1, H3, H4's rate and frequency, and H5 passed their checks **as coded**; H2 does **not**
> meet its declared tolerance once it is evaluated on what the contour solver located rather than on
> Newton's fixed point; several checks could not have failed. *Scientific outcome:* H6 failed as declared;
> **annulus stability is not established, and both archived annulus runs carry a growing m = 2 mode.**
>
> The owner found the first of these. I then had the stage audited, five lenses with adversarial
> verification, and re-verified the consequential findings myself. **The repaired cold-ring eigenmode
> results — the m = 2 root and the full-state seeded test — stand. Most of what this report said about the
> annuli does not.** Nothing below changes an archived number; it changes what the numbers mean.
>
> 1. **The sampler draws a different population from the one the solver constructs** (owner). In
>    (R, L, v_r) the measure is `f dR dL dv_r dtheta` with no factor of R; `sample()` multiplies its node
>    weights by `2 pi R`, so it draws p(R) proportional to R times the intended distribution. Confirmed: the
>    sampler's exact node marginal gives mean radii 0.89671 and 0.86673 against an intended 0.89498 and
>    0.86364, the owner's analytic estimates to five digits. Every annulus *simulation* here evolved a
>    mis-weighted population.
> 2. **"The first evidence on temperature" is withdrawn.** Every statistic in `_drift` is an azimuthal
>    average, second order in any mode amplitude, so it could not see what was happening. Replaying the
>    archived runs bit for bit and recording the field's azimuthal harmonics: the **warm** run's m = 2
>    component of the inward pull grows from 5×10⁻⁴ at 5 T0 to 6×10⁻³ at 25 T0, **e-folding 8.2 to 9.1 T0,
>    inertially stationary** — the signature of the cold ring's mode — and the **cold** run's grows with
>    e-folding 5.3 to 5.7 T0. The cold run's "+14% dispersion" is that mode: rms v_r rises 12.8% while the
>    **residual dispersion changes by −0.6%** (frozen control −1.0%) and the m = 2 streaming rms triples.
>    It is the mistake stage 5 had already corrected once. The contrast is not "14% against 0.2%"; it is an
>    m = 2 growth rate near 0.030 against 0.018, still confounded by item 1 and by the next item.
> 3. **The cold/warm comparison was never temperature-only** (owner). The archived narrow populations have
>    mass 0.026785 and 0.068643 and alpha 1.77354 and 0.99890: the warm one carries 44% more total writing.
> 4. **"10% support" is not what these populations feel, and the label is not even well defined.** It was
>    evaluated at the argmax *grid node* of the density, beside the field's maximum where the gradient
>    passes through zero and changes by 0.05 to 0.065 per node — half the target. `support_at(alpha)` is
>    therefore a sawtooth with several roots (narrow cold: alpha = 1.37, 1.56, 1.77), the secant returned
>    whichever its path visited, and the archived alpha moves 10 to 35% with the radial grid. "Support
>    0.10000" was false precision. The **mass-weighted** support is 13.4% (cold) and 20.3% (warm), with the
>    memory force running from −16% to +46% and from −28% to +78% of Newtonian across the central 90% of the
>    mass; the archive's own `support_fraction_max` is 0.83 to 1.30 and this report never mentioned it.
> 5. **The wide annuli are not equilibria of the declared distribution function.** `circular_energy()`
>    inverts `L_c(r)` by sorting and interpolating, and in these potentials that map is not monotonic
>    (kappa² < 0 over dozens of grid radii; three circular radii for some L), so E_circ(L) zig-zags between
>    branches, up to 5.5 dE from the true minimum, for 7 to 9% of the wide annuli's mass, and is clamped at
>    the low-L end. They are a narrow core plus a Kepler-like tail — rms widths 0.237 and 0.278 against the
>    0.135 and 0.174 FWHM the table shows. The narrow annuli are affected only beyond 4.5 widths.
> 6. **H5 shows a discrete fixed point, not an accurate equilibrium** (owner). `consistency_residual()`
>    reuses the solver's own quadrature. At fixed alpha the discretization error is about 7×10⁻⁴, not
>    10⁻¹¹. And a Gaussian factor in energy is not by itself a finite-mass guarantee: the implementation
>    relied on five-width windows and a bounded radial domain that were never declared as the
>    distribution's support, and with them the most extended populated orbit reaches r = 2.9 to 3.0, past
>    the simulator's limit at 2.25.
> 7. **"An equilibrium of the simulator … stationary to a few parts in 10⁴" is withdrawn as evidence.**
>    Window averages of r, its spread and rms v_r in a static axisymmetric field are close to orbit
>    invariants for *any* start. A sample with every body kicked outward by a full sigma_r shows drifts of
>    0.06%, 0.08% and 0.01% — just as "stationary" — and the mis-weighted sampler of item 1 passed unnoticed.
> 8. **"Stage 5's spectra were complete" and "exactly one unstable mode" are false below the region's
>    edge.** The certified region starts at Re s = 0.001, which excludes e-folding times beyond 159 T0, and
>    that edge was a constant in the code, not a declared limit. Below it the shipped two-stage determinant
>    has unstable roots at every m from 1 to 8 — m = 1: 5.583×10⁻⁴ + 1.04878i (e-folding 285 T0, a *static
>    lopsided* mode, pattern speed 2×10⁻⁶ Omega, on a branch continuation cannot reach) and 2.173×10⁻⁴ +
>    0.00649i (732 T0); m = 2: 1.038×10⁻⁴ + 0.00627i (1533 T0); m = 3: 1.868×10⁻⁴ + 2.47038i (852 T0) — each
>    with |det| below 10⁻¹³ and winding number 1.0000 on a circle wholly in Re s > 0. Beyn's own moments
>    located them and the rectangle filter threw them away. H3 passes only because both lists were
>    filtered by the same edge; it fails at m = 1 as soon as the edge moves to 5×10⁻⁴.
> 9. **Stage 5's continuation returns a duplicated root.** For the two-stage ring two of the four branches
>    collapse onto the same root at every m from 1 to 15, and the partner lost is the unstable one. The
>    stage 5 archive's m = 1 row reads "growth −1.4×10⁻⁶, e-folding infinite".
> 10. **H2's 1.8×10⁻¹⁶ measures Newton's polish, not the contour solver.** All three solves were polished
>    on the same determinant before being compared. Unpolished, the m = 2 root this stage rests on moves
>    5.8×10⁻⁶ under doubled quadrature and 4.5×10⁻⁶ under the shifted contour, **both above the declared
>    10⁻⁶**; the polish moved it by 1.2×10⁻⁵. The root is in no danger of being missed, only of being
>    mislocated without the polish, but the gate did not test what the protocol declared. H1's 10⁻¹⁷ is the
>    polish too (raw: 5×10⁻⁹), and its transcendental case kept 5 of the true 6 eigenvalues with
>    multiplicity while checking only the distinct count.
> 11. **Protocol deviations I did not name.** H3's tolerance is declared as 10⁻⁶ and coded as 10⁻⁵ (the
>    outcome is unchanged at 2.6×10⁻¹⁵, and 8 of the 9 two-stage rows compare an empty list with an empty
>    list). H6 was declared on "mean radius, **support** and radial dispersion" and evaluated on mean
>    radius, *radial spread* and rms v_r — and the omitted observable is the one that moved: the support on
>    the ring fell 22% in the warm run and rose 21% in the cold one over 30 T0. The protocol declared
>    "annuli of three widths"; two were run.
> 12. **H4's amplitude-scaling check could not fail**: it read the mode amplitude one step after release,
>    which is the imposed displacement, and the no-history seed scales by the same 10.0000. (A meaningful
>    version passes on the archived data: the fitted amplitudes at the start of the measurement window
>    scale as 9.990.) The archive also holds a two-term Prony estimate that this report never mentioned and
>    that misses by 10% at the larger amplitude; stage 5's seeded gate was judged on Prony and this one on
>    the line fit, and that change of estimator should have been stated. H4's rate and frequency agreement,
>    0.014% and 0.002%, stand, as does the no-field-history control.
> 13. Smaller: `run_ring_field` returns no status, so a terminated run would be read as complete (none of
>    the archived runs terminated); `prime()` zeroes the field inside r = 0.3; the FWHM column is quantised
>    to whole grid cells and low by 3 to 12% (interpolated: 0.099, 0.139, 0.141, 0.184); and the one-stage
>    ring has one located root at m = 1, not two.
>
> **Stage 7 (protocol-rut7.md) replaces the annulus construction, the sampler, H6 and the eigensolver
> gates, each with a negative control — a deliberately broken input the check must reject.**

What follows is the report as first written. Read it with the corrections above; claims they withdraw are
marked where they are made.

Part 1 is the two repairs the owner's review of ceb0c86 named as preconditions, plus the construction of
the state part 2 will test. **No configuration is claimed stable here.**

## R1. The spectrum is now located, not counted — ~~and stage 5's spectra turn out to have been complete~~ (withdrawn: corrections 8 to 10)

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
| one-stage | every m from 1 to 8, two roots per mode for m ≥ 2 and one at m = 1 | 4.451 T0 |
| **two-stage** | **m = 2 only, one root** | **8.579 T0** |

So the two-stage cold ring has exactly one unstable mode **with an e-folding time shorter than 159 T0**,
and continuation had found it. ~~Stage 5's spectra were complete; the only thing broken was the counter
that was supposed to certify them.~~ **Withdrawn (corrections 8 and 9):** slower unstable roots exist at
every m from 1 to 8, continuation misses most of them, and its root lists contain a duplicate.
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

## R3. A supported population, self-consistent to the stated residual on the declared discretization

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

~~**And it is an equilibrium of the simulator, not just of the equations.** Sampled to 64 bodies and evolved
with the field held fixed, the narrow warm annulus changes its mean radius by 0.02%, its radial spread by
0.06% and its dispersion by 0.08% over 20 periods.~~ **Withdrawn as evidence (correction 7):** those
window averages are near-invariants of any start in a static axisymmetric field, and a deliberately
kicked sample passes the same test.

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

## ~~The first evidence on temperature~~ — withdrawn (corrections 1 to 3)

The same two runs for the narrow **cold** annulus (σ_r = 0.094 against the warm annulus's 0.167):

| narrow cold, 20 T0 | radial spread | radial dispersion | angular momentum |
|---|---|---|---|
| field written | 0.03906 → 0.04448 (**+13.9%**) | 0.09071 → 0.10341 (**+14.0%**) | 0.99769 → 0.99570 |
| field frozen | 0.03905 → 0.03908 | 0.09077 → 0.09078 | unchanged |

~~Under the written field the colder annulus heats by 14% in twenty periods while the warmer one changes by
0.2%, and both are stationary when the field is frozen. That is the direction one expects if velocity
dispersion suppresses the collective instability, and it is the first thing in this programme that bears
on the owner's question.~~

**Withdrawn.** "Radial dispersion" in these tables is rms v_r, which a two-lobed stream raises without any
random motion changing. Both runs carry a growing, inertially stationary m = 2 field mode — e-folding 8.2
to 9.1 T0 in the warm run, 5.3 to 5.7 T0 in the cold — that azimuthal averages cannot see; the cold run's
residual dispersion changes by −0.6%. The warm annulus is not quiet. It is unstable on the cold ring's own
timescale, and its mode had simply not yet grown large enough in twenty periods to move an average.

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

**A note on the exit status, as the owner has now ruled.** `rut6.py` exits on whether it reproduces its
archive. The owner approved that for a reproduction job and corrected the reasoning: I presented a false
either/or, withhold a failed result or reword its gate. There is a third option, and it is the right one —
archive the outcome while reporting **three separate statuses**: *reproduction* (does the run reproduce
its archive), *numerical verification* (do identities, samplers, integrators and solvers pass their
correctness checks — a failure blocks every dependent physical conclusion), and *scientific outcome*
(does the model do what is hoped — a failure is a result). An incorrect sampler must not become acceptable
because its wrong output reproduces. Stage 7's driver implements all three; this one is left as it ran.

## Reproduce

```bash
python research_work/results/path-memory/rut6.py    # about ten minutes; also the suite job
```
