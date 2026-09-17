# RUT-1 stage 5: which equilibria the written-track response supports, and what the m = 2 mode actually is

Declared before any stage 5 computation. **This stage has its own protocol file**, and stage 4's archive,
code and protocol stay frozen: `protocol-rut1.md`, `formation.py`, `longrun.py`, `rut3.py`, `rut4.py` and
`rut4-series/` are not edited by this stage, because `rut4_checks.py` pins their hashes to what ran. New
diagnostics and new source variants get new files.

## Why this stage exists

Stage 4 measured a large-scale m = 2 instability of a cold, mature, supported ring and reported it as a
property of the two-stage response. **That attribution was not controlled.** The owner's review supplied
the missing control: a *zero-delay* field with the same kernel and the same static gain,

    C_instantaneous(x, t) = tau_keep * sum_i q_i exp[-|x - X_i(t)|^2 / (2 w^2)]

which is not "memory off" — it keeps the extra attraction and removes only the lag. Their exploratory
linearization of the 32-writer, w/R = 0.2 configuration finds that this control is *also* m = 2 unstable,
with an amplitude e-folding time near 0.191 initial reference periods. If that reproduces here, then the
correct statement is not "the two-stage response is unstable" but "collective attraction destabilises this
cold ring, and the two-stage lag slows it by a large factor".

So stage 5 does not try to eliminate the mode. It asks **what the mode is, and which equilibria this
response can support.**

## The declared ladder

Every rung shares the kernel, the writing normalization and the geometry; they differ only in how a
disturbance becomes force-producing. In Laplace form the linear field response to a source perturbation is
`C_hat = H(s) S_hat`:

| rung | H(s) | note |
|---|---|---|
| no attraction | 0 | Kepler ring; every mode must be neutral |
| instantaneous | `tau_keep` | the owner's matched control: same static gain, no lag |
| one-stage | `tau_keep/(1 + s tau_keep)` | `dC/dt = S - C/tau_keep` |
| two-stage | `tau_keep/[(1 + s tau_keep)(1 + s tau_form)]` | the stage 3-4 candidate |

A fifth rung is declared as an identification test: **pairwise Newtonian attraction** between the same
bodies, at the strength that matches the instantaneous rung's equilibrium inward support. If the cold ring
behaves the same way there, the mode is the classical instability of a discrete ring of co-orbiting
attracting bodies and has nothing to do with memory.

## The linear method, declared in full

Base state: `N` equal writers on a circle of radius `R0 = 1`, equally spaced, rotating rigidly at `Omega`.
The memory force on body `i` is the exact history integral

    a_mem,i(t) = sum_j q_j integral_0^inf h(u) grad G(X_i(t) - X_j(t-u)) du,     G(d) = exp(-|d|^2/2w^2)

with `h` the inverse Laplace transform of `H`: `tau_keep*delta(u)` instantaneous, `exp(-u/tau_keep)`
one-stage, `tau_keep/(tau_keep - tau_form) [exp(-u/tau_keep) - exp(-u/tau_form)]` two-stage. `Omega` is
solved self-consistently with the radial balance `Omega^2 R0 = GM/R0^2 - a_mem,radial`.

The base state is only *quasi*-stationary for the memory rungs: the body's own wake exerts a tangential
force, so the ring loses angular momentum slowly. **That force is measured and reported for every
configuration**, and the frozen-`Omega` assumption is declared rather than assumed silently.

Linearising about that state in the rotating frame, with displacements `u_j` and the discrete Bloch
decomposition `u_j = u_hat exp(2*pi*i*m*j/N) exp(s t)`, the mode condition is

    det[ s^2 I - 2 Omega s J - Omega^2 I - D - A + B_m(s) ] = 0,   J = [[0,1],[-1,0]]

with `D = diag(2GM/R0^3, -GM/R0^3)` the central term, `A = grad grad C^0` at the body (from the same
history integral) and

    B_m(s) = q sum_j integral_0^inf h(u) exp(-s u) exp(2*pi*i*m*j/N) T_j(u) Rot(-Omega u) du,
    T_j(u) = grad grad G(P_0 - Rot(-Omega u) P_j).

For the instantaneous rung `B_m` is independent of `s` and the condition is an exact quartic; for the
memory rungs it is a nonlinear eigenvalue problem, solved by Newton iteration from the instantaneous root
and reported with its residual. Growth rate is `Re(s)`; pattern speed is `Omega - Im(s)/m` in the inertial
frame.

**Matching.** The owner's control matches the static gain. That leaves the rungs with different *net*
support, because a memory field includes the body's own wake and an instantaneous one cannot. Both are
therefore reported: matched-gain (primary, as declared by the owner) and matched-support, where `q` is
rescaled so the equilibrium inward acceleration equals the two-stage rung's.

## Declared computations

1. **The ladder**, all five rungs, at the two stage 4 configurations (16 writers w/R = 0.1, 32 writers
   w/R = 0.2), every mode `m = 0 .. N-1`, both matchings.
2. **A writing-strength scan** at fixed geometry: label fraction 0, 0.025, 0.05, 0.10 (the stage 4 value),
   0.20, for the instantaneous and two-stage rungs. This asks whether a supported cold ring exists at all,
   and if so below what support.
3. **A discreteness scan** at fixed total writing rate: N = 8, 16, 32, 64, 128. If the growth rate falls
   steeply with N the mode is a property of the discrete writers; if it converges it is a property of the
   supported ring itself.
4. **Nonlinear verification** of the linear prediction by seeding an eigenmode at two amplitudes.
5. **Archive diagnostics on the frozen stage 4 series** (read-only): instantaneous source-pattern speed and
   growth rate at the saved checkpoints from the analytic `S` and `dS/dt`, and a separation of coherent
   streaming from residual dispersion by a low-order azimuthal fit.

## Gates, with thresholds fixed now

| gate | requirement |
|---|---|
| G0 base state against stage 2 | one writer, one-stage: radial and tangential base force equal `rut1.self_force` to 1e-6 relative |
| G1 response matrices | `A` and `B_m` reproduce centred finite differences of the full history force to 1e-6 relative |
| G2 integrator | the stage 5 integrator reproduces `formation.run` (two-stage, 3 T0, same step and grid) to 1e-6 relative in position |
| G3 eigenmode, instantaneous | measured growth rate and pattern speed within 2% of the linear prediction over 3 T0, at seed amplitudes differing by 10x |
| G4 eigenmode, two-stage | within 10% of the linear prediction over 10 T0 |
| G5 no-attraction control | `max Re(s)*T0 < 1e-8` over all modes |
| G6 quadrature | doubling the quadrature nodes and the history range moves `Re(s)` by < 1e-4 relative; Newton residual `\|det\| < 1e-10` |
| G7 pattern-speed recovery | a synthetic rigidly rotating pattern returns its known speed to 1e-6 relative |
| G8 streaming decomposition | total variance equals coherent plus residual to 1e-12, and a synthetic flow with known coherent and random parts is recovered to 1% |

A gate that fails is reported as a failure; no threshold is relaxed after seeing a number.

## What this stage cannot establish

It is a **linear** analysis about one family of equilibria: a cold, equally spaced, rigidly rotating ring of
equal writers. It says nothing about the saturated nonlinear state beyond the seeded verification, and a
linear growth rate is not a lifetime. It does not test finite-width annuli, velocity distributions,
unequal masses, three dimensions, or overlapping populations — those are stage 6, declared here and not
run. It does not change the source law, identify the dissipation reservoir, or add spatial causality. The
archive diagnostics are limited by what stage 4 saved: five checkpoints of particle state, from which an
*instantaneous* source-pattern speed can be computed but not the continuous phase history of the
force-producing field.

## Declared for stage 6, not run here

Finite-width annuli and specified velocity distributions built in equilibrium with the declared potential;
the mode map around those states; the reciprocal accounting carried through growth and saturation;
complex mode coefficients of `S`, `E` and `C` saved per checkpoint; and full restart states including the
fields, so a run can be continued rather than restarted from particles alone.
