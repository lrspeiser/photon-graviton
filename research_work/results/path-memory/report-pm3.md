# PM-3: the memory does the settling, not the support — and a field memory is blind to steady rotation

[protocol-pm3.md](protocol-pm3.md), declared in 9f57c99 before this run. PM-3 reproduces the owner's toy as a control and then asks what a memory belonging to one shared field could contribute in a galaxy. It solves no galaxy, fits nothing and scores no rotation curve.

**Both stages pass.**

## Stage 1: the toy reproduces, and its stability is not luck

| Gate | Requirement | Measured |
|---|---|---|
| Equilibrium identity | j² = GM·r_c + K·r_c² to 10⁻¹² | **0.0** exactly (r_c = 1) |
| Energy ledger E + Q − W | conserved to 10⁻¹⁰ over 40 reference periods | **7.6×10⁻¹³** |
| Settling, five disturbed starts | \|r/r_c − 1\| < 10⁻⁶ in the fortieth period | **2.2×10⁻⁸** |
| The B = 0 control | must **not** circularize | radius 0.801–1.300, spread 0.499 |
| Routh–Hurwitz across τ | stable at every positive finite τ | max Re λ < 0 at all nine τ |
| Two preparation histories | same released state, different acceleration | **0.999576 vs 1.001591** |

The energy statement is exact, not numerical: differentiating E = ṙ²/2 + j²/(2r²) − GM/r + K·ln(r/r_ref) + B·ln(cosh z) along the flow cancels every term against r̈ except **dE/dt = −(B/τ)tanh²(z) ≤ 0**. Carrying Q with dQ/dt = +(B/τ)tanh²z closes it. The memory can only remove disequilibrium energy; where Q goes physically remains a hypothesis, and tracking it only keeps the accounting honest.

**The stability is structural.** The linearization gives τλ³ + λ² + τ(κ²+b)λ + κ² = 0 with κ² = GM/r_c³ + 2K/r_c² and b = B/r_c², and Routh–Hurwitz reduces to τ(κ²+b) > τκ², i.e. **b > 0**. Candidate C's cubic put its extra term in the *constant* coefficient and so required q₀ < 0; here it sits in the *linear* coefficient. That single placement is the whole difference between "unstable for every positive memory time" and "stable for every positive memory time". Decay is fastest near τ ≈ 1 (max Re λ = −0.0679) and vanishes at both ends — −2.0×10⁻⁴ at τ = 10⁻³ and −9.1×10⁻⁵ at τ = 10³ — so stability does not imply prompt settling.

**Two histories, same present state.** Two runs begin identically, carry a *prescribed* excursion at t = 1 or t = 4, and are released at the same r, ṙ and θ to 1.4×10⁻¹³. Their memories differ (χ = +0.00106 against −0.00398) and so do their inward accelerations, by 0.20%. The driver's work is carried in the ledger, which closes to 2×10⁻¹³. Prescribing the excursion amplitude matters: with it free, the shooting collapses to the trivial undisturbed orbit and the "history dependence" is vacuous. The owner's own numbers were 0.971150 and 0.925548, a larger effect; their preparation driver was not specified, so this reproduces the structure and not that particular pair.

## Stage 2: what a memory can do in a steady state — nothing

**The theorem.** Let a memory obey τ·dQ/dt = F[fields] − Q for *any* functional F of the instantaneous fields. In a steady state the fields are time-independent, so Q relaxes to F[fields] and is itself time-independent. The memory then contributes a term fixed by the instantaneous fields — nothing a static constitutive law of those same fields could not also produce. **Memory buys transients. In a steady state it buys nothing.**

**Corollary 1, inside the toy itself.** On the settled orbit z = 0, so the extra force is exactly K/r whatever B and τ are. The identity is algebra and holds at 0.0; numerically the residual memory force falls to the integrator's floor as the orbit settles:

| B | memory force at 64 periods | at 128 periods | static part |
|---|---|---|---|
| 0.0 | — (never circularizes, r = 0.959 at 128P) | — | — |
| 0.2 | +2.7×10⁻⁸ | +6.0×10⁻¹⁴ | 0.800 |
| 0.4 | −2.2×10⁻¹⁴ | −4.0×10⁻¹⁴ | 0.800 |
| 0.8 | +5.6×10⁻¹⁴ | −1.4×10⁻¹⁴ | 0.800 |

**The memory produces the settling; the flat-curve term K/r is assumed and is what the settled orbit rides on.** The owner says this too — "the equilibrium K/r term is assumed in this toy" — and it is worth stating as a measurement rather than an aside, because it locates exactly what the toy does and does not supply.

**Corollary 2, for the proposed field extension.** For τ·dQ/dt = p − Q with p = ∇ψ, entering through η(p − Q): a smooth axisymmetric galaxy in steady rotation has a **stationary** density, hence stationary ψ and p, hence Q → p with |p − Q| decaying as e^{−t/τ}. The memory flux vanishes and the field equation returns **exactly to Completion I** — at every τ and every η. The extension is blind to steady rotation by construction.

This is the owner's own caveat made exact, and PM-1 already measured the density version of it: a rigidly rotating axisymmetric source is indistinguishable from a static one to 3.7×10⁻¹¹.

## What follows, and it is constructive

Two things, and PM-3 states both rather than choosing between them.

**PM-2A's stage C is already the right test of the steady-state content.** If smooth steady memory is equivalent to a static law, the static completions span that whole space. The memory extension is not a shortcut around the disk comparison; it is a different question about transients.

**Steady rotation-dependent gravity has to come from somewhere else.** Either a coupling to a quantity that is stationary but nonzero for a rotating disk and zero for a static one — the mass current ρv, or a velocity stress — or from genuine non-steadiness, which in a galaxy means the discreteness of stars rather than a smooth fluid. That second route is what the walking-droplet analogy actually describes: a droplet remembers *its own past positions*, and a single body's r(t) is genuinely time-dependent, which is precisely why the toy works and a smooth steady disk would not.

Note what the first route implies: if a stationary current coupling is what distinguishes a rotating disk, then the distinguishing ingredient is a **coupling, not a history**, and no memory is needed for it at all. That is declared as stage 3 and deliberately not run, with its obligations fixed in advance — independence of how the same matter is divided into numerical particles, momentum conservation with a central-force limit, and reduction to the verified static equation as the current vanishes.

## What this does and does not establish

**Shown.** The owner's toy reproduces: bounded memory settles disturbed orbits, its energy ledger closes exactly, its stability is a Routh–Hurwitz theorem rather than a lucky trajectory, and two histories ending in the same state respond differently. Candidate C's instability was a property of that constitutive law, not of delayed response in general — that obstacle is removed. And any relaxational memory is equivalent to a static law in a steady state, which bounds what the field extension could contribute before it is built.

**Not shown.** No galaxy, no rotation curve, no shared medium, no path memory. One test body around a fixed centre, with no spatial propagation. The toy does not derive K, and its memory contributes exactly zero to the orbit it settles into. A numerical relaxation time is never a physical memory time.

## Reproduce

```bash
python research_work/results/path-memory/pm3.py
```
