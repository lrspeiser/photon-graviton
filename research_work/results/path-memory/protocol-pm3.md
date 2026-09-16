# PM-3: what a shared-field memory can and cannot do

Declared before execution, 16 September 2026. Baseline: `main` at d621ebf.

**The question.** The owner's review supplies a toy in which bounded memory makes disturbed orbits settle, stably, with a closed energy ledger. It removes a false obstacle: our earlier feedback equation's instability was a property of that constitutive law, not a theorem against delayed response. The question PM-3 asks is the next one — **does a memory belonging to one shared field produce rotation-dependent gravity in a galaxy, and if so, what must the field be told about the source?**

**Why it is run before the field solver.** The owner's own caveat is the crux: "Remembering a field is not automatically the same as detecting orbital motion. A perfectly smooth, stationary axisymmetric density does not change merely because its constituent stars rotate." PM-1's candidate A already measured the density version of that — a rigidly rotating axisymmetric source is indistinguishable from a static one to 3.7×10⁻¹¹. Before building a propagating memory field, PM-3 settles by proof what such a field could contribute in a steady state, because the answer determines what the field must couple to and therefore what is worth implementing.

**What this is not.** PM-3 solves no galaxy, fits nothing, and scores no rotation curve. It reproduces a declared toy, proves a scope statement about it, and states what the next construction has to do. Nothing here bears on whether extra attraction exists, only on which mechanisms could supply it in a steady state.

## Stage 1: reproduce the owner's toy as a control

Exactly the declared system, per unit test-body mass, with `z = ln(r/r_ref) − χ`:

    g_inward = GM/r² + [K + B·tanh(z)]/r
    r̈ = j²/r³ − g_inward
    τ·dχ/dt = tanh(z)
    θ̇ = j/r²

with K > 0, 0 < B < K, τ > 0. The response is bounded into (K−B, K+B)/r, so the total force stays attractive and the memory is a restoring correction rather than a ratchet.

**Verifications, each gated:**
1. **Equilibrium identity.** j² = GM·r_c + K·r_c², to 10⁻¹².
2. **Energy ledger.** With E = ṙ²/2 + j²/(2r²) − GM/r + K·ln(r/r_ref) + B·ln(cosh z) and dQ/dt = +(B/τ)tanh²z, the quantity E + Q − W is conserved, where W is the work done by any preparation driver. Gated at 10⁻¹⁰ over 40 reference periods. dE/dt = −(B/τ)tanh²(z) ≤ 0 is verified symbolically by the same cancellation that produces it.
3. **Settling.** Five disturbed starts reach |r/r_c − 1| < 10⁻⁶ during the fortieth reference period.
4. **The control that must fail.** The same law with B = 0 must *not* circularize; its final-period radius range is reported.
5. **Stability, analytically.** The linearization gives τλ³ + λ² + τ(κ²+b)λ + κ² = 0 with κ² = GM/r_c³ + 2K/r_c² and b = B/r_c². Routh–Hurwitz reduces to τ(κ²+b) > τκ², i.e. **b > 0** — so every positive finite τ is stable, in contrast with candidate C, where the extra term sat in the constant coefficient instead and forced q₀ < 0. The decay rate is reported across τ, since it vanishes at both τ → 0 and τ → ∞ and stability does not imply prompt settling.
6. **Two preparation histories.** Two runs that begin identically and are released at the **same** r, ṙ and θ, differing only in when a prescribed excursion occurred, must give different inward accelerations, with the driver's work carried in the ledger. The excursion amplitude is prescribed, not fitted, so the solver cannot satisfy the match by removing the excursion.

## Stage 2: the steady-state scope of any relaxational memory

**The theorem.** Let a memory state obey τ·dQ/dt = F[fields] − Q for any functional F of the instantaneous fields. In a steady state the fields are time-independent, so Q relaxes to F[fields] and is itself time-independent. The memory then contributes a definite term fixed by the instantaneous fields — **nothing that a static constitutive law of those same fields could not also produce.** Memory buys transients; in a steady state it buys nothing.

**Corollary 1, inside the owner's own toy.** On the settled orbit z = 0, so the extra force is exactly K/r whatever B and τ are: the memory contributes identically zero to the circular orbit it produced. Verified at B = 0.2, 0.4, 0.8. The memory does the **settling**; the flat-curve term K/r is assumed, as the owner states. This is a scope statement about the toy, not a defect in it.

**Corollary 2, for the proposed field extension.** For the owner's candidate — τ·dQ/dt = p − Q with p = ∇ψ, entering the field equation through η(p − Q) — a smooth axisymmetric galaxy in steady rotation has a stationary density, hence stationary ψ and p, hence Q → p with |p − Q| decaying as e^{−t/τ}. The memory flux vanishes and the equation returns **exactly** to Completion I, at every τ and every η. The extension is therefore blind to steady rotation by construction. This is the owner's caveat made exact.

**What follows, and it is constructive.** Two things, and PM-3 states both rather than choosing between them:
- **Stage C is already the right test of the steady-state content.** If smooth steady memory is equivalent to a static law, then PM-2A's static completions span that whole space, and the disk comparison is testing the correct object. The memory extension is not a shortcut around it.
- **Steady rotation-dependent gravity must come from somewhere else.** Either a coupling to a quantity that is stationary but nonzero for a rotating disk and zero for a static one — the mass current ρv, or a velocity stress — or from genuine non-steadiness, which in a galaxy means the discreteness of stars rather than a smooth fluid. That second route is the one the walking-droplet analogy actually describes, since a droplet remembers its own past positions.

## Stage 3: what the field must be told, declared and not yet run

The minimal well-posed successor, declared here so its result cannot be chosen after the fact: a source coupling that sees the mass current. A rigidly rotating axisymmetric disk has stationary ρ *and* stationary J = ρv ≠ 0, while a static disk has J = 0, so a term in J distinguishes them without any memory at all — which is itself the point, because it shows the distinguishing ingredient is a coupling, not a history. Its first obligations are the ones that killed earlier candidates: it must not depend on how the same matter is divided into numerical particles, it must conserve momentum with a central-force limit, and it must reduce to the verified static equation when J → 0.

PM-3 does not run stage 3. It declares what would make stage 3 informative, so that a later run cannot present a coupling chosen for its fit as a prediction.

## What would make this informative either way

If stage 1 reproduces and stage 2's theorem holds, the programme gains a proven scope boundary: memory is a transient mechanism and the steady galaxy problem belongs to the static completions plus, possibly, a current coupling. If stage 1 failed to reproduce, the toy's claim would be in doubt and that would be reported as such. Neither outcome is a statement about whether extra attraction exists.

A settling toy is not a galaxy theory, a numerical relaxation time is not a physical memory time, and this protocol establishes neither path memory nor a shared medium.

## Files

`pm3.py` runs stages 1 and 2 and writes `pm3-results.json`; it is registered as its own suite job. No file under any other experiment's directory is modified.
