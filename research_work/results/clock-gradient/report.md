# CG-0 report: the clock gradient is not a usable rolling frame on its own

**Result.** Using one scalar's gradient as both the rolling frame and the galaxy force does not survive the transition to Newtonian gravity.
- **F1 (the model contract as written).** The gradient turns spacelike (X < 0) at every SLACS Einstein radius, inside 2.0 kpc in the Milky Way, at 5.6% of SPARC radii (23 of 149 galaxies), and throughout the Solar System.
- **F3 (a limiting gradient).** X stays positive but falls to zero exactly there. The force then saturates at c²α, 7×10⁴ times the Cassini bound.
- **F2 (excess only).** The frame stays well defined everywhere (X/X₀ ≥ 0.986). But the clock field's source is built from the Newtonian field, so F2 is a two-field (bi-potential) structure, and with the simple ν its Solar System excess exceeds the Cassini bound.

**No formulation passes the declared rule.** An independent timelike vector (as in AeST) or a second field is needed. The contract's one-field proposal (C2–C3) is withdrawn in that form.

Protocol: [protocol.md](protocol.md), declared in bc4a64b before execution. Script: [timelike.py](timelike.py). Results: [timelike-results.json](timelike-results.json).

## Why the transition meets X = 0

- **The link.** C2 ties the static response to the rolling: with χ = ln n and a_χ = c²|∇δχ|,
  - X/X₀ = 1 − (a_χ/c²α)²,
  - c²α = 7.25×10⁻¹⁰ m/s².
- **Where the two scales sit.** The archived transition scale is a* = ξc²α with ξ = 0.118. So X reaches zero at a_χ = c²α = 8.47 a*, only a factor 1/ξ above the transition. At a_χ = a*, X/X₀ is still 0.986.
- **Consequence for F1.** If the field carries the total acceleration, every system with g > 8.5 a* crosses X = 0.
- **Why no single-scalar version avoids it.** A single scalar with a monotone flux law has only two options as the source grows:
  - its gradient grows without bound (F1);
  - it approaches a limiting value (F3).

  It cannot fall back to zero. Carrying only a vanishing excess (F2) needs a source built from a second field.

## Results

| System | F1: X/X₀ | F2: X/X₀ | F3 |
|---|---|---|---|
| SPARC, 3,152 radii with R, V > 0 | ≤ 0 at 5.6% of radii in 23 galaxies; minimum −701 | ≥ 0.986 | → 0 where F1 < 0 |
| Milky Way, model I, 0.5–25 kpc | X = 0 at R = 2.02 kpc; −29.7 at 0.5 kpc | ≥ 0.987 | → 0 inside 2 kpc |
| SLACS Einstein radii, 6 lenses × 2 IMFs, g_N/a* = 9.3–13.6 | −0.45 to −1.94: negative at every lens | 0.988 | → 0 |
| Solar System, 0.39–40 AU | −3×10¹⁵ (Mercury) to −2.6×10⁷ (40 AU) | 0.986 | → 0 |

Anomalous radial acceleration at Saturn, isolated Sun:

| | F1 | F2 | F3 | Bound (Cassini, reviewed in [Turyshev & Toth 2010](https://arxiv.org/abs/1001.3686)) |
|---|---|---|---|---|
| m/s² | 8.6×10⁻¹¹ | 8.6×10⁻¹¹ | 7.25×10⁻¹⁰ | < 10⁻¹⁴ |

- **F1 and F2 share the simple ν.** Its excess tends to the constant a* in strong fields, which is a known problem for that interpolating function independent of this check. A ν whose excess falls fast enough would pass the Solar System bound but was not tested.
- **The Galaxy's external field is not included.**

## What this means for the proposal

- **Lensing (C3).** It relies on the rolling frame, which F1 cannot define at the lenses.
- **Keeping one field responsible for the redshift.** The galaxy-scale response needs one of:
  - an independent unit timelike vector that defines the frame, with the scalar supplying the force (the AeST structure);
  - a bi-potential construction (F2) in which the clock field is sourced through the Newtonian potential. This is a second field in all but name.
- **Before CG-1.** Any CG-1 protocol must name its second field or vector. It must test Milky Way radial and vertical forces first, as the review asks, and state a Solar System–safe interpolating function.

## Not claimed

- That AeST-type or bi-potential completions work.
- Any fit. Nothing was adjusted.
- The data are exposed archived samples; the SLACS accelerations come from RPG-1's lens models.
