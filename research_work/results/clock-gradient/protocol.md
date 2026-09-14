# Clock-gradient check CG-0: is the rolling frame well defined?

Declared before execution, 14 September 2026. Baseline: `main` at 6ce4c7a.

**Why this check comes first.** The project owner's review asked for it before any clock-field gravity model (CG-1) is built.

## Question

The unified proposal in the [model contract](../../../research_plan/model-contract.md) (C2–C3) uses the clock field's gradient as a rolling frame:
- u_μ = ∂_μχ/√X, with X = (∂_tχ/c)² − |∇χ|².

This is only defined where X > 0, and it breaks down as X → 0. Does X stay positive, and bounded away from zero, across the transition from low-acceleration to Newtonian gravity, in the systems the proposal must describe?

## Declared setup

- **Rolling.** χ = ln n, so ∂_tχ = ṅ/n, the redshift rate. This check uses c·(ṅ/n) = c²α = 7.2496×10⁻¹⁰ m/s², the value behind RPG-1's link a* = ξc²α with ξ = 0.118.
- **Static response.** Around matter, χ acquires a static perturbation δχ(x). Matter feels a_χ = c²|∇δχ| (the normalization in which C2 ties a* to c|∂_tχ|). Then
  - X/X₀ = 1 − (a_χ / c²α)²,
  - with X₀ = (∂_tχ/c)² the cosmological value.
- **Three formulations of the static response:**
  - **F1, the contract as written.** The clock field's static equation is RPG-1's AQUAL equation for the total potential, so a_χ = g, the total acceleration.
  - **F2, excess only.** The clock field carries only the departure from Newtonian gravity: a_χ = g_N[ν(g_N/a*) − 1] with the simple ν. This is a bi-potential (QUMOND-like) structure whose source is built from the Newtonian field.
  - **F3, a limiting-gradient kinetic function.** It keeps X > 0 by making the field's response stiffen as X → 0. The static force then saturates at c²α near every mass.
- **Data.** All are exposed and archived; none is new.
  - **SPARC.** 149 galaxies, with g = V_obs²/R and g_N from the archived Υ_disk = 0.5 and Υ_bulge = 0.7.
  - **The Milky Way.** Baryon model I, spherically averaged, from 0.5 to 25 kpc.
  - **The six SLACS lenses.** RPG-1's g_N/a* at the Einstein radii, 9–14.
  - **The Solar System.** Newtonian g from 0.39 to 40 AU. The bound on an anomalous constant radial acceleration at Saturn is below 10⁻¹⁴ m/s², from Cassini ranging, as reviewed in [Turyshev & Toth 2010](https://arxiv.org/abs/1001.3686).

## Outputs

- X/X₀ along each system under F1–F3.
- The fraction of SPARC radii, and the Milky Way radius, where X ≤ 0.
- The anomalous Solar System acceleration under F2 and F3.
- A statement of whether an independent vector or second field is needed.

## Pass rule (declared)

A formulation passes if both hold:
- X/X₀ ≥ 0.1 at every data point;
- its anomalous acceleration at Saturn is below 10⁻¹⁴ m/s².

## Not claimed

- A gravity model.
- That structures not tested here fail. Examples are an independent timelike vector (as in AeST) and a kinetic function with a second scale.
- Any fit. Nothing is adjusted.
