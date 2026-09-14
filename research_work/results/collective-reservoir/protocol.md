# Driven collective reservoir (CR-1): a counted seed under a derived photon source

Declared before execution, 13 September 2026. Baseline: `main` at 660a4f3. This is the third of three branches opened by the project owner.

CR-1 asks whether a coherent reservoir, started from a counted seed and driven by a source derived from the project's photon conversion law, grows, disperses, overheats or collapses. The Gross–Pitaevskii–Poisson equations are established mathematics. Their use for a photon-fed reservoir, and every channel below, is a hypothesis.

Before this declaration, only development runs of the solver on the archived analytic Schrödinger–Poisson ground state were made. No seed or source run has been computed.

## What the repository already fixes

- **Occupation.** Any coupling to the density |ψ|² conserves the number of quanta ([bound-cloud-exchange](../bound-cloud-exchange/report.md)). A source must therefore break that symmetry and name its channel.
- **Stimulated balance.** Creation and removal are both Bose-enhanced, so the steady occupation of incoherent bound modes is finite ([pair-production-balance](../pair-production-balance/report.md)).
- **Supply.** Converting all present luminosity for 10 Gyr falls short of the deposited mass-energy by a median factor of about 5,400 ([recovered supply report](../../../companion_causal_test/report.md)). The local conversion law falls short by about 10⁹.
- **Shape.** A single stationary cloud has a fixed shape and an eventually declining exterior speed ([companion-self-binding](../companion-self-binding/report.md), [oscillating-field-profile](../oscillating-field-profile/)).

## Postulates

- **Q1. The reservoir.** One coherent, nonrelativistic scalar condensate: iħψ_t = [−ħ²∇²/2m + m(Φ_self + Φ_b)]ψ, with ∇²Φ_self = 4πG m|ψ|².
  - There is no self-interaction in this version (g = 0).
  - The geometry is spherical with ℓ = 0.
  - Φ_b is the spherically averaged Milky Way variant I baryonic potential.
  - Two constituent masses are used: m = 1.34×10⁻²⁴ eV, the archived self-binding value (a galaxy-sized soliton), and m = 1×10⁻²² eV (the compact regime).
- **Q2. The counted seed.** N0 = M_seed/m quanta with M_seed = 1.0×10⁹ Msun (1% of the variant I baryons), in the stationary ground state of the combined potential.
- **Q3. The derived source.** The conversion law sets the mass-rate density q(r) = αc u_γ(r)/c², with α = 2.488993×10⁻⁴ Mpc⁻¹.
  - u_γ is the Milky Way's starlight: optically thin, with emissivity proportional to the spherically averaged stellar density and L_bol = L[3.6] = 1.79×10¹¹ Lsun (an upper bound). The microwave background's 4.17×10⁻¹⁴ J/m³ is added.
  - Each created quantum costs mc²; binding (|μ|/mc² < 10⁻⁶) is neglected.
  - Stimulated creation needs the occupied mode, so the source acts only inside r₉₉, the radius enclosing 99% of the current condensate mass. This also keeps the source independent of the numerical box.
- **Q4. Channels.**
  - **A. Products take the photon's momentum.** With no third body, created quanta carry at least 1/c of momentum per unit of created rest energy. They are relativistic and unbound, so bound growth is zero. This is derived analytically, and the seed is evolved without a source.
  - **B1. Collective recoil, stimulated into the occupied mode.** The whole reservoir takes the recoil, and all mass converted inside r₉₉ enters the ground mode without changing its shape. This is an optimistic bound. The recoil-free fraction for optical photons on a kiloparsec condensate is essentially zero, and no mechanism for it is proposed.
  - **B2. Collective recoil at the local phase.** Mass converted inside r₉₉ is added where it is converted, with the local condensate phase. Its shape follows the photons, not the mode.
  - **C. Required-rate diagnostic (not derived).** B1 and B2 with q scaled by the factor that, at the seed's initial collection rate, would give the archived exact-third reference inventory of the Milky Way fiducial in 10 Gyr. This assumes the supply problem is solved and tests only the formation dynamics.

## Outcomes, per mass and channel

Each run lasts 10 Gyr, or until the solution leaves the resolved or nonrelativistic domain.

- **Growth.** N(t)/N0 and M(t).
- **Dispersal.** Mass removed by the absorbing boundary.
- **Overheating.** The excess energy E_exc = E − E_0(N) above the ground state at the current N.
- **Collapse.**
  - r_half(t) and the central density;
  - M against the Kaup limit 0.633ħc/(Gm), and r_half against the Schwarzschild radius;
  - runs stop when r_half < 10 grid cells, and the self-gravitating ground-state sequence r_half = 3.925ħ²/(GMm²) then gives the endpoint analytically.
- **Energy ledger.** Injected mass-energy against the change in field mass plus the escaped mass.
- **Rotation.** The reservoir's contribution to the Milky Way rotation curve, added to the archived baryon speeds at the 38 Eilers radii. This is reported, not scored.

## Validation (tolerances declared now)

- **V1. Ground state.** The Schrödinger–Poisson ground state in units ħ = m = G = M = 1 must match the archived E = −0.0542564, μ = −0.1627692 and r_half = 3.92510 within 1e-4 relative, with a virial residual below 1e-4.
- **V2. Stationarity.** Evolved in real time for 100 dynamical times, the ground state must conserve energy to 1e-8 relative, and its density must change by less than 1e-3.
- **V3. Absorbing boundary.** An outgoing free packet must be at least 99.9% absorbed, with reflection below 1e-3 of its norm.
- **V4. Source.**
  - Growth into the mode at a fixed rate must reproduce N(t) = N0 + Ṅt to 1e-8. At a slow rate, the ground-state shape must be preserved: excess energy below 1e-3 of |E|.
  - The ledger must close to 1e-6.

## Assessment rule

A channel supports a collective reservoir only if all of the following hold:
- it grows the seed to the required inventory within 10 Gyr;
- it does not overheat (E_exc/|E| < 0.1);
- it does not disperse (escaped mass below 10%);
- it stays in the nonrelativistic domain.

Channel C tests only dynamics and cannot support formation on its own.
