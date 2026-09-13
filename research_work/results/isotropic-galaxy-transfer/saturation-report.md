# Capacity-limited companion settling

Three receiving-capacity rules were tested while keeping the one-third capture law fixed. Incoming packets fill available bound states inward of their original location; blocked packets try farther-out cells within their permitted migration interval. Any remainder stays in an explicitly counted traveling reservoir. No missing packet energy is renormalized into the bound profile.

## Proposed rules and formula provenance

Use mass-equivalent capacity rho_cap=u_cap/c^2. The E=mc^2 conversion and shell bookkeeping are known mathematics. These companion capacity rules are proposed, not derived physics:

| Rule | Bound-state capacity density |
|---|---|
| Constant | rho_cap=C |
| Stellar density | rho_cap=zeta*rho_stars |
| Baryonic gravity | rho_cap=C*g_b/(1e-10 m/s^2) |

Stellar density and acceleration are estimated with a known spherical Hernquist form, M_star(<r)=M_star*r^2/(r+a_star)^2 and a_star=Re/1.8153. This is a declared approximate proxy with the original companion-fit stellar mass, not a measured three-dimensional stellar distribution or an NFW-derived driving field. No halo target mass was used to set capacity or supplied inventory. The mass inputs themselves were fitted to prior stellar data, so they are not wholly independent measurements.

For each source packet, the allowed interval is exp(-b)*r to r. The algorithm fills the innermost available cell, continues outward if full, and records unbound overflow. b is a phenomenological migration-rate-times-duration parameter; no universe age is imposed. This is an allocation rule, not a solved time-dependent drift equation or a gravitational force law. Source packets are processed inside-out; reversing order is reported as a sensitivity check.

The supplied inventory is the original deposit mass inside diagnostic 5 Re, redistributed according to each conditional stream profile. R is not a measured halo boundary. Cumulative bound mass is divided by SUPPLIED mass, so its endpoint can be below one. The target is still a normalized halo profile; this does not verify the absolute halo mass, original photon supply or projected lensing. Traveling overflow also gravitates, but its field and stresses are not solved here; these comparisons concern the bound component only.

## Shared fits and transfer

Two parameters per rule were fitted on a declared grid. Each omitted-galaxy result fits the other five and freezes the rule for the omitted target. All targets were previously inspected. Numbers below are RMS cumulative-profile errors in percentage points, not velocity/lens errors or significance.

| Geometry | Source | Capacity | No migration | Shared fit | Omitted-galaxy transfer |
|---|---|---|---:|---:|---:|
| companion_regular | distant | constant | 26.69 | 22.18 | 24.64 |
| companion_regular | distant | stellar_density | 26.69 | 22.18 | 28.08 |
| companion_regular | distant | baryonic_gravity | 26.69 | 22.18 | 25.46 |
| companion_regular | near_1.1R | constant | 31.66 | 23.02 | 27.34 |
| companion_regular | near_1.1R | stellar_density | 31.66 | 22.69 | 27.90 |
| companion_regular | near_1.1R | baryonic_gravity | 31.66 | 22.51 | 23.62 |
| standard_flat_FLRW | distant | constant | 23.51 | 23.51 | 27.27 |
| standard_flat_FLRW | distant | stellar_density | 23.51 | 23.08 | 30.31 |
| standard_flat_FLRW | distant | baryonic_gravity | 23.51 | 23.28 | 29.69 |
| standard_flat_FLRW | near_1.1R | constant | 26.18 | 23.90 | 26.85 |
| standard_flat_FLRW | near_1.1R | stellar_density | 26.18 | 23.31 | 30.09 |
| standard_flat_FLRW | near_1.1R | baryonic_gravity | 26.18 | 23.64 | 29.82 |

All population cases, parameter boundaries, per-target errors, original supplied inventories, blocked fractions, refined curves and sensitivity checks are in saturation-results.json. Grid bounds were not expanded after seeing results. A lower bound-component error achieved by leaving energy unbound is not by itself a successful total-gravity prediction.

## Outcome

None of these capacity prescriptions beats the previous partial-migration model on the retained-geometry omitted-galaxy comparison. Distant-input errors are about 24.64, 28.08 and 25.46 points for constant, stellar and gravitational capacities. The prior partial-migration value is about 20 points. Grid-matched comparisons to those archived predictions are included in the JSON; small changes from earlier numbers reflect radial sampling.

The best distant-input shared fits choose high enough capacities to behave essentially like unsaturated contraction. Thus their improvement over no migration is not evidence for saturation. The near-source gravitational-capacity branch does actively redistribute packets and improves its own baseline to 23.62 points, but still does not outperform the previous partial model. These are already-inspected data and the capacity amplitudes are phenomenological.

Some stellar-capacity omitted fits fail to bind about 28% of one target inventory. That energy remains in the traveling account; its gravitational field cannot be ignored in a physical model. Source ordering can change a selected bound cumulative profile by as much as 43.8 percentage points, much larger than the radial-grid effect (under 0.8 points). This is a serious assembly-history dependence of the allocation prescription, not numerical proof of a unique stable halo.

The next physical requirement is a time-dependent transport and source-history rule allowing occupied states, incoming packets and traveling overflow to interact consistently. The present test gives no basis to claim that a universal saturation threshold has solved the halos, or to modify the one-third capture exponent.

## Numerical and physical limits

Packet accounting closes to 2.22e-16; no bound cell exceeds its capacity. Doubling radial cells changes a selected cumulative curve by at most 0.00783, and reversing packet order changes it by at most 0.438. These differences quantify unresolved numerical or assembly-order dependence rather than disappearing into a fit. The shared/omitted scores in the table use the declared 160-cell grid; refined scores are archived.

The bound-plus-traveling packet energy is conserved. Gravitational work released during settling, pressure/support, momentum exchange, self-gravity and time-dependent supply still require a physical completion. Capacity does not switch gravity off: filled bound cells and traveling energy remain gravitational sources. No research goal, new observational lens fit, or derivation of saturation is claimed complete.
