# SE-LF1: continuum local energy flux verified

20 September 2026. Protocol8a1f543; executablee89ca77. All24 resolved Fourier
fixtures pass. Maximum scaled local continuity residual is2.08e-12; maximum
absolute integrated time derivative is1.73e-12. Two directional-difference
steps were used, with source hash and first results in local-energy-flux-v1.

For K=P-lambda A-epsilon(curl A) cross Q, the conserved Hamiltonian's flux is

    F_i = -Pi dot partial_i A - c_Q^2 K dot partial_i Q
          + epsilon [Pi cross (Q cross K)]_i,
    partial_t H + div F = 0.

The direct coupling enters through K. Its lack of spatial derivatives means
it adds no separate gradient-flux term. The final cross-product term accounts
for the curl-dependent coupling and cannot be omitted from an energy ledger.
Canonical Hamiltonian calculus and this derivation method are established
mathematics; no historical novelty is claimed.

This identity allows a continuum surface-energy balance to be formulated.
It is not an exact local balance for the production code's mixed centered-curl
and nearest-neighbor gradient discretization. A compatible discrete face flux
and its volume/surface check are still needed before reporting transported
energy from those simulations. Existing production equations and live jobs
were not changed. No angular-momentum flux, material fuel rate, sustained
gravity or observational prediction has been established by these fixtures.
