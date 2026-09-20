# SE-DF1: discrete control-volume energy balance

20 September 2026. Protocol92b376e; numerical executable8d4e178. All24
preregistered random-state fixtures pass. Maximum scaled local residual is
4.75e-11; maximum scaled region-balance residual is4.18e-11. Both finite-time
directional-difference steps, spherical/cubic masks, periodic cancellation,
face/divergence agreement and total-energy reconstruction meet their gates.
First results and source hashes are preserved in discrete-flux-v1/results.json.

The compatible outgoing face flux is

    F_j = -(Dplus_j A) dot Pi_plus - c_Q^2 (Dplus_j Q) dot K_plus
          + epsilon/2 [Pi cross S_plus + Pi_plus cross S]_j,
    K=P-lambda A-epsilon(curl A) cross Q, S=Q cross K.

Here plus means the neighboring site in direction j. The site energy assigns
each forward gradient square to its starting site, matching the production
Hamiltonian. Taking its derivative with the unchanged LR3 RHS gives minus
the backward divergence of this face flux. Summing over a selected set of
sites telescopes to its crossing faces, with outward power

    P_out=h^2 sum_(sites,j) (mask-mask_plus) F_j.

This addresses a specific missing diagnostic: continuum flux alone was not
an exact local balance for the production code's mixed centered-curl and
nearest-neighbor-gradient stencil. The identity now uses that stencil
directly. The local assignment of link energy is a bookkeeping convention;
the staircase boundary represents the chosen grid mask, not an exact smooth
spherical surface. Spatial convergence of flux measurements still matters.

No production evolution was changed or replayed. These are instantaneous
semi-discrete checks, not measured time-integrated outgoing energy. The next
step is to accumulate flux at every Runge-Kutta stage and verify that the
region's energy change plus the integrated outward power stays within the
declared tolerance. Angular-momentum surface flux requires its own derivation.

The unchanged model remains a finite internal-excitation subsystem without
ordinary-matter generation, physical fuel, joint light/matter dynamics or
galaxy/cluster validation. No dark matter, expansion or old-gravity prediction
is used. Hamiltonian calculus and discrete control-volume identities are
established mathematics; this implementation is not a historical novelty claim.
