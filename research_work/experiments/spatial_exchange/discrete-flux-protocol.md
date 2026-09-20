# SE-DF1: exact semi-discrete local energy flux

Declared 20 September 2026 before calculation. Do not alter the pinned LR3/LR4
equations. Their local energy density assigns forward-link gradient squares
to the lower-index site. With d_j=Dplus_j A, e_j=Dplus_j Q, S=Q cross K,
and a plus subscript indicating the next site in direction j, propose the
outgoing flux on that face:

    F_face,j = -d_j dot Pi_plus - c_Q^2 e_j dot K_plus
               + epsilon/2 [Pi cross S_plus + Pi_plus cross S]_j.

Nearest-neighbor Laplacian terms combine with time derivatives of forward
gradient squares to give a backward divergence. The symmetrized cross-product
face term exactly accounts for centered-curl exchange. Thus

    dH_site/dt + sum_j (F_face,j-F_face,j_minus)/h = 0.

For a binary region mask m, outward power is

    P_out = h^2 sum_i,j (m_i-m_i_plus) F_face,j.

This is a discrete control-volume identity for this particular allocation of
link energy, not a unique microscopic definition of energy location. It tends
to the continuum flux for smooth fields. Finite-difference product identities
and Hamiltonian calculus are established mathematics, not new gravity laws.

Declare24 fixtures: n8/12,L8,epsilon0/2,lambda0/1, three random states each,
seed20261005, independent normal components scaled.03. Use unchanged LR3 RHS.
Compute site energy at state plus/minus delta*RHS for delta1e-6 and5e-7.
Require both local continuity residuals scaled by max(1,max(abs(div F)))<1e-7.
Require both region balances for a centered sphere radius1.5 and cube with
max(abs(x))<2 to have absolute residual /max(1,abs(P_out))<1e-7. Compare
explicit crossing-face power against summed divergence to1e-12 scaled.
Periodic total face power must cancel to1e-12. Verify summed site energy
against the independently implemented LR3 Hamiltonian to1e-12 scaled.

This is instantaneous semi-discrete conservation, not a time-integrated
production flux measurement. A later evolution must accumulate face power
at integrator stages and test the full volume change plus integrated outward
power. No angular flux, physical source fuel or observational solution follows
from these checks. Existing live simulations remain unchanged.
