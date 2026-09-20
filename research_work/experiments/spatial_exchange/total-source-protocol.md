# SE-T: put the generated-wave source in the full scalar equation

Declared 20 September 2026 before evaluating this diagnostic. Apply to every
completed SE-1 and SE-2 final archive, including controls and refinements. Derived
partial summaries can be regenerated; original states remain immutable.

At fixed canonical state, vary U uniformly and differentiate the coupling
energy, keeping the phi restoring potential separate. Define

    D_field = 2 a sum_all Pi^2 - 2 beta dot J_all
              + mass2 (chi r^2-k_U X r),
    D_matter = sum_i W_i [alpha(E_i+z |p_i|^2/E_i)+2 beta dot p_i]/dV.

Integrate both over volume, and separately integrate the X/Y subset D_wave.
The local scalar momentum equation is

    Pi_phi_dot = discrete_laplacian(phi) - omega^2 phi
                 - g (D_field+D_matter) - gamma phi_dot.

Report D_wave / D_total_coupling, absolute D_total_coupling, matter and other
field contributions, and the separate restoring integral omega^2 integral(phi)/g.
Do not label these mass, force, total gravity or lensing multipliers. Source
geometry, restoring terms and evolving propagation all matter.

Two checks: independently reconstructed full Hamiltonian under uniform phi
shifts +/-epsilon/g at epsilon=1e-7 and 5e-8, subtracting the exact restoring
derivative, must agree with D_total_coupling to 1e-7+1e-5 abs(D_total_coupling).
Integrating the implemented scalar RHS, adding back sponge loss and restoring
terms, must reproduce D_total_coupling to 1e-10+1e-10 abs(D_total_coupling).
Periodic finite differences have exactly telescoping total Laplacian here;
pre-boundary qualification is still required for physical interpretation.

Use independent full-grid spherical source weights for matter and energy
reconstruction. Include the emitter-angle combination when reconstructing mass.
Archive state/source hashes, all errors, counts and explicit partial status.
This is accounting of the candidate's equations, not an observational or
nonlinear stability gate. Credit Hamiltonian differentiation as established math.
