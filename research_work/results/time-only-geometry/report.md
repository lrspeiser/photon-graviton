# Pure clock-rate geometry needs a gravitational completion

## Question and provenance

Can changing only local time rates, with exactly flat unchanged physical space, constant local light speed and zero shift, be sourced by the positive-energy radiation/companion candidate under ordinary Einstein gravity?

This is a conditional comparison, not an assumption that the final fictional-universe theory must use GR. The metric, curvature equations, Einstein constraints and null energy condition are established physics. Their application and symbolic verification here are not claimed unique. The standard 3+1 constraint framework is described in [Gourgoulhon, 3+1 Formalism and Bases of Numerical Relativity](https://arxiv.org/abs/gr-qc/0703035).

## Exact result for the planar diagnostic

Established metric representation applied to the stipulated restricted geometry:

ds^2 = -c0^2 q(t,x)^2 dt^2 + dx^2 + dy^2 + dz^2, q>0.

Established clock/ray relation: fixed-position observers have d tau=q dt, and null rays measured with their local physical rulers have speed c0. This implements constant local speed kinematically. It does not yet specify a source for q.

Computed here from the established curvature definitions, using c0=1:

G_tt=G_tx=G_xx=0,   G_yy=G_zz=(partial_x^2 q)/q.

The symbolic calculation retains arbitrary time dependence in q. No q_t term supplies normal-observer energy density in this exact geometry. The general 3+1 interpretation is that intrinsically flat spatial slices with zero shift and no spatial-metric evolution have both zero spatial curvature and zero extrinsic curvature. The established Hamiltonian constraint then sets the total normal energy density to zero when the cosmological constant is zero.

Thus this exact ansatz cannot be sourced by the positive total radiation-plus-canonical-field energy used in the prior Hamiltonian tests, unless additional gravitational/matter assumptions are changed. Treating radiation as a negligible test signal on a prescribed background is a useful approximation but cannot demonstrate a self-consistent source-funded geometry. The constraint concerns total gravitating energy, not the algebraic photon-to-field transfer ledger alone.

## An illustrative stress check

Proposed example profile for this diagnostic only, not a derived void solution or novelty claim: q=1-0.1 exp(-x^2), in dimensionless length units. Established Einstein normalization gives transverse pressure proportional to q_xx/q, with proportionality c0^4/(8 pi G).

At x=0 that dimensionless factor is +0.222222; at x=2 it is -0.0256889. Since the normal energy density is zero, a null direction along y or z samples a negative stress-energy contraction at the latter point. This violates the usual null energy condition there. It is a planar profile with infinite transverse extent, not an actual finite galaxy void. Its positive q has no coordinate horizon. The source has not been supplied.

A positive-energy minimally coupled canonical scalar and ordinary radiation satisfy the relevant null-energy inequality; simply naming them as the source does not generate this prescribed stress pattern. Nonminimal couplings or different gravitational equations require a fresh analysis, not automatic use of this conclusion.

## Consequence without importing expansion

The failure is specific to exactly flat fixed spatial slices, zero shift, ordinary Einstein gravity, zero cosmological constant and the proposed positive-energy sources. Nonexpansion does not imply flatness. A nonexpanding but spatially curved geometry may have a nonzero spatial-curvature contribution to the energy constraint and must be solved separately. General shifts, local spatial evolution, or an interaction affecting light differently from a universal spacetime clock also change the assumptions and require explicit treatment.

No expanding scale factor or expansion-based redshift formula was used. This calculation does not refute the user's broader premise; it identifies why kinematic clock/ruler ratios alone are insufficient. The source and gravity equations must make the geometry and measured clocks compatible.

## Verification and next work

The script computes all Christoffel symbols and Ricci/Einstein components directly using SymPy, checks the full matrix against the expected expression, and evaluates the profile stress signs. Run `python research_work/results/time-only-geometry/run.py`. This is an exact symbolic consistency check, not a numerical redshift fit or observational validation.

Next formulate the minimum curved, nonexpanding or explicitly modified-gravity completion with the same matter and companion action. Preserve fixed published distances as the user's adopted observations without equating them automatically to flat-coordinate path lengths in a curved model. Derive clocks, null propagation, energy support and lensing together. The existing empirical redshift fit, new-data quarantine, uncertainty model and all broader requirements remain unfinished.
