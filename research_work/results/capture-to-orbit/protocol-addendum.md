# RB-1 consistency revision: one incident field, anisotropic drag, closed ledger

Declared before execution, 13 September 2026, after review of `main` at 21147a1. This revises how the archived RB-1 candidate is evaluated. It does not add a capture branch. The archived results stay in place as the failed first evaluation.

## Findings being corrected

1. **Mismatched spectra.** Bound production was validated with a smooth spectrum extending below threshold. The efficiency and drag constants used a spectrum cut at threshold. A receiver moving at speed beta Doppler-shifts companions, so a threshold-cut spectrum roughly halves slow-product production (review controls: 0.514 at beta=0.01c with v_esc=0.03c; 0.501 at 200/500 km/s).
2. **Drag assumed isotropic.** The ledger used the isotropic coefficient everywhere, although the attenuated field I = exp(-tau(r, mu)) is anisotropic.
3. **Open ledger.** The escaping energy omitted the receivers' own energy change.

## Revised evaluation

- **One incident distribution per control.** Define I(E, n) = S(E) A(r, n): S is a flat number spectrum per unit bath-frame energy, and A is the reference attenuation. The same I(E, n) feeds bound production, absorbed power, receiver force, receiver energy change and recoil. Two labeled controls share the same density per unit energy:
  - S1, threshold-cut: E in [E_th, 2 E_th].
  - S2, extends below threshold: E in [0.5 E_th, 2 E_th].

  Results for S1 and S2 are reported separately and never combined.
- **Bound production.** Computed by deterministic quadrature over receiver-frame product speed and incident direction, at the actual receiver speed and local escape speed. The isotropic, smooth-spectrum case must reproduce the uniform-ball rate v_esc^3/3 per unit spectral density. Exact-kinematics Monte Carlo checks both spectra at scaled speeds.
- **Receiver force.** First order in beta: f = Sigma_P [F_A - (1+chi) P_A.beta - u_A beta], with u_A, F_A and P_A the energy density, flux and pressure tensor of A, and chi the power-weighted d ln sigma/d ln E. Tangential motion therefore feels kappa = 1 + (1+chi) p_perp, with p_perp = P_perp/u_A. This formula is established radiation-force mathematics generalized to the threshold cross section; its use for companions is hypothetical. It is validated against exact reaction-weighted kinematics on four control fields:
  - isotropic
  - an isotropic bath boosted uniformly relative to the receiver
  - two fields with equal u and F but different pressure tensors

  The attenuated field keeps its second moments, and each receiver site uses its own p_perp.
- **Closed ledger.** Per retained rest energy, per control:
  - absorbed incident energy
  - bound rest plus kinetic energy
  - receiver energy change, from the same force
  - escaping energy = absorbed - bound - receiver change

  An aggregate test with exact kinematics must close channel by channel (bound, escaping, receiver) at the event level.

## What stays fixed

The receivers, potentials, orbit machinery, fixed retained inventory, frozen comparison rows and the archived population predictions. Only a per-site bound-production factor could change the population shape. If that factor varies across sites by less than 1% within a control, the archived shape is reused and the variation is reported; otherwise the populations are recomputed. Nothing is fitted.

## Wording correction

The first report's general statement, that an external isotropic bath cannot build an extended reservoir through ordinary-matter receivers, went beyond what one tested interaction establishes. The revised conclusion is limited to the tested RB-1 prescription.

## Promotion

No promotion is possible from this revision; it only corrects RB-1's evaluation. The self-illumination pilot is declared separately before it runs.
