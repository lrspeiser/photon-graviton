# Annulus sampling correction: implementation and verification protocol

Declared against main at 72aef4456aae63bbb3482706a5ecda6af77422c8, before executing the new checks.

## Scope

Correct the phase-space measure used to sample an already constructed WarmAnnulus. Add a new, opt-in sampler and a WarmAnnulus subclass for follow-on work. Leave the historical equilibrium.py, rut6.py, all previous protocols, results, and manifest entries unchanged so the archived stage-6 experiment (including failed H6) remains reproducible. This is a numerical correction, not a new force law or a stability result.

The independent package lives under research_work/annulus_sampling, outside the frozen results tree. Its checks do not require the user's local databases, network access, or external data. Registration in the full local research suite and integration with database-backed work follow the handoff; no full-suite success is claimed by the small package's checks.

## Measure

For L = R*v_theta:

    dM = f(E,L) R dR dtheta dv_r dv_theta
       = f(E,L) dR dtheta dv_r dL.

The R factor cancels. Positive-v_r quadrature nodes represent both velocity signs; angle is uniform. Nodal masses are 4*pi*dR*dL*w_vr*f, with radial trapezoidal weights. Preserve the archived density routine's equal rectangular L weights rather than silently changing its quadrature. The deterministic radial marginal must agree with 2*pi*R*Sigma(R)*dR, and the total must agree with integration of that same density. State clearly that matching this discretization is not proof of its continuum accuracy or global finite-mass support.

## Verification before use

1. Radial marginal: maximum absolute probability difference <= 1e-12 against an independently called density routine on the same quadrature.
2. Total mass and radius moments: relative error <= 1e-12 against the density integral (use an absolute criterion for an exactly zero moment).
3. Positive control: deliberately multiply probabilities by R in an extended analytic fixture. The verification must reject it, and its mean-radius shift must reproduce Var(R)/Mean(R).
4. Monte Carlo smoke test: fixed-seed draws reproduce the deterministic mean R, mean L, mean energy and mean v_r^2 within six predicted standard errors. Repeated seeds reproduce identical arrays, and recovered Cartesian L and v_r^2 match the selected nodes. This is a sampler test, not an equilibrium or stability test.
5. Input validation: nonfinite, nonpositive or nonmonotone radii, invalid widths, malformed source arrays, an unsupported nonuniform L grid, zero/nonfinite probability mass and invalid particle counts fail explicitly; sampling must not mutate source fields or parameters.
6. Integration check: test the adapter on the existing WarmAnnulus code without modifying it. Optional full four-annulus audits re-solve the old matched-support examples, report achieved mass/coupling/support, and refuse to hide construction failures or reduced sample counts.

No thresholds will be loosened to obtain a pass. A defect found before release is fixed and checks rerun; any scientific result remains not evaluated.

## Status and output policy

Report numerical verification separately from scientific conclusions. These checks exit nonzero on a numerical failure even if an old output would reproduce. A diagnostic report is written only to an explicitly requested new output path; existing outputs are refused. Historical H6 remains failed and is not replaced by this protocol.

## Remaining local work

After sampler verification: explicitly define bound-energy/domain cutoffs; independently refine the equilibrium integrals and domain; re-run corrected cold/warm populations with matched live/frozen initial states, common observation times and coherent-streaming diagnostics; then perform body-count scans and warm-population mode analysis. Fixed-support examples are not fixed-mass/fixed-coupling temperature experiments. No new claim of stability is licensed by correcting a sampler.
