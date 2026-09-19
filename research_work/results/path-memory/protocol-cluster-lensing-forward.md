# CL-F1: cluster optics from the frozen written-field response

Declared 19 September 2026, before running this experiment. Base: merge
`37dab0641ddd0af99b762d4a23cf24e8efba0e79`. This is a new, bounded follow-up;
the six-item review mentioned in the handover was not present in the checkout,
so this is not represented as its completed CL-2 stage 2 protocol.

## Hypothesis and boundaries

Use this project's existing written field, C = sum Lambda (rho_b * K_w), and
its explicitly conditional light rule L1 (Phi = Psi = Phi_b - C). Use only
ordinary gas and stellar sources. No halo, MOND interpolation, fitted extra
mass, object-type switch, or new light-bending multiplier is introduced.
Lambda and w are copied from CL-2's **clusters-alone pressure fit**, frozen
before any optics calculation. A mathematically solvable lens is not a
validated explanation of cluster observations. L1, formation, photon funding,
time evolution, and cross-scale agreement remain open.

## Six calculations

1. Reproduce CL-1 and CL-2's existing check jobs, preserving their archives.
   Pin the input and response archive hashes. Keep stage 8's failed status
   and stage 9's unrun declaration unchanged.
2. Project each of the 12 X-COP baryon sources as nonnegative, piecewise
   constant-density spherical shells. Each shell's projected mass and surface
   density are analytic. Convolve that projection with the frozen Gaussian
   widths using CL-1's independently checked ring kernel. The source ends at
   CL-2's r_out; that boundary is explicit, not an inferred invisible envelope.
3. Export physical bend alpha, projected baryon mass, signed equivalent field
   mass, Sigma and DeltaSigma over 1 kpc to max(5 r_out, 8 w_max). Equivalent
   field density is a diagnostic of the potential, not particle matter;
   negative values must not be clipped. Compute kappa, shear, reduced shear,
   both Jacobian eigenvalues, signed magnification, and bracketed critical
   radii for a separately labeled **fictional geometry demonstration**:
   D_l = 200 Mpc, D_s = 1000 Mpc, D_ls = 800 Mpc. These are not adopted
   distances for the actual clusters. Export geometry-free curves so users
   can supply measured/adopted distances without changing the field.
4. Solve the signed spherical lens equation for a point source displaced by
   5 arcseconds in that demonstration. Return roots, residuals, parities,
   and search bounds; a sign-bracket search cannot certify tangent roots,
   images below the minimum radius, or images outside its range. Aligned
   sources are handled as critical rings, not a finite image list.
5. Transfer the same frozen spectrum to CL-1's two Coma baryon brackets.
   Score the six already exposed weak-shear bins using one nonnegative
   inverse-critical-density nuisance per bracket, exactly as a shape test.
   Report baryons-only alongside. No absolute Coma amplitude or new holdout
   is claimed: source-distance information is still missing. Report the
   reduced-versus-weak correction as a diagnostic, not a second fit.
6. Save strict JSON, per-cluster CSV, a plot and a readable report in a fresh
   generated directory. Add a portable numerical check job without touching
   the frozen CL-1/CL-2 results.

## Numerical gates (fixed before execution)

- Analytic uniform-sphere projection, total mass, and zero-response limit:
  relative error < 1e-10. A doubled projected mass must fail the bend check.
- Gaussian ordinary-matter source plus Gaussian written footprint: analytic
  projected mass, bend and surface density vs the implementation, normalized
  maximum error < 2e-4; doubling source and projection grids must improve the
  error. A missing sqrt(2 pi) w projection factor must be rejected.
- Three-dimensional force integral vs two-dimensional projected deflection
  for a finite, smooth source and a mixture: normalized error < 2e-4.
- Lens-equation images and critical curves vs an independently analytic
  point mass: relative radius error < 1e-8; lens-equation residual < 1e-9
  arcsec. Invalid distances, negative shell mass and invalid widths rejected.
- Actual-source projection refinement (1500/4000 to 3000/8000 source/projected
  samples) on A1795 and A2319: bend and DeltaSigma normalized maximum change
  < 3e-3 over 10 kpc to the exported outer edge. Inner <10 kpc is explicitly
  extrapolated and not certified. Failure blocks corresponding numerical
  claims and is recorded without moving the gate.

No observational success gate is added: the deliverable is a verified forward
calculation and an honest transfer result, not a guaranteed empirical fit.
