# Orbital-condition-constrained combined fits

13 September 2026. Written before executing this refit.

Keep all data, geometry, capture parameters, lens-calibrated stellar masses, gradient form and radial anisotropy form from gradient-orbits-protocol.md. Fit only inner stellar motions; report the unfitted conditional outer bin. Keep the same numerical grids and parameter bounds.

Enforce gamma(r)>=2 beta(r) on 4097 logarithmic points, 1e-6<=r/Re<=100, using the same Abel-deprojected tracer profiles, quadrature order 512. Also enforce the exact central asymptotic limit: for positive Sersic components with at least one n>1, the steepest central cusp has gamma(0)=1-1/n_max. This follows from the Abel integrand near zero, nu proportional to r^(1/n-1); n=1 is logarithmic and subdominant. Thus beta0<=gamma(0)/2 is added. This extrapolation assumes the adopted central light model persists to zero and is not a direct central measurement.

Write f=r^2/(r^2+Re^2). At each outer beta, enforce beta0<=min[(gamma/2-beta_infinity*f)/(1-f)] over the grid, also capped by 0.45 and the asymptotic limit. Parameterize beta0 between -2 and this upper bound with a unit-interval coordinate. This algebra introduces no extra physical freedom. Use the previous 20 starts plus the unconstrained combined optimum mapped into the allowed interval. Compare against the feasible prior radial-only solutions. Positivity of a distribution function and stability remain unproved: this is only a necessary condition for the specified separable augmented-density class, not a universal sufficient orbital constraint.

Report all fits, parameter boundaries, cost relative to unconstrained fits, central constraint and finite-radius margins. Verify the resulting margins on a doubled, independently spaced 8193-point grid, and compare original reproduction limits. Do not use the outer bins to select parameters. Retain failed or boundary-dependent outcomes; do not promote a new companion law on this exposed-data test.
