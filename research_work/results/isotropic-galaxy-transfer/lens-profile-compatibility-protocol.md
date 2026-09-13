# Lens/profile compatibility diagnostic

Declared before execution. Keep exact q=1/3, all capture constants, optical geometry, photometric proxies and the total deposit inventory fixed. Reuse the six lens systems and both population proxies. This is a deliberately target-fitted inverse problem, not a predictive model or a newly successful test.

Allow nonnegative mixtures of seven conservative dilations s=[0.1,0.25,0.5,1,2,4,10], with weights summing to one. Each component has M_s(<r)=M0(<r/s), so each has the original total inventory. Finite mixture profiles are smooth and positive. This spans more possibilities than previous shared two-parameter migration, but is not every possible profile.

Fit all observed stellar bins with their full covariance. Compare (a) original profile with stellar mass/beta free; (b) original profile with its lens equation forced to the catalog angle; (c) mixture with that same exact-angle constraint. The lens equation at fixed angle is linear in stellar mass and mixture weights; eliminate stellar mass analytically. Preserve the existing broad stellar mass and beta bounds. Multi-start the simplex optimization and check the exact coefficient calculation after interpolation-assisted fitting. Record mass and beta boundaries.

The catalog angle has no full source-image likelihood here. Exact agreement is a diagnostic constraint, not a precision measurement or a prediction. There is no calibrated compatibility p-value; compare motion residuals and maximum standardized marginal residual, while retaining the full covariance score. If flexible profiles still fit poorly, that constrains only this basis and these fixed inputs, not all redistribution or all gravity laws. A successful flexible fit establishes only algebraic/Jeans compatibility, not a physical phase-space distribution, stability, work budget or universal migration law.

Report per-system weights and cumulative mass changes at 0.5,1,2,5 effective radii so the required locations are explicit. No universe age, size, external amplitude, one-third exponent, or light-bending multiplier is adjusted.
