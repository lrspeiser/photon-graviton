# Singular-kernel integration benchmark

Before calculation: compare direct tensor Gauss source quadrature with exact radial integration around an interior test point, followed by angular quadrature. Use uniform parallel straight trajectories spanning a square transverse area[-1,1]^2 and longitudinal range[-1,1]. Evaluate at(0.137,0.219,0). No softening, removed central region, density fit or modified physics.

Use the known age-integrated Newtonian passage kernel. Direct source orders6,8,12,24,48,96 are retained. For the reference, split polar angle at rays through each square corner and integrate radially exactly, comparing64/128 Gauss nodes per angular sector. Require relative potential and force differences<1e-10. Independently check force=-gradient potential by centered finite differences at h1e-4 and5e-5, requiring refined discrepancy<1e-7 of max(norm force,0.01). These are code/derivation checks, not a companion-theory validation.

The benchmark tests an integration remedy for known straight streams. Real bar trajectories are curved and their direction-age map may fold. Applying a similar method there requires locating close passages and integrating the actual source measure; success here does not establish convergence of the Galactic force coefficients.
