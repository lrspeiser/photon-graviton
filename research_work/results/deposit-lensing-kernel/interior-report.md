# Interior-ray lensing benchmark

The previous goal turn made progress by checking an exterior-ray connection between the deposit force solver and conditional weak-field lensing. This follow-up covers light paths through a source, where the projected point-source kernel is singular but the extended-source deflection is finite.

## Calculation and independent reference

Use a uniform cube [-1,1]^3 with total mass one, represented by six tetrahedra of mass 1/6. Set G=c=1 for numerical coefficients. The same explicitly conditional equal-potential metric as in the earlier report gives direction change 2 times the transverse acceleration integrated along the ray. Its physical applicability to companions remains unproved.

For an independent calculation, integrate the volume's point-source ray kernel analytically over source height, then use polar coordinates in the projected square, centered at the ray's impact position b. This removes the apparent projected singularity with the polar area factor. For symmetric ray endpoints at z=+-L, L>1, cube density 1/8, the result is

\[
\Delta\boldsymbol n_\perp=
\frac12\int_0^{2\pi}\boldsymbol e(\theta)
\left[H(r_{\max},L+1)-H(r_{\max},L-1)\right]d\theta,
\]

\[
H(r,a)=\frac12\left[r\sqrt{r^2+a^2}+a^2\operatorname{asinh}(r/a)\right].
\]

Here r_max is the distance from b to the square boundary along the unit vector e(theta). Split the remaining angular integral at the four corner directions. These are conditional integrations of the known Newtonian kernel and known weak-field lensing relation, not new physical equations or a theory of energy conversion.

## Results

Impacts (0.137,0.219), (0.8,0.7) and (0.99,0.97) are each tested with L=5 and L=20. All six comparisons pass the declared 1e-7 relative criterion. Maximum disagreement between the tetrahedral-force ray integral and independent polar-volume result is 2.00e-15. Tightening the polar reference tolerance from 1e-8 to 1e-11 changes results by at most 4.79e-16. These errors characterize this benchmark; they are not observational error bars or universal accuracy bounds.

The cube ray integration explicitly splits at its outer faces z=+-1. No mass is removed or softened. No cell gain or physical parameter is fitted.

## Remaining scope

This verifies six finite, interior cube rays, including near-edge impacts. It does not prove convergence for the folded trajectory-generated density, exactly edge-coincident rays, arbitrarily thin cells, time-varying fields, shear or observed cluster lensing. The companion metric, source energy normalization and self-gravitating evolution remain unresolved. All nine goals remain open.

The two full-sphere trajectory preparations and their waiting comparison coordinator remain separate live work. This benchmark does not restart, replace or certify them.

Run `python research_work/results/deposit-lensing-kernel/verify_interior.py` from the repository root. All inputs and outputs are retained in `verify_interior.py` and `interior-results.json`. The [earlier report](report.md) gives the standard-lensing provenance and explicit metric assumptions.
