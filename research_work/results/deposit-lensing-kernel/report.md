# Conditional connection from deposited gravity to lensing

While the second full-sphere trajectory preparations continue, this calculation implements and checks the mathematical connection needed for goal 5. It does not replace the outstanding source-resolution tests.

## What must be specified

A Newtonian potential determines slow-star acceleration, but by itself does not specify all the geometry experienced by light. A proposed weak-field completion is

\[
ds^2=-(1+2\Phi/c^2)c^2dt^2+(1-2\Psi/c^2)d\boldsymbol{x}^2.
\]

This two-potential notation is known mathematics. The model still needs field equations specifying both potentials. In these conventions slow-star acceleration is approximately \(-\nabla\Phi\), while the null condition gives a coordinate optical index \(n\simeq1-(\Phi+\Psi)/c^2\). Therefore the first-order change in propagation direction is

\[
\Delta\boldsymbol n_\perp=-\frac1{c^2}\int\nabla_\perp(\Phi+\Psi)\,dz.
\]

For the **conditional equal-potential branch** \(\Psi=\Phi\), this becomes \(\Delta\boldsymbol n_\perp=2c^{-2}\int\boldsymbol a_\perp dz\). These are standard weak-field lensing relations, not novel equations. The static equal-potential metric and factor of two are described in [Bartelmann and Maturi, sections 2.1–2.2](https://arxiv.org/html/1612.06535). Only the local lensing relations are used here, without an expansion-based distance conversion.

The physical direction-change vector points toward a positive-mass lens. Some lens-equation conventions define the deflection vector with the opposite sign; our numerical output is the actual change along the incoming propagation direction.

## What was actually checked

Reuse the existing exact, unsoftened constant-density tetrahedron force solver. Integrate its transverse acceleration along four finite straight rays, with two impact positions and two segment lengths. Compare with a separate volume integral using positive tetrahedral Gauss quadrature and analytically integrated point-source ray kernels. No deposited galaxy distribution, fitted gain or observed lensing mass is used.

All four exterior-ray comparisons pass the declared 1e-7 relative tolerance. The largest disagreement is 1.12e-14; changing the independent volume quadrature from 16 to 32 nodes per transformed coordinate changes the result by at most 6.27e-10. Three finite-segment point-lens checks also pass, with maximum relative error 2.23e-16.

For the point lens, the known finite-segment result used as a check is

\[
\Delta n_x=-\frac{4GM L}{c^2 b\sqrt{b^2+L^2}},
\]

for impact parameter \(b>0\) and symmetric endpoints \(z=\pm L\). It tends to the standard inward deflection of magnitude \(4GM/(bc^2)\). Numerical coefficients set G=M=c=1; no physical source amplitude is supplied by that choice.

## Why this is conditional

Equal potentials are an explicit response assumption for this test, not a derived property of the companion/time field. A field with significant directional stresses, a modified coupling, or additional optical effects can require a different relation. Changing only clock rates does not automatically provide the same lensing as changing both potentials. The full model must derive its response and apply the same fields to clocks, stars, light and gravitational waves.

The tests cover exterior projected rays through a static weak field. They do not validate interior-ray singular integration, a moving or self-gravitating population, shear or time-delay predictions, strong-field propagation, physical energy funding, or any observed galaxy/cluster. The unconverged source map remains unconverged. All nine goals remain open.

## Reproduction

Run `python research_work/results/deposit-lensing-kernel/verify.py` from the repository root. `results.json` records all cases and errors. The new implementation is a conditional numerical connection to the existing gravity kernel; it is not a new first-principles theory of time or gravitons.

Subsequent check: [six interior cube rays](interior-report.md) now pass an independent polar-volume comparison. This extends numerical coverage while retaining the original metric and physical limitations.
