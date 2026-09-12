# Quadratic cell geometry: early improvement, late failure retained

The completed second mesh still failed 14 of 30 spatial force gates. This test asks whether a curved positional representation can improve the source map using the actual trajectories already available, before committing to another gravity calculation.

## Known mathematical construction

For each first-nested triangle, use its three vertex trajectories and three actual edge-midpoint trajectories from the second mesh. Let b_i be the barycentric coordinates obtained by intersecting the launch ray with the parent triangle plane. The positional approximation is

\[
\boldsymbol X(b,t)=\sum_{i=1}^3 b_i(2b_i-1)\boldsymbol X_i(t)
+4b_1b_2\boldsymbol X_{12}(t)
+4b_2b_3\boldsymbol X_{23}(t)
+4b_1b_3\boldsymbol X_{13}(t).
\]

These are standard degree-two Lagrange triangle shape functions, not new physics or a novel formula; see [DefElement's explicit basis](https://defelement.org/elements/examples/triangle-lagrange-equispaced-2.html). The implementation checks nodal identity, partition of unity at the declared test coordinates, and correspondence of each midpoint to its parent edge. All trajectory caches are hash-checked. No new orbit integration or fitting is involved.

The coefficients can be negative. They interpolate positions; they must not be reinterpreted as positive source-mass weights. A gravity implementation would require a separate mass-preserving integration over the reference coordinates, accounting for folding. None is claimed here.

## Same-information comparison

Compare against the piecewise-affine second mesh, using the same actual trajectory information and the same exposed source8 reference paths/weights. RMS errors below are percent of launch radius.

| Source radius | Age | Piecewise affine | Quadratic |
|---|---:|---:|---:|
| 1 kpc | 0.05 | 1.45 | 0.69 |
| 1 kpc | 0.1 | 6.13 | 5.34 |
| 1 kpc | 0.25 | 23.35 | 25.63 |
| 3 kpc | 0.05 | 1.45 | 0.68 |
| 3 kpc | 0.1 | 4.23 | 2.54 |
| 3 kpc | 0.25 | 40.26 | 40.52 |

Age is in kpc/(km/s); the final age is approximately 244 million years. Early-time curvature is better represented, but the late-time source is not rescued. At final age, 93.07% of inner reference weight and essentially all outer reference weight still exceed the 2%-of-radius diagnostic. Maximum final position errors rise to 1.271 and 2.681 kpc.

The targeted inner probes improve to [2.32%, 7.40%, 2.36%, 7.05%], while outer probes give [88.31%, 16.53%, 42.55%, 74.46%]. Their mixed behavior and the worse overall final RMS remain in the outputs. They are development probes, not astronomical holdouts.

## Decision and remaining work

Do not adopt this quadratic map as a successful cure for the late-time density or run a gravity calculation expecting it to establish convergence. It identifies a useful early-time representation but shows that the late map needs additional resolved trajectory information or a different controlled representation; simply increasing interpolation order on this mesh is insufficient.

No gravitational coefficient, source amplitude, capture postulate or physical field equation changed. No mass was reassigned. No new force prediction follows. All nine goals remain open, including physical energy funding, self-gravity, companion metric and joint observed redshift/lensing tests.

Run `python research_work/results/bar-cell-geometry/check_quadratic_geometry.py` from the repository root. All seven ages and every reference/probe error are saved in `quadratic-geometry.json`. No pipeline process is running from this test.
