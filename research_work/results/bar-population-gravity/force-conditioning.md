# Why local force needs more source resolution

The source6/source8 comparison already showed stable potentials but21 failed force gates. This diagnosis uses every retained contribution; none is removed, capped, softened or refitted.

## Measured cancellation and concentration

Let a_i be a direction's already source-weighted and age-integrated vector acceleration, including its four symmetry partners. Define

\[
C=\frac{|\sum_i\boldsymbol a_i|}{\sum_i|\boldsymbol a_i|},\qquad
D=\frac{\max_i|\boldsymbol a_i|}{|\sum_i\boldsymbol a_i|}.
\]

These are known sum-conditioning diagnostics, not new physical laws or uncertainty estimates. Small C means strong cancellation. D can exceed one without violating conservation: other directions oppose the largest contribution. A large D does not, by itself, identify a bad trajectory.

For source8 at the final epoch, the3kpc source's central diagnostic point (0.1,0,0) has C=0.1817 and D=1.0297. Its1kpc-source counterpart has C=0.6165 and D=0.2358. For the1kpc source evaluated above the plane at (1,0,1), C=0.9921 and D=0.0569. All directions and the three largest contributions at each point/time are retained in concentration.json. These statistics describe these quadratures, not an established continuum density.

Thus an accurately computed contribution can still produce an inaccurate total if the continuum's opposing directions are inadequately sampled. The diagnosis is consistent with the observed sensitivity near the center; it does not establish cancellation as the only source of all force discrepancies.

## Local passage calculation

**Known Newtonian kernel, derived here for a numerical example.** Approximate a short passage by a straight path X(t)=(b,0,vt), with closest approach b>0 to a test point at the origin and times -L/v to L/v. This approximation is only for diagnosing integration; it does not replace the bar trajectories. The time-averaged potential and transverse acceleration per G times particle mass are

\[
\langle\Phi/(GM)\rangle=-\frac{\operatorname{asinh}(L/b)}{L},\qquad
\langle a_x/(GM)\rangle=\frac{1}{b\sqrt{b^2+L^2}},\qquad
\langle a_z\rangle=0.
\]

They follow by integrating1/sqrt(bÂ²+vÂ²tÂ²) and b/(bÂ²+vÂ²tÂ²)^(3/2), respectively, and dividing by the passage duration2L/v. The acceleration points toward the path. For b much smaller than L, potential grows only logarithmically as b decreases, while force grows as1/b. Therefore point-force quadrature can be substantially harder than potential quadrature even with the same trajectories.

Independent adaptive time integration at b/L=0.1,0.01,0.001 agrees with both formulas to better than2e-16 relative error (passage-checks.json). This verifies the mathematical example, not a new physical law.

A continuum of nearby paths must be integrated over its actual distribution. For a locally smooth distribution in a two-dimensional impact plane, the area element b db dphi makes the1/b magnitude locally integrable; opposing directions can cancel further. Cold focusing or a singular source distribution requires a separate analysis. A finite direction quadrature is not a physical collection of massive isolated streams, so a close sampled path is not evidence for a new physical force spike.

## Next numerical action

The source12 calculation is complete and still fails19/30 force gates. Use local source-cell refinement or a singularity-aware integration scheme, retaining all angular cells and their weights. The straight-stream benchmark in ../stream-force-quadrature/report.md verifies one local-coordinate remedy, which still needs adaptation to the bar trajectories. Direct the error estimate at vector force, including cancellation; a stable potential cannot serve as its substitute. Validate any new scheme against analytic passage/continuum integrals and against the existing outputs before using it for galaxy fits.

Do not add a velocity dispersion, a softened physical force or permanent retention merely to make the numerical test pass. Those would be additional physical assumptions requiring their own rationale and observational predictions. No lensing or source-energy closure follows from the present diagnostic. All nine goals remain incomplete.
