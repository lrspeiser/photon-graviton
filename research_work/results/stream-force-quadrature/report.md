# Close-passage integration without modifying gravity

This numerical benchmark holds the Newtonian force and a smooth source fixed. It demonstrates why increasing a direct source grid can give a stable potential while its force remains inaccurate, and tests an integration method that preserves the original force law. It is not a new companion model or a Milky Way fit.

## Geometry and known formulas

Uniform parallel trajectories span the square transverse area[-1,1]^2 and longitudinal interval[-L,L], with L=1. Their age distribution is uniform. Evaluate in the middle plane at transverse position(0.137,0.219). Divide all quantities by G times total source mass. If b is the transverse separation vector and r=|b|, the known time-integrated Newtonian kernels are

\[ K_\Phi=-\operatorname{asinh}(L/r)/L,\qquad \boldsymbol K_a=\frac{\boldsymbol b}{r^2\sqrt{r^2+L^2}}. \]

The direct method applies tensor Gauss quadrature in the transverse square. The second method centers polar coordinates on the test point, splits angles at the four corner rays and integrates the radial coordinate exactly. With the square boundary at radial distance R(theta), its radial primitives are

\[ P(R)=-\frac{R^2\operatorname{asinh}(L/R)+L\sqrt{R^2+L^2}-L^2}{2L},\qquad \boldsymbol A(R,\theta)=\boldsymbol e_r\operatorname{asinh}(R/L). \]

Here e_r is the outward unit vector(cos(theta),sin(theta)) along the ray. Integrate these expressions over theta and divide by area4. These are standard coordinate substitution and Newtonian integration, derived for this benchmark; no novelty claim applies. The integrable singularity is included analytically, not softened or cut out.

## Results

| Direct order per axis | Nodes | Potential error (%) | Vector-force error (%) |
|---:|---:|---:|---:|
| 6 | 36 | 1.286221 | 335.978 |
| 8 | 64 | 1.579823 | 383.175 |
| 12 | 144 | 0.115539 | 58.190 |
| 24 | 576 | 0.079098 | 17.449 |
| 48 | 2304 | 0.011015 | 14.618 |
| 96 | 9216 | 0.003429 | 8.690 |

The angular64/128 reference changes by7.39e-15 in potential and6.31e-15 in vector force, relatively. Its force is[-0.07079463149110038, -0.11522296313787489]. An independent finite-difference check of minus the potential gradient agrees to6.14e-10 relative scale. All declared reference checks pass. Direct-grid errors are measured against that verified numerical reference, not an astronomical measurement.

## Consequence for the actual calculation

The real captured trajectories are curved, with a nonuniform direction-to-position map that can fold. The square-stream primitive cannot simply replace their field integral. The next method must locate close passages and account for the actual source measure, treating the difficult local part analytically or with appropriate local coordinates, while integrating the remainder and verifying convergence. Do not declare a Gaussian grid sufficient because the potential looks stable.

This provides a tested numerical approach without adding screening, softening, a velocity distribution or a new physical parameter. It does not resolve capture energy, self-gravity, redshift or lensing. All nine project goals remain incomplete. Reproduce with run.py; protocol.md records the criteria before calculation.
