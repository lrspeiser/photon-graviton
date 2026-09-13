# Feedback from companion gravity into reservoir exchange

This test updates the exchange-driving Newtonian potential using ordinary stars AND the current compact and extended companion deposits. It also updates the occupancy penalty from the current density. The one-third capture law, original finite inventory and baryonic gravity prescription are unchanged. This is a coupled snapshot closure, not a simulated source history or complete relativistic calculation.

## What is coupled

Known spherical potential: Phi(r)=-G[M(<r)/r+integral_r^infinity dM(rprime)/rprime]. Stellar Phi uses the original Hernquist proxy. The prior proposed rate fractions are recomputed in the combined field; solve f=F[Phi(f),rho(f)] with damped iteration. Compact packets move to s times their ORIGINAL radius, not repeatedly inward on each numerical pass. The fixed-point iteration has no assigned physical duration.

Both populations conserve their combined mass-equivalent inventory. The baryonic density is held fixed; its orbits do not respond here. The rate laws and contraction distance remain phenomenological, and final support has not been solved. Thus this is more complete feedback than the initial-field test but still not everything responding dynamically to everything else.

## Fixed parameters first, then shared refits

Frozen tests preserve previous shared and omitted-target parameters for every geometry/population/source branch. Primary refits use Chabrier/distant inputs, separately for both geometries, on the declared 25x25 grid. Nonconvergent trial closures are excluded and their counts recorded. The targets have been previously inspected; five-to-one transfer is not blind evidence.

| Geometry | Rule | Previous omitted error | Frozen with feedback | Refit with feedback, omitted error |
|---|---|---:|---:|---:|
| companion_regular | binding | 20.34 | 21.69 | 20.04 |
| companion_regular | released_binding | 20.15 | 22.66 | 20.59 |
| companion_regular | occupancy_penalty | 21.85 | 25.63 | 20.86 |
| standard_flat_FLRW | binding | 24.16 | 26.96 | 24.56 |
| standard_flat_FLRW | released_binding | 24.14 | 28.35 | 24.51 |
| standard_flat_FLRW | occupancy_penalty | 23.47 | 27.74 | 22.93 |

Errors are RMS cumulative-profile differences in percentage points, not velocities, lens residuals or significance. Refit grid spacing differs from the older grid; small changes must not be overinterpreted. Standard geometry is a comparison only. All frozen variations, failed trial counts, parameters, per-target curves and convergence diagnostics are in the JSON.

## Traveling companions

The instantaneous traveling inventory is unknown. At frozen parameters, separate scenarios use traveling/stored ratios 0, 0.01 and 0.1, with a spherical profile proportional to capture density divided by opacity. These are sensitivity assumptions, not measured supply or inferred cosmic ages. Its prescribed energy density contributes to the Newtonian driving potential; total enclosed output adds the traveling component explicitly. It is never normalized away. Relativistic pressure, directional stresses, photon gravity and a covariant field equation are not computed, so this does not determine the true gravity of radiation-like companions. The baseline stellar proxy already includes stellar mass; no separate photon inventory is inferred from it.

## Outcome

With prior parameters frozen, deposited self-gravity increases retained-geometry omitted errors: binding 20.34 to 21.69 points, released binding 20.15 to 22.66, and occupancy penalty 21.85 to 25.63. Ignoring the deposits in the driving field was therefore a material approximation.

Refitting recovers approximately the previous best performance: binding gives 20.04, released binding 20.59 and occupancy penalty 20.87 points. The earlier partial-migration benchmark was 20.01. The small binding difference is below the profile-refinement scale and is not evidence of an improvement. Standard-geometry occupancy feedback gives 22.93 points versus 23.47 previously, but this comparison is not adoption of expansion or a new lensing success.

The shared retained-geometry binding fit uses s=0.1732 and v_ex=568.37 km/s, with compact fractions about 33.3-42.5%. All frozen runs and refit grid trials converged; all selected alternate seeds converged to essentially the same solution. This supports uniqueness within the tested seeds and parameter families, not a proof of global uniqueness or dynamical stability.

Prescribed traveling energy is not the missing cure in these scenarios. At the old retained-geometry omitted parameters, increasing its ratio from zero to 0.1 changes deposited RMS from 21.69 to 21.86 for binding, 22.66 to 22.67 for released binding and 25.63 to 24.71 for occupancy. These compare deposit response only; additive traveling gravity is separately recorded and remains a relativistic modeling limitation.

Feedback is now represented for the deposited component, but the hypothesis still lacks evolving baryons, physical arrival histories, angular/radiation stresses, and a demonstrated support-and-energy-loss mechanism. It has not become a complete all-energy gravitational model.

## Energy and verification

Potential energy before and after settling is evaluated for stars plus deposited self-gravity in the zero-travel case. The reported released specific binding energy is a required work/radiation channel, not extra stored mass or a completed energy budget. Changes in baryonic motion, support and radiation would have to close that ledger.

Direct shell potential checks pass. Turning companion feedback off reproduces prior curves within 0.00118. Refined selected profiles change by at most 0.001615; the largest difference between all-extended and all-compact seeds is 2.498e-08. A converged fixed point is numerical consistency, not a dynamical stability proof. Nonconverged and alternate-seed flags are retained. No new stellar-motion/lensing likelihood fit or completed research goal is claimed.
