# Growing collecting area with finite source energy

## Result

Coupling the calculated area growth to a shared finite radiation supply changes capture timing, not the total available energy. For three conversion-to-initial-capture rate ratios, growth shortens the time to capture 90% of an emitted pulse by about 2.3%, 19.3% and 30.9%. Both fixed-area and growing-area models ultimately deposit the same one injected energy unit.

| alpha / beta_initial | Fixed-area time to 90% | Growing-area time to 90% | Time reduction |
|---|---:|---:|---:|
| 0.1 | 24.584 | 24.021 | 2.3% |
| 1 | 4.421 | 3.569 | 19.3% |
| 10 | 2.949 | 2.039 | 30.9% |

Time units are 1/(beta_initial c), not years or an assumed universe age. beta_initial=n_receiver sigma_initial. These are synthetic supply histories and rate ratios, not measured cosmic parameters. The slow-conversion case gains little because photons have not yet supplied enough companions for a larger area to collect.

## Coupled equations

Use a homogeneous population of identical receivers and express photon energy X, traveling companion energy Y and deposited energy-equivalent mass m in common per-receiver units. Let h(m)=sigma(m)/sigma_initial be interpolated from the preceding A=1 immobile-storage calculation, at m=0,0.1,0.5,1. The known coupled transport equations become

    dX/dtau=j-a X
    dY/dtau=a X-h(m)Y
    dm/dtau=h(m)Y,

where a=alpha/beta_initial. Use j=1 for 0<=tau<=1 and zero thereafter, with initially empty radiation/deposit reservoirs. The explicit ledger is

    X+Y+m=integral j dtau=min(tau,1).

Increasing h drains traveling companions faster. It cannot supply extra energy after the source is exhausted. The last saved time is tau=300, by which deposited energy is within 1e-8 of unity in all cases. A continuing source would change the outcome but still has to pay for all stored energy.

## Inherited assumptions

The area curve inherits the spherical Newtonian immobile-storage model; support, gravitational binding energy and stress corrections remain absent. The per-receiver energy ledger is therefore not a complete relativistic gravitational-energy calculation. This does not repair the prior support problem.

The population is homogenized, radiation is isotropic and the normalized spatial capture profile is assumed to depend on accumulated mass, not on temporal bath fluctuations. Transport within each receiver is quasi-static. Real source/receiver correlations, propagation delays, evolving baryons and competition among unequal receivers require a more general model. The finite-area calculation's radial and boundary uncertainties are retained. This is a conditional finite-supply extension, not a self-consistent universe or observational halo fit.

## Verification

An implicit Radau solver integrates source-on and source-off intervals separately. The initial explicit-solver sensitivity run produced a small negative tail (-1.2e-9 in energy units) for one linear-interpolation case; it was not accepted as a physical solution. The final implicit runs satisfy the nonnegativity tolerance, monotone deposited energy and full injection ledger checks without clipping the radiation states. Only the area lookup is restricted to its computed mass interval [0,1].

For fixed area, an independent augmented matrix exponential verifies four times at each conversion rate to absolute 1e-9. Replacing shape-preserving cubic area interpolation with linear interpolation changes the growing-area 90% times by -0.0052%, -0.0864% and -0.1816%. This numerical sensitivity does not quantify physical uncertainty in the area law. Input hashes, area knots, all six primary trajectories, interpolation comparisons and ledger errors are retained.

The conclusion is conditional but useful: positive capture feedback and finite energy conservation can coexist. Larger collecting area can accelerate acquisition, while permanent deposits cannot exceed the supplied energy. Actual source histories, field support and observational motion/lensing predictions remain necessary. All six objectives remain open; no final holdouts were opened.

Run `python research_work/results/cluster-acceleration-capture/finite-bath.py` after feedback.py.
