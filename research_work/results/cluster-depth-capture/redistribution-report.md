# Deposition followed by radial motion

The push-only preceding turn verified that main was synchronized but made no new scientific progress. This continuation takes the next available physical step: calculating a stored distribution from a specified post-capture motion rule. All six research objectives remain open.

## Optional rule and scope

Assume captured energy creates massive, nonrelativistic particles at rest at the capture location. After creation they move collisionlessly in the fixed ordinary Plummer well. At-rest birth gives zero angular momentum, so trajectories run through the center and return. Random spatial directions of birth give random orbit axes; uniform orbital phases give the phase-mixed limit. This is a deliberately explicit comparison branch, not a derivation of capture microphysics. The receiving system must absorb the incoming momentum and account for the energy needed for this initial state; that step has not been modeled. No claim of a complete energy-conserving conversion process follows from the orbit calculation.

The birth-radius distribution is the previously calculated depth-capture distribution, w(s)=4 pi s^2 q(s)/sigma. It is not fitted to a desired halo. No universe age is imposed. Full phase mixing is a limiting assumption; at any finite age, sufficiently distant, long-period orbits will not have sampled all phases. A finite-history model must replace this limit before comparison with data.

## Known orbital mathematics applied to our birth rule

Use GM=a=1, Psi(r)=1/sqrt(1+r^2), and let s be the birth radius, also the orbit's apocenter. Conservation of orbital energy gives

    v(r|s)=sqrt(2[Psi(r)-Psi(s)]),
    T(s)=integral_0^s dr/v(r|s).

T is time from apocenter to center. The fraction of time spent within radius R is

    F(R|s)= integral_0^min(R,s) dr/v(r|s) / T(s).

The accumulated mass fraction is then

    M_mixed(<R)/M_total=integral_0^infinity w(s) F(R|s) ds.

All weights and residence times are nonnegative. This supplies an explicit anisotropic orbit population and therefore does not require the negative isotropic populations found in support-report.md. It does not establish collective stability or the formation process. These are standard Newtonian orbit and probability identities; the model-specific assumptions are the companion birth distribution and at-rest capture rule.

## Results

| Capture strength A | Birth half-mass radius / a | Mixed half-mass radius / a | Reduction |
|---:|---:|---:|---:|
| 0.1 | 2.30862 | 1.56622 | 32.2% |
| 1 | 2.67757 | 1.83684 | 31.4% |
| 3 | 3.31865 | 2.31296 | 30.3% |
| 10 | 4.63936 | 3.29593 | 29.0% |

At A=1, the mass fraction within a rises from 0.1290 at birth to 0.2991 after phase mixing. At A=10 it rises from 0.0110 to 0.1273. Thus capture concentrated outside the center need not imply that the persistent gravitational source remains concentrated there.

The formal fixed-potential central limit is rho_mixed(r) proportional to r^-2, with enclosed mass proportional to r. Each radial orbit crosses the center at finite speed, giving finite probability per unit radius; dividing by a shell's r^2 area gives the cusp. The coefficient is a positive integral over w(s)/[4 pi T(s) v(0|s)], evaluated by the script. This is not an observed-core prediction: for any finite deposit mass, its central acceleration eventually dominates the ordinary Plummer acceleration, invalidating the fixed-potential approximation there. A self-consistent calculation or nonzero angular momentum is required before interpreting the central density or its lensing signature physically. The prior frozen-profile lensing curves cannot simply be retained for this branch.

## Verification and limits

redistribution.py writes redistribution-results.json. The substitution r=s sin(theta) cancels the turning-point singularity exactly. Capture weights integrate to unity within 3.2e-8. Doubling the orbital quadrature from 64 to 128 nodes changes the tested mass fractions by less than 1e-7 (recorded differences are below 9e-14). The predicted linear central enclosed-mass asymptote agrees at r=0.001a to within 2.7e-6 relative.

Independent differential-equation integrations from apocenters 0.3a, a, 3a and 10a reproduce falling times and three residence fractions per orbit within 5e-12. Specific orbital energy drift is below 4e-11 relative. These checks validate motion after the assumed initial capture, not the capture energy and momentum ledger.

The result advances the storage calculation: the birth profile now leads to a distinct calculated mass profile through a concrete motion rule. It is still a tracer-limit diagnostic. Capture recoil, self-gravity, angular momentum, finite-time evolution, stability, and comparison with real motion and lensing remain outstanding. No observational holdouts were accessed.
