# RB-1: outgoing waves can supply the radial shape, but not the budget

20 September 2026. Analytic necessary-condition screen of the weak, massless, eta=lambda=0 SR-1 branch. These are illustrative inputs, not fitted observations. The derivative coupling of wave energy to the scalar potential has coefficient 1+s relative to slow matter; it is not a new adjustable charge.

The time-averaged outgoing energy density is u=L/(4 pi r^2 c). Its gravity has the required asymptotically flat radial shape inside the wave front:

v_wave^2(r) = (1+s) G L/c^3 * (1-r0/r).

That shape alone is insufficient. Sustained luminosity must be L=v_flat^2 c^3/[(1+s)G], and its maximum fuel-limited lifetime is epsilon (1+s) G M/(v_flat^2 c).

## Results

- 36 independent shell-quadrature force/energy comparisons pass; maximum force relative error 2.22045e-16, energy error 3.33067e-16.
- 54 source-budget fixtures evaluated; 0 last the declared 10 billion years.
- Across fixtures, even the longest lifetime is 2821896.27 years, using all available rest energy.
- For s=1, v=200 km/s and M=2e41 kg: required luminosity 8.07395605e+45 watts. Converting the entire source rest energy lasts only 70547.4067 years. Allowing only 1e-6 of rest energy reduces that lifetime by 1e6.
- The same example would require an active-source coefficient 283497.31 to last 10 billion years with all rest energy, versus the derived coefficient 2. This is a diagnostic discrepancy, not permission to insert a fitted multiplier.
- Storing the wave population out to 6e20 m, with r0=6e18 m, requires 1.59975001e+58 joules, equivalent to 1.77996194e+41 kg or 0.889981 of the illustrative source rest mass. Trapping eliminates rapid escape but does not remove that positive-energy inventory.

## Decision and scope

Reject sustained, freely escaping weak-wave self-gravity as the long-lived explanation for every declared fixture. This does not reject all swirl models: coherent near-field constitutive response, nonperturbative states, or a different derived matter/field coupling are outside this calculation. Each would need a new Hamiltonian, stability and source-budget audit. Slower massive waves require rederiving their source coefficient and dispersion, not simply replacing c in this formula.

The conceptual advance is precise: an outward wave population can produce an extended 1/r acceleration, but the present Hamiltonian makes that route far too expensive for the declared targets. A useful next candidate must avoid depending on continuous relativistic energy escape and demonstrate an affordable stored or constitutive response. No such candidate is validated here. A finite front is essential; an infinite stationary 1/r^2 wave-energy distribution would have infinite total energy. Source depletion can invalidate the stationary approximation on the crossing timescale; that further weakens, rather than rescues, this route.

## Attribution and reproduction

Spherical flux conservation, Poisson integration and mass-energy accounting are established mathematics. The coefficient 1+s is derived from SR-1, whose shared metric optics is attributed in its protocol. Constants: [NIST G](https://www.nist.gov/how-do-you-measure-it/how-do-you-measure-strength-gravity) and [BIPM c](https://www.bipm.org/en/si-base-units/metre). Fixture masses/speeds are hypothetical. No cosmological distances, expansion or dark-matter component enter.

Run radiation_budget.py in a clean checkout without radiation-budget-v1 to reproduce; preserve the first archive. The manifest pins source and protocol hashes and source commit. Numerical identity checks do not constitute empirical validation.
