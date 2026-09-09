# Capturing energy also transfers momentum

## Why this is needed

Our receiver model can store energy, and the support calculation shows that its location matters. Now we ask whether filling it pushes it away or changes its orbit. Without this check, a model could balance energy while ignoring a force that prevents the intended galactic distribution from forming.

This pass assumes ordinary massless companions with momentum E/c, special-relativistic absorption and ordinary weak-field forces. A different companion dispersion law needs its own momentum derivation. No receiver population, capture strength, expanding universe, dark matter or final microscopic law is adopted.

Radiation-induced orbital forces have an established context in the [general relativistic Poynting–Robertson study by Bini, Jantzen and Stella](https://arxiv.org/abs/0808.1083). Here we independently derive a **mass-increasing absorption** comparison. We do not borrow constant-mass reradiation formulas or claim that their model establishes our deposit mechanism.

## 1. Exact local four-momentum accounting

In units c=1, a receiver of rest mass M and velocity v initially has four-momentum (gamma M, gamma M v). Deliver a radiation bundle of energy E and net spatial momentum xi E n, where n is a unit radial direction and -1<=xi<=1. Such a bundle can be formed from oppositely directed massless packets with appropriate energy weights. Full absorption gives

    P_final = P_initial + (E, xi E n)
    M_final^2 = M^2 + 2 E gamma M(1-xi n.v) + E^2(1-xi^2)
    v_final = (gamma M v + xi E n)/(gamma M+E).

Here xi describes momentum actually delivered to the receiver. Xi=0 does not automatically follow from an isotropic ambient radiation field: a moving absorber samples that field with Doppler shifts, aberration and a possibly frequency-dependent cross section. Those properties require an actual transport calculation.

The code checks 200 random cases, reconstructs each final four-momentum from its rest mass and velocity, and separates rest-energy gain from kinetic-energy change. If slowing the receiver releases some initial kinetic energy into its increased rest energy, the rest-energy gain can exceed the incident energy; the kinetic-energy decrease balances it. There is no extra energy source.

For a receiver initially moving tangentially at 200 km/s, an instantaneous bundle with energy equal to its initial rest energy illustrates the distinction. A balanced bundle roughly halves its tangential speed; a wholly radial bundle also gives a large radial speed. These large-bundle examples are exact local kinematics, not a slow galactic loading history or a prediction that such packets exist.

## 2. Continuous loading ties force to stored mass growth

Let P be delivered energy per unit local lab time. Differentiating the exact expressions gives the radiation contribution alone:

    dv/dt = P/(gamma M c^2) [xi c n-v]
    dM/dt = gamma P/c^2 [1-xi n.v/c].

The numerical finite-packet differences check these derivatives independently. In the slow-velocity limit, dM/dt is approximately P/c^2 and

    a_radiation approximately xi (P/Mc) n - (P/Mc^2) v.

The second term reflects increasing inertia. With radial delivered momentum, total receiver angular momentum about the center is conserved during absorption, but its angular momentum per unit rest mass decreases. One must not confuse that with a conserved orbit or with loss of total angular momentum through reradiation.

If some energy escapes through the shelving or leakage channels, both its energy **and momentum** must be subtracted from the delivered four-momentum. The earlier 81% storage example cannot simply be inserted into this full-absorption force formula without specifying the directions and frames of its outgoing radiation.

## 3. A quantitative directional-force requirement

In a specified galaxy potential, the gravitational acceleration at radius r is v_c^2/r. Define the fractional receiver mass-growth rate mu=(dM/dt)/M. For slow motion and complete local absorption, the ratio of radial radiation force to gravity is approximately

    a_radial/g = xi mu c r/v_c^2.

This is a local diagnostic; r, the potential and mu are held fixed for the comparison. For r=10 kpc, v_c=200 km/s and exponential receiver mass growth that doubles it in 10 billion years, mu=ln(2)/(10 billion years):

| Diagnostic | Value |
|---|---:|
| Outward force / gravity for xi=1 | 5.08 |
| Maximum absolute xi for radial force below gravity | 0.197 |
| Maximum absolute xi for radial force below one-tenth gravity | 0.0197 |
| Growth e-folding time giving equal force for xi=1 | 73.3 billion years |

These are chosen conditions, not observed mass growth or an adopted cosmic age. The force comparison does not evolve the trajectory or prove escape. Inward delivery would add inward force instead, so its spatial consequences also need calculation.

The amount of momentum per delivered energy is independent of how that energy is divided into ordinary massless packets. Smaller packets can reduce event-to-event noise but do not remove the mean force for a fixed total energy flow and direction.

## 4. Balanced delivery can still change an orbit

For exactly balanced delivered momentum, no external radiation torque, slow loading and an approximately circular orbit in a fixed external potential, total angular momentum is approximately M r v_c(r). Conserving it gives:

    Flat circular speed:       r proportional to 1/M
    Fixed Kepler potential:   r proportional to 1/M^2.

Doubling receiver mass would therefore halve its circular-orbit radius in the first comparison and reduce it to one-quarter in the second. These are adiabatic scaling laws, not integrated self-consistent galaxy predictions. Slow loading relative to the orbit and a fixed dominant potential are necessary assumptions. If the stored energy significantly changes the potential, source depletion and the evolving gravitational field must be solved together.

The point is not that every reservoir must shrink in precisely this way. It is that canceling the radial shove does not automatically preserve the receiver's phase-space distribution. A finite supported store needs a momentum source or redistribution process consistent with its energy ledger.

## Implication for the next unified calculation

Capture geometry, storage efficiency and orbital evolution must be computed together. Required inputs now include the actual angular/spectral companion field seen by the receiver, all emitted channels, the resulting rest-mass growth, finite storage capacity and angular-momentum exchange. A gravity-well preference cannot substitute for these equations.

The results constrain this full-absorption, ordinary-momentum comparison. They do not rule out balanced/coherent fields, momentum transfer to a larger supported system, or a different companion law; those alternatives need their own explicit receiving sectors and support. The complete conversion, capture, motion/lensing and 32-observation program remains unfinished.

Run `python -X utf8 research_work/results/capture-momentum/check_capture_momentum.py`. Saved results are in `capture-momentum-results.json`. These mathematical and local-force checks do not establish astronomical agreement.
