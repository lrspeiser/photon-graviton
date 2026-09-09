# Environmental suppression without automatically blocking companion waves

## Purpose and assumptions

The cumulative-time candidate needs a physical field configuration, including what happens inside the galaxies containing sources and detectors. This calculation tests an optional restoring interaction that makes the companion/time field small in a matter-rich region while allowing it to evolve in a void. It asks whether doing so necessarily prevents companion energy from entering that region.

The companion sector remains the agreed recipient of photon energy. This is not ordinary gravitational time dilation. Using local matter density to control the restoring term is a new candidate assumption, not a derivation from low gravity: matter density and gravitational potential are not interchangeable. No dark-matter population or expanding geometry is introduced.

## 1. Positive restoring energy and a slow field profile

Write n=1+phi and consider the field part of an effective action

    L=K/2[(phi_t)^2-v^2*(phi_x)^2-m(rho)^2*phi^2],  K>0.

The new term gives positive potential energy K*m^2*phi^2/2 and favors phi=0. In this slab comparison, m=0 outside and m>0 inside a region of half-width a. A microscopic density-dependent matter action is not supplied; the mass profile is a prescribed environment. Momentum forces on that environment must be included in a dynamical completion.

With negligible photon driving inside the region and equal boundary values phi(-a,t)=phi(a,t)=A(t), the field equation is

    phi_tt-v^2*phi_xx+m^2*phi=0.

For constant or linearly changing A, the exact interior solution is

    phi(x,t)=A(t)*cosh(m*x/v)/cosh(m*a/v).

Thus the central field and its linear drift are suppressed by

    f=sech(B),  B=m*a/v.

Six independent initial-value ODE calculations verify this profile, including strong-suppression parameters. The profile requires compatible initial/boundary data; a suddenly imposed boundary value also launches transients. It is not a full self-consistent void-plus-galaxy solution including photon backreaction.

As a clock-response illustration only, the earlier q(n)=n^-2 scaling gives q_center=(1+f*A)^-2. When f is small, even an order-one exterior field produces a small interior change. For sources and detectors with the same f, the observed factor would be J*[(1+f*A_emitted)/(1+f*A_received)]^2 under that particular matter response. Actual clock constraints and density coupling are not derived here; shielding cannot be declared adequate merely because f was chosen small.

## 2. The same term has a frequency-dependent effect on waves

Perturbations outside the slab obey omega=v*k; inside they obey omega^2=v^2*q^2+m^2. Waves below m are evanescent inside the material, while waves above m propagate. Continuity of phi and phi_x across both interfaces gives a lossless transmission fraction T.

For w=omega/m and slab thickness 2a:

    w<1: T={1+sinh^2[2B sqrt(1-w^2)]/[4w^2(1-w^2)]}^-1,
    w=1: T=(1+B^2)^-1,
    w>1: T={1+sin^2[2B sqrt(w^2-1)]/[4w^2(w^2-1)]}^-1.

Reflection is R=1-T. The calculation independently solves the complex interface-matching equations for 21 moderate-strength/frequency cases, including the threshold limit. Both transmission formulas and flux conservation agree to numerical precision. Large-B illustrations use the analytic expressions; their tiny transmission values are not claimed as independent high-condition-number matrix measurements.

At w=2, the worst possible interference phase still gives T>=48/49, approximately 97.96%. At higher frequencies this lower envelope rises further. Consequently strong suppression of a slowly changing field does not imply opacity to all companion waves. Which regime actually matters depends on the companion spectrum, which must come from the photon-to-field interaction rather than be selected to ensure transmission.

## 3. Examples and scale

| Center suppression f | Required B | Transmission at omega=0.1m | Transmission at omega=m | Transmission at omega=2m |
| --- | ---: | ---: | ---: | ---: |
| 0.001 | 7.6009 | 1.15e-14 | 0.0170 | 0.9823 |
| 0.000001 | 14.5087 | 1.32e-26 | 0.00473 | 0.999999 |
| 1e-12 | 28.3242 | 1.75e-50 | 0.00124 | 0.9909 |

The values above threshold oscillate because of reflection at both interfaces; they are not monotonic in thickness. The general lower bound at w=2 is more robust than a favorable transmission resonance.

For an illustrative half-width a=10 kpc and v=0.5c, millionfold center suppression requires m about 7.05e-12 radians per second, equivalent to hbar*m about 4.64e-27 eV. This is a selected size and model scale, not a measured particle mass. The earlier empirical effective stretching rate is about 2.42e-18 per second, much slower than this illustrative restoring frequency. That comparison shows a possible scale separation; it does not fit a void coupling, derive companion frequencies or demonstrate observational compatibility.

## 4. Energy conservation and what this does not capture

The field energy density and flux are

    e=K/2[phi_t^2+v^2*phi_x^2+m^2*phi^2],
    F=-K*v^2*phi_t*phi_x.

For the exact linearly driven profile,

    dE_inside/dt=2K*v*m*A*A_dot*tanh(B)
                =incoming boundary energy flux.

Six volume-integral/boundary-flux checks verify this equality. A field growing inside the slab acquires energy from the surrounding field, not from an unrecorded source. In a full model the exterior field and its photon-driven energy balance must be evolved too. Prescribing a boundary amplitude is not permission to supply free energy.

The wave-scattering model has **zero absorption**. A transmitted wave has reached the region but is not captured; an evanescent or reflected wave is not a permanent deposit. A restoring potential alone cannot be credited with the required irreversible deposition or extra gravity. A receiving bound state or other capture interaction, with its complete energy and momentum balance, is still needed.

## 5. Result and next work

This candidate shows that a slowly evolving time/companion field can be strongly reduced inside a region without blocking sufficiently high-frequency companion disturbances. The frequencies, environmental coupling and spatial scales must be jointly constrained rather than independently tuned.

Next derive a matter coupling that produces the restoring profile and its forces, compare the interaction-generated companion spectrum with its transmission cutoff, and then test a real capture channel. The field must still yield the required cumulative photon stretch, a stable source/detector response and compatible gravitational motion and lensing. The result is an optional mechanism comparison, not an adopted theory or a passed astronomical test.

```sh
python research_work/results/environmental-screening/check.py
```

Requires NumPy and SciPy. The standard offline diagnostic suite includes the calculation.
