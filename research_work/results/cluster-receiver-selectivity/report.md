# Does loading ordinary matter equally favor the galaxy locations?

Under the narrow additive receiver model examined here, **no**: equal positive
loading of stars and gas preserves the ordering of the two published aperture
masses. Selective stellar loading can reverse it, but that is a requirement to
explain, not a prediction we have obtained from photon capture.

## Conditional model, not a fitted lens map

Suppose deposited energy stays with ordinary receivers and contributes ordinary
positive gravitating mass. Let B_star and B_gas be total effective mass divided
by original ordinary mass for stellar and gas receivers. Both are at least one.
Assume each factor is common to both apertures of the system. The effective
aperture mass is then

    M_eff = B_star M_star + B_gas M_gas.

This is an effective postulate using established additive mass accounting; it
is not a new gravity law. It does not cover an independently moving field,
nonlocal gravity response or spatially varying capture history.

Write B for the galaxy-centered aperture and P for the plasma-centered aperture.
When gas_P > gas_B and star_B > star_P, elementary algebra gives

    M_eff(B) > M_eff(P)
      iff B_star/B_gas > (gas_P-gas_B)/(star_B-star_P).

This is a threshold for **aperture ordering only**. Two aperture totals cannot
determine a lensing peak, and this inequality is not a necessary condition for
every lensing geometry or alternative gravity theory.

## Published inputs and numerical result

We use Table 2 of [Clowe et al. (2006)](https://arxiv.org/abs/astro-ph/0608407),
with its reference aperture and mass calibration. The original stellar conversion
assumes I-band mass-to-light ratio 2; the text discusses 0.5 to 3. Quoted stellar
errors omit that systematic. Gas masses also require an emission/deprojection
model. These masses are not raw, theory-free observables.

| Component | Gas B | Gas P | Stars B | Stars P | Required B_star/B_gas |
|---|---:|---:|---:|---:|---:|
| Main | 5.5 | 6.6 | 0.54 | 0.23 | >3.55 |
| Subcluster | 2.7 | 5.8 | 0.58 | 0.12 | >6.74 |

Mass entries are in 10^12 solar masses under the source calibration. The
published lensing entries use different profile subtractions at the two
positions; we do not treat their ratio as an unsubtracted mass ratio.

Changing a common stellar mass normalization, while keeping the tabulated gas
normalization, gives this limited sensitivity check:

| Stellar mass scale relative to table | Main threshold | Subcluster threshold |
|---|---:|---:|
| 0.25 | >14.19 | >26.96 |
| 1 | >3.55 | >6.74 |
| 1.5 | >2.37 | >4.49 |

These are conditional algebraic results, not confidence intervals. No
independence of aperture errors is assumed and no covariance is invented.

## Dependence on the proposed universe's distance and light law

Unlike the earlier purely angular benchmark, this calculation imports published
mass calibration. It is therefore NOT an expansion-independent quantitative
bound. If a revised geometry/brightness analysis rescales gas masses by q_g and
stellar masses by q_s, the threshold rescales by q_g/q_s, provided each scale is
common to both apertures. Changing deprojection or angular aperture selection
can require more than a uniform rescaling. We do not assume that the alternative
photon law preserves standard luminosity or X-ray mass conversion.

## Energy and momentum consequence

In the reference subcluster comparison, even setting B_gas=1 requires B_star>6.74.
If the extra mass is internal stored energy in the stellar receivers and ordinary
mass-energy equivalence applies, this means stored energy greater than 5.74 times
their original rest energy. That is not an inferred reservoir measurement or a
funded energy supply. It also increases inertia in this branch; it cannot be
inserted only into the gravity calculation while leaving receiver mechanics fixed.

The earlier [capture momentum derivation](../capture-momentum/derivation.md) and
[isotropic absorption calculation](../isotropic-capture/derivation.md) show how
delivery affects velocity and momentum. The present comparison does not evolve
that loading. Nor does it establish that real stars could store such energy
without changing their structure or emission.

Thus universal loading per unit ordinary mass does not resolve this particular
aperture contrast. A stellar-selective receiver would need independently
calculated cross sections, exposure histories, capacity and subsequent dynamics.
An independent bound field avoids being tied to stellar internal energy, but
needs its own transport and stress equations. Neither alternative is selected
by matching this threshold.

`run.py` verifies the threshold against the direct two-aperture mass calculation
above and below each crossing. That verifies algebra, not the model's physical
or observational validity. No reserved observations were used; cluster collision
and common-source lensing predictions remain unfinished.
