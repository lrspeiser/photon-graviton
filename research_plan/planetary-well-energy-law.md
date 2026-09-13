# One retained-energy to well-depth law across Solar System bodies

User request: translate solar energy into gravitational well depth and apply the same law to Earth, Moon and other bodies. This calculation separates a consistent known gravity law from the unproved hypothesis that incoming solar energy explains the existing wells.

## Universal conditional law

Define positive well depth W=-Phi with Phi(infinity)=0 for an isolated spherical body. W is energy per unit test mass (J/kg), not watts. For deposits confined inside the evaluation radius R, with negligible stress/binding corrections, the **known weak-field energy/mass and Newtonian relations** give

    W(R,t)=G/R * [M_ordinary + E_retained(t)/c^2],
    g(R,t)=W/R,
    v_escape=sqrt(2*W),
    dot W=G*P_retained/(c^2*R)     (fixed R and ordinary mass).

Thus a retained joule contributes G/(c^2*R) J/kg of depth. Its effect depends on where it is stored and where the potential is measured. There is no universal conversion from joules to a number of newtons without geometry and a test mass. Lifting a mass m at radial speed v against gravity requires mechanical power m*g*v; that power is not the same quantity as a well depth or deposited luminosity.

For nonspherical or extended weak-field sources the known generalization is

    Phi(x)=-G integral [rho_ordinary(x')+u_retained(x')/c^2]/|x-x'| d^3x'.

Ordinary matter already gravitates. Energy already included in a measured mass must not be added to it again. Traveling radiation, pressure and significant binding require the complete stress-energy description rather than this confined-energy approximation.

## Solar feeding specialization

Let a be heliocentric distance and R the receiving body's effective geometric radius. **Known inverse-square interception geometry:**

    P_intercepted = L_sun*R^2/(4*a^2).

Define eta as the fraction of this incoming original solar power that becomes permanently retained deposit energy. The **conditional feeding postulate and derived growth law** are

    P_retained=eta*L_sun*R^2/(4*a^2),
    Delta W=eta*G*L_sun*R*T/(4*c^2*a^2),
    Delta g=eta*G*L_sun*T/(4*c^2*a^2).

These assume constant distance, radius, luminosity and retention over T. Use time integrals otherwise. The eta=1 case is a hypothetical all-energy ceiling; actual absorption followed by thermal re-emission does not build such a reservoir indefinitely.

To apply our current **postulated gradual photon conversion** rather than full incoming energy, assume companions remain in the outward solar beam, no intervening capture, and geometric interception by the body:

    eta = eta_capture * {1-exp[-alpha*(a-R_sun)]},
    alpha=0.0002488993286382367/Mpc, 0<=eta_capture<=1.

Capture over a larger region, focusing or nonradial transport would require a derived effective area replacing pi*R^2. The table does not silently assume such enhancement. External non-solar companions are not included.

## Same calculation for every body

Masses, diameters and heliocentric distances use the rounded [NASA planetary fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/). The Moon uses Earth's approximate heliocentric distance, not the geocentric 0.384-million-km entry. R is half the listed diameter. These spherical estimates omit rotation, oblateness, other bodies and orbital variability; giant-planet effective surface gravity need not equal GM/R^2. Catalog masses derive partly from dynamics, so reproducing their gravitational potential is a consistency calculation, not an independent discovery.

| Body | Existing self-well depth W (J/kg) | Added W per year if ALL intercepted solar energy retained (J/kg) | Added W per year at fitted conversion, perfect companion capture (J/kg) |
|---|---:|---:|---:|
| Mercury | 9.029e6 | 1.632e-9 | 7.531e-25 |
| Venus | 5.371e7 | 1.159e-9 | 1.005e-24 |
| Earth | 6.247e7 | 6.391e-10 | 7.677e-25 |
| Moon | 2.804e6 | 1.741e-10 | 2.091e-25 |
| Mars | 1.262e7 | 1.465e-10 | 2.686e-25 |
| Jupiter | 1.772e9 | 2.646e-10 | 1.660e-24 |
| Saturn | 6.290e8 | 6.591e-11 | 7.610e-25 |
| Uranus | 2.267e8 | 6.974e-12 | 1.612e-25 |
| Neptune | 2.749e8 | 2.724e-12 | 9.921e-26 |
| Pluto | 7.304e5 | 7.638e-14 | 3.638e-27 |

The last two columns are different hypotheses, not additive contributions. At Earth's surface the current self-well depth means about 62.5 million joules per kilogram to escape the isolated, nonrotating Earth. The lunar figure is about 2.80 million J/kg. The same law yields spherical accelerations 9.795 and 1.614 m/s^2 respectively. These are distinct from the Sun's contribution to the potential at both locations.

## A direct consistency test of the solar-origin hypothesis

At essentially the same solar distance, equal eta and duration predict

    Delta W_Earth/Delta W_Moon = R_Earth/R_Moon = 3.67079,
    Delta g_Earth/Delta g_Moon = 1.

The existing mass/radius values instead give W_Earth/W_Moon=22.27880 and g_Earth/g_Moon=6.06921. Thus a hypothesis that builds their entire gravity from the same solar illumination history and common geometric capture efficiency fails this elementary ratio comparison. Increasing a shared duration cannot change those ratios. An Earth efficiency-times-duration 6.069 times the lunar value could algebraically fit them, but would require an independently specified physical reason, not a separate number selected for each body.

If solar deposits are only a correction to the existing ordinary mass, there is no such requirement to reproduce the entire well from illumination. In that interpretation the derived increments are small at the retained conversion rate. A new coupling that greatly amplifies the deposited energy's gravitational response is a separate hypothesis and needs a common field law and tests; it is not inferred from electromagnetism being stronger than gravity.

## Conservation and applicability to the Sun

The Sun obeys the same potential rule, but does not intercept its own isotropic luminosity through L*R^2/(4*a^2) at a=0. Its retained self-generated companions must use the [radial capture calculation](solar-system-photon-feeding.md). Solar radiated energy reduces the Sun's mass-energy; retained planetary or circumsolar deposits redistribute that supply rather than create it. A body's known existing mass includes its historical retained energy, so model it as M_ordinary+E_retained/c^2 or use its observed total mass, not both simultaneously.

The well-depth relation is consistent across all objects under its spherical weak-field assumptions. The tested solar feeding law is a conditional addition, not a demonstrated common cause of their existing gravity. No cosmic age or size is imposed and no per-body fit is made.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/planet-well-energy.py`. All input values, powers, mass equivalents and ratios are saved in `planet-well-energy-results.json`. The internal relation v_escape^2=2W is verified; no planetary ephemeris or anomaly validation is claimed.
