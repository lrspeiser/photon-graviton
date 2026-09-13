# Solar photons feeding a companion deposit

This is a user-requested quantitative specialization of our hypothetical conversion/capture model. It uses nominal solar luminosity, an approximate solar photon count, the retained fitted alpha, and an unspecified capture efficiency. It is not an observed detection of companions or a derivation of graviton production. No Solar System fit has been done.

## One photon: energy versus force strength

Electromagnetic interactions being stronger than gravity does not mean that a photon gains or loses a universal multiplier of energy when converted. **Known energy relations:** E_gamma=h*nu_gamma and, for a hypothetical graviton quantum, E_g=h*nu_g. If a fraction f of a photon's energy transfers into N_g gravitons with one frequency,

    N_g*h*nu_g = f*h*nu_gamma,
    N_g = f*nu_gamma/nu_g.

This is conservation bookkeeping, not a microscopic conversion law. A spectrum of graviton frequencies requires summing their energies. There is no established count ratio without their frequency distribution and a production mechanism. For illustration only, converting the full energy of an average solar photon into 100-Hz quanta would provide about 3.25e12 such quanta. 100 Hz is an arbitrary example, not a solar prediction. Momentum, angular momentum and the assisting interaction must also be specified before this can be an allowed physical process.

**Known stored-energy relation under ordinary gravity:** a confined, weakly gravitating deposit with negligible additional stresses has mass equivalent E_deposit/c^2. Its complete field depends on total stress-energy and geometry; this is not a rule assigning a fixed rest mass to a traveling graviton. No enhancement of gravitational coupling is assumed. If we propose one, it needs a separately defined/tested field equation.

## Solar photon output

Use IAU nominal L_sun=3.828e26 W and effective temperature T=5772 K ([IAU 2015 B3](https://arxiv.org/abs/1510.07674)). These are conversion constants, not exact instantaneous measurements.

**Known Planck-spectrum result applied as an approximation:** mean photon energy is

    E_mean = [pi^4/(30*zeta(3))]*k_B*T = 2.70118*k_B*T
           = 2.15260e-19 J = 1.34355 eV.
    photon_rate = L_sun/E_mean = 1.77832e45 photons/s.

This estimates photons escaping the Sun; internal absorption/re-emission events are not counted as additional net energy. The Sun is not an exact blackbody. Planetary reflected sunlight must not be counted again as an independent power supply. This is a Sun-dominated calculation, not a complete inventory of all internally generated planetary emission or externally arriving radiation.

## Main feeding formula

**Our postulated fractional conversion law, with its empirically fitted coefficient:**

    f_conv(R) = 1-exp[-alpha*(R-R_sun)],
    alpha = 0.0002488993286382367/Mpc, R>=R_sun.

It assumes the same unscreened coefficient locally; this is an extrapolation of the large-distance fit, not a locally measured rate. Use consistent distance units. Define eta_cap(R) as the fraction of converted energy generated on this outward path that is retained inside R. Then the **conditional derived deposit formula** is

    dot E_deposit(<R) = eta_cap(R)*L_sun*f_conv(R),
    dot M_deposit(<R) = eta_cap(R)*L_sun*f_conv(R)/c^2,
    M_deposit(<R,T)-M_initial(<R) = eta_cap*L_sun*f_conv*T/c^2.

The final line assumes constant luminosity, efficiency and retained deposits over duration T. Otherwise integrate the time-dependent power. It sets no age of the universe or Sun. eta_cap lies between zero and one; none of the current solar observations fixes it here. Perfect capture gives an upper bound, not a selected fitted value.

| Outer radius, measured from Sun | Fraction converted | Converted power (W) | Maximum retained mass equivalent per year (kg) |
|---|---:|---:|---:|
| 1 AU | 1.20109e-15 | 4.59776e11 | 161.44 |
| 30 AU | 3.61953e-14 | 1.38556e13 | 4,865.05 |
| 100 AU | 1.20664e-13 | 4.61903e13 | 16,218.58 |
| 1000 AU | 1.20669e-12 | 4.61922e14 | 162,192.61 |

Multiply the final column by eta_cap for a specified retention efficiency. Conversion is gradual energy loss, not necessarily disappearance of a fraction of photons. At 100 AU its transferred energy equals the full energies of 2.14579e32 average solar photons per second; this is an energy-equivalent count, not a prediction of that many destroyed photons.

As a separate scale, converting and retaining all solar luminosity would supply L_sun/c^2=4.25922e9 kg/s, or 1.34411e17 kg/year. Our fitted local path conversion is vastly smaller. The full-conversion number is not compatible with unchanged observed sunlight and is not our proposed solar prediction.

## A capture rule that can replace the free efficiency

**Known linear transport framework applied to our hypothetical species:** let x=r-R_sun, and let L_gamma and L_c be outward photon and companion luminosities. For a nonnegative capture opacity kappa(x),

    dL_gamma/dx = -alpha*L_gamma,
    dL_c/dx = alpha*L_gamma-kappa*L_c,
    dP_deposit/dx = kappa*L_c,
    L_gamma(0)=L_sun, L_c(0)=P_deposit(0)=0.

These preserve L_gamma+L_c+P_deposit=L_sun. They assume radial outward companion transport, no return rays, no escape from deposited states and a stationary source. Under spherical symmetry, after sufficient time for illumination,

    dot rho_deposit(r) = kappa(r)*L_c(r)/(4*pi*r^2*c^2).

For constant kappa!=alpha the known solution is

    L_c(x) = alpha*L_sun*[exp(-alpha*x)-exp(-kappa*x)]/(kappa-alpha).

At kappa=alpha its limit is alpha*L_sun*x*exp(-alpha*x). P_deposit=L_sun-L_gamma-L_c determines eta_cap. kappa=0 gives zero deposit even though conversion occurs. No solar kappa is asserted; it must follow a shared capture law and be tested against planetary dynamics, not be chosen to manufacture agreement. External incoming companions require an additional boundary flux, not a multiplier on solar luminosity.

## How deposits change the well

**Known Newtonian spherical-field approximation for the retained source:**

    g(r) = G*[M_sun(t)+M_deposit(<r,t)]/r^2 + other ordinary sources,
    Delta g_deposit(r) = G*M_deposit(<r,t)/r^2.

The force is inward and follows the deposited distribution. Material deposited outside a planetary orbit does not increase that orbit's radial acceleration in exact spherical symmetry; it can change the potential's constant offset there. Thus the total retained energy alone does not predict every planet's motion.

Solar emission removes energy from the Sun: its photon-related mass-loss rate is -L_sun/c^2. Capturing some of that energy nearby keeps it within the system; it does not add the same energy on top of an unchanged Sun. In a long-lived steady outward flow, this photon budget gives dot M_sun+dot M_deposit=-(1-eta_cap*f_conv)*L_sun/c^2, apart from the traveling inventory and other solar channels. The source's original mass-energy, emitted photons and deposits must not be counted simultaneously as three independent supplies. The Sun's ordinary gravitational well also does not require continual photon feeding merely to exist.

## Reproduction and status

[Executed calculation](../research_work/results/isotropic-galaxy-transfer/solar-photon-budget.py) and [numerical output](../research_work/results/isotropic-galaxy-transfer/solar-photon-budget-results.json) give the inputs and stored values. Photon-energy reconstruction and conversion bounds are checked. This establishes a quantitative energy ledger and conditional Solar System scaling, not a gravitational anomaly, allowed microscopic decay, capture solution or complete theory.
