# Energy requirement and outer mass of the fitted capture rule

The previous turn tested the same source against stellar motion and lensing. This calculation connects that source's fitted amplitude to incoming energy and integrates its full predicted mass. No parameter is fitted again and no cosmic age or size is imposed.

## Exact capture area

For the postulated kappa(r)=k0/[1+(r/a)^2]^2, a complete ray at impact parameter b has the known line integral

    tau_chord(b)=T/[1+(b/a)^2]^(3/2), T=pi k0 a/2.

Its probability of being captured is 1-exp(-tau_chord). Integrating across the impact plane gives an effective absorption cross section:

    sigma_eff=2pi integral_0^infinity b[1-exp(-tau_chord(b))] db
             =pi a^2 {T^(2/3) gamma_lower(1/3,T)-1+exp(-T)}.

gamma_lower denotes the unregularized lower incomplete gamma function. This is a derived integral of known transport mathematics for the project's stipulated opacity, not a new physical interaction or an originality claim. An independent direct chord integral agrees with the expression to 5.4e-15 relative over the 149 actual galaxy scales.

The thin limit is sigma_eff approximately pi^2 k0 a^3. At large T, sigma_eff approximately pi Gamma(1/3) a^2 T^(2/3), hence proportional to k0^(2/3) a^(8/3). This differs from a fixed-radius perfect absorber's area law because this opacity has a declining infinite tail: increasing a changes both area and optical depth, and moves the effective absorbing region outward. It is not evidence that energy or geometric area is created.

## Total deposits and conservation

Known incoming-energy accounting for the stationary isotropic bath gives

    M_d,total=(sigma_eff/c) integral u_external(t) dt
             =(C/k0) sigma_eff,
    integral c u_external(t)dt=(C/k0)c^2.

The last expression includes consistent conversions for mass per area. It is an all-direction energy-flux convention, not the one-sided flux through a plane; for an isotropic bath that one-sided flux is c*u/4. Do not multiply a full sphere's surface by c*u and double-count directions. Companions captured by one galaxy are not available for capture by a second galaxy unless explicitly re-emitted.

The fitted C=2.71737e7 Msun/kpc^3 and k0=0.122551 kpc^-1 require

    C/k0=221.734 Msun/pc^2,
    integral c u_external dt=4.16189e16 J/m^2.

This is a model requirement, not a measured radiation field. For a constant assumed companion energy density u, the corresponding duration scales as

    duration=4.39912e14 years * (1e-14 J/m^3 / u).

The reference density is only a units/scaling example. It is not an adopted cosmic background, a measured companion density, an imposed universe age, or an exclusion of the hypothesis. A brighter or longer-lived bath changes the supply; conservation still requires a causal source history and competition among receivers to pay for it. Extending age alone also requires modeling opacity, illumination and retention throughout that history.

## Consequence using the actual galaxy scales

With each of the 149 observed disk scales retained, the median predicted total deposited mass is 1.853e11 Msun. The median fraction outside the last measured rotation radius is 76.40%. These are predictions of the tail and fitted source model, not observed halo masses. The individual mass requirements and covered-radius fractions are recorded in energy-results.json.

Outer shells do not affect interior spherical circular speeds, so good interior rotation behavior alone cannot validate that large unmeasured mass fraction. Lensing, more distant tracers and physical source/support calculations must constrain it. This is particularly relevant because the six-system lensing transfer is not yet satisfactory. The integrated finite mass is mathematically well-defined; its astronomical existence is not thereby demonstrated.

## What is resolved and what remains

The fitted deposition amplitude now has an explicit energy meaning and a calculable absorption area. It can no longer be described as an unspecified conversion factor that evades conservation. But the incoming companion density and stellar history have not been calculated, so this is not proof of adequate or inadequate supply. Microphysical capture, support, shared time/brightness behavior and the remaining observational goals are unresolved. All six goals stay open.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/energy.py`. The calculation keeps all fitted parameters fixed and reports every galaxy rather than selecting favorable examples.
