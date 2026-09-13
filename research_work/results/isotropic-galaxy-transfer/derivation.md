# From incoming companion rays to a deposited profile

Postulate a fixed spherical absorption coefficient, tied to the observed disk scale:

    a=s R_disk; kappa(r)=k0/[1+(r/a)^2]^2.

This is a phenomenological capture rule, not a microscopic graviton interaction or a new mathematical form. It declines continuously rather than ending at an assumed physical edge. Incoming companions propagate at c and retain their energy until absorption. Treat the external isotropic bath and capture coefficient as stationary during the integrated exposure. Straight rays neglect gravitational focusing.

Known linear transport gives I(r,mu)=I_external exp[-tau(r,mu)]. With x=r/a, t=x mu and B^2=1+x^2(1-mu^2), the exact upstream optical depth is

    tau=k0 a {t/[2 B^2(B^2+t^2)] + [atan(t/B)+pi/2]/(2 B^3)}.

This comes from integrating (B^2+t^2)^-2 along the incoming ray from minus infinity. The sign convention for mu can be reversed without changing the angular average. At the center tau=pi k0 a/4 in every direction. The mean surviving intensity is J(r)=1/2 integral_-1^1 exp(-tau) dmu.

Known energy accounting gives local absorbed power per volume q=4pi kappa I_external J=c kappa u_external J. If all absorbed energy is retained as ordinary gravitating mass equivalent, then

    rho_d=(1/c^2) integral q dt
         =C J(r)/[1+(r/a)^2]^2,
    C=(k0/c) integral u_external dt.

C is fitted exposure, not a claim that the stellar energy supply has been calculated. No maximum cosmic age or size is imposed. An incoming ray's energy is attenuated once along its path; the same incident energy is not deposited again at every interior point. A bath from infinity is a stationary boundary idealization, not a derivation of a finite-duration causal supply.

The spherical deposited mass and motion follow known Newtonian formulas M_d=4pi integral_0^r rho_d r'^2 dr' and v_pred^2=v_b^2+GM_d/r. The transparent approximation has J=1 and analytic mass

    M_d=2pi C a^3 [atan(x)-x/(1+x^2)],
    M_d(infinity)=pi^2 C a^3.

The attenuated density is nonnegative and no larger than the transparent profile, so its total mass is also finite for finite C,a. None of this supplies orbital support: retention is still assumed. Since opacity is prescribed by stellar scale rather than recalculated from total depth, this model also does not complete capture feedback. These limitations must accompany any observational fit.
