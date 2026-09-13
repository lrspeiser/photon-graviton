# A more concentrated capture profile

Before fitting, test one fixed structural change to the retained capture model:

    kappa(r) = k0/[1+(r/a)^2]^3, a=s*R_disk
    rho_deposit(r) = C*J(r)/[1+(r/a)^2]^3
    C = (k0/c) integral u_external(t) dt

The previous exponent is 2. Fix the new exponent at 3; do not fit it or add per-galaxy parameters. The physical hypothesis is a faster decline of capture susceptibility outside the galaxy, reducing the extended deposit tail. This profile is a postulate, not a microscopic derivation. Its radial tail is r^-6 instead of r^-4 before attenuation. Both source and attenuation must use the same exponent to preserve the common incident exposure relation.

Known radiative-transfer mathematics gives J=(1/2) integral exp(-tau) dmu. Set y=r/a, t=y*mu and B^2=1+y^2*(1-mu^2). The upstream line integral for exponent 2 is

    I2 = t/[2 B^2 (B^2+t^2)] + [atan(t/B)+pi/2]/(2 B^3).

The standard integration recurrence gives the new exact integral

    I3 = t/[4 B^2 (B^2+t^2)^2] + 3 I2/(4 B^2),
    tau = k0*a*I3.

Compute enclosed deposited mass by spherical integration and predict v^2=v_baryon^2+G*M_deposit(<r)/r, using known Newtonian mathematics. Keep ordinary matter inputs, distances, stellar mass-to-light assumptions, objective, bounds and training partitions unchanged. Fit shared C,k0,s on 89 training galaxies, then freeze for 29 validation and 31 test galaxies. Also fit the transparent control with the same new shape. All partitions have previously been examined, so this is exploratory transfer, not blind validation.

Use the inherited equal-galaxy logarithmic speed loss and report both logarithmic and km/s prediction errors. Adopt no change solely because training improves. Preserve old outputs. Record numerical drift and optimizer status. Compare frozen rotation results before proposing additional lens or Milky Way transfers. This profile alone supplies neither a retention mechanism nor proof of sufficient stellar energy supply.
