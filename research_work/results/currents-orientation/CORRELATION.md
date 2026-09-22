# JR-6 addendum: the exact connection from overlap to mean conversion

A separate analytic exploration after the primary current toy. This is not a new evolved galaxy solution or an inferred astronomical parameter.

Let an overbar denote angular averaging at fixed radius, height and time. For net exchange Q=kp*P-km*C, the mean is exactly

    mean(Q) = mean(kp)*mean(P) - mean(km)*mean(C)
              + mean(delta_kp*delta_P) - mean(delta_km*delta_C).

The last two terms retain spatial correlation between material response and illumination/companion population. They vanish in the primary axisymmetric-rate toy. In clumpy gas their values can depend on orientation and relative passage. They are not extra freely fitted galaxy multipliers.

For an instantaneous illustration, take

    kp(phi)=k0[1+b*cos(2phi)]
    P(phi)=P0[1+a*cos(2(phi-delta))].

The angular means remain k0 and P0. Their mean forward transfer is

    mean(kp*P)=k0*P0[1+(a*b/2)*cos(2delta)].

For the illustrative a=b=0.8, aligned patterns give 1.32 times the unstructured forward-transfer rate, axes offset 45 degrees give 1.00, and perpendicular axes give 0.68. Rates and densities stay positive. The parameters were not fitted to observations. Independent 4096-angle quadrature at offsets 0,15,30,45,60,75,90 degrees agrees with this expression to maximum absolute error 2.23e-16; the means of rate and radiation remain one to roundoff.

These values are NOT 32 percent more or less gravity: reverse transfer and the subsequent spatial response must also be calculated. Energy transferred from P enters C; correlation does not create new total energy. Uniform sampling of a full relative rotation averages this illustrative contribution to zero unless sustained alignment, residence, nonuniform illumination history or nonlinear feedback prevents it.

This supplies a specific next physical link: couple the moving distributed source to nonaxisymmetric gas and its reverse-conversion pattern. Determine the correlations dynamically rather than choosing a separate gas-gravity coefficient for each galaxy. Density-only or angle-averaged source descriptions discard the terms that could distinguish enhancement from suppression.
