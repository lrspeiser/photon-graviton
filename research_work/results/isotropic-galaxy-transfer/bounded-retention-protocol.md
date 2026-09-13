# Two bounded retention revisions

Declared before fitting: change the physical role of the environmental response. Keep the capture opacity kappa=k0/[1+(r/a)^2]^2 and incoming attenuation unchanged; multiply only the retained source density by a probability eta(X)=X^q/(1+X^q). Fit C_half so rho_deposit=2*C_half*eta*J/[1+(r/a)^2]^2. Thus the source before retention is C_source=2*C_half, and eta is always between zero and one. q=0 reproduces the original profile exactly after normalization.

Two separate candidates:

1. Radiation: X=(L_3.6/10^9 L_sun)/(R_disk/kpc)^2.
2. Ordinary-well proxy: X=v_baryon^2(R_disk)/(1000 (km/s)^2), using the inherited interpolation/outer continuation. This is a circular-speed depth proxy, not the full potential or a dynamically inferred missing-mass field.

The pivot X=1 is fixed as specified, not adjusted after seeing predictions. q in [-2,2]. Fit C_half,k0,a/R_disk,q on 89 training galaxies; freeze for 29 validation and 31 test galaxies. Keep previous initializations/objective/input masses/distances. Record bounds, convergence, both log and km/s errors and per-object predictions. These partitions are already exposed.

This is a new capture-then-retention postulate using standard probability and transport mathematics. Fraction 1-eta must leave the retained sector, not disappear or be added again to deposited energy. The single-pass calculation assumes that released energy escapes the modeled system in a channel not recaptured there. A complete source equation must specify that channel and its outgoing flux. No global energy supply or microscopic release mechanism is derived. Within each fit C_source/k0 corresponds to the shared external exposure before retention; C_half/k0 alone is no longer that full exposure.

A radiation probability may represent promotion or disruption of storage depending on the fitted sign; the depth alternative tests stronger ordinary binding. Do not claim a mechanism merely from the sign. Compare to both the original interception model and the previous radiation-dependent opacity model. A change in the proxy's reference scale would change the saturation position and is not an innocuous unit change unless the pivot changes consistently.
