# Direction-dependent disk capture: bounded comparison

13 September 2026, specified before execution. Same 149 exposed SPARC galaxies,
3150 radii, baryonic contributions and original 89/29/31 split. No new catalog
or observed vertical motions. Preserve exact-third eta and reference constants.

Replace spherical opacity by kappa=k0 W(n_z)/(1+m squared/a squared)^2,
m squared=R squared+z squared/q squared. Incoming intensity is isotropic at the
outer boundary. Top and bottom remain symmetric. Test four specified shapes:
q=0.5 with W=1; q=0.5 with W=3 n_z squared; q=0.5 with
W=1.5(1-n_z squared); and q=0.25 with W=1. W has angular mean one and represents
capture-direction sensitivity, not an observed anisotropic radiation sky.
Attenuation uses that same W, so direction changes both absorption and screening.

For each, rho=2 C0 eta f(m) <W exp(-tau)> over incident directions. Compute an
axisymmetric density and its equatorial Newtonian force, not spherical enclosed
mass gravity. Use even Legendre multipoles and test numerical convergence.
Compare raw capture-derived normalization and normalization to the reference's
mass over the same finite domain. The latter is a redistribution diagnostic, not
a derived transport solution; no per-galaxy amplitude is fitted to observations.
Report tail bounds and raw inventory changes.

For each of eight branches fit one shared mixture f_mix in [0,1] between the
unchanged spherical reference and the directional endpoint, using training-only
equal-galaxy mean squared log10 speed. Scan 101 values and polish minima and
endpoint neighborhoods. Freeze f_mix for validation/test. Endpoints and s=0
identity control remain explicit. Reject nonpositive total circular speed squared
rather than clipping an outward deposited force to zero.

Report signed inner/middle/outer errors at R/Rdisk boundaries 1 and 3, split
speed and log RMS, existing MOND control, raw versus fixed inventories, parameters
and numerical limits. Refine geometry/force predictions and check frozen-parameter
scores before selecting any useful candidate. Multiple tried families and exposed
samples preclude blind significance claims. Geometry here is an idealized capture
region tied to Rdisk, not a measured full disk/bulge opacity or a physical creation,
support or lensing model. Do not promote a branch based only on average bias.
