# PF5: five photon-written responses in a fixed nonexpanding universe

Declared 19 September 2026 before new computations. Baseline: 55a4b89.
User authorization: run all five proposed hypotheses thoroughly; no dark matter
or expanding universe in the tests; publish checkpoints to main regularly.

## Scope and common rules

No dark-matter source, halo benchmark, expanding distance calculation, scale
factor, or cosmological redshifting is used by this campaign. Geometry is
Euclidean and static. Previously adopted catalogue distances are retained as
numbers with their provenance. Previously exposed measurements are not called
fresh holdouts. This campaign does not run the old joint suite whose active
comparisons include excluded premises; its known RUT-8 failure stays recorded.

All added energy states start at zero. Photon donation, stored field energy,
heat, escape and boundary flux must be counted where that dynamics is modeled.
Conditional response coefficients are not demonstrated microphysics. A
prescribed radiation map is not proof of a self-consistent photon interaction.
Numerical verification, mechanism screening, and observational validation have
separate statuses. No parameter is selected on a declared transfer target.

### Data sufficiency

Use archived ordinary stellar photometry, X-COP gas density/thermal-pressure
profiles, six Coma shear bins and the six lens galaxies' light and stellar
motions. Do not read fitted halo masses as sources. Bolometric radiation maps,
stellar population histories, and measured merger entropy-production histories
are incomplete. Explicit morphology/emissivity proxies and controlled merger
simulations are therefore allowed, labeled as such. Their success cannot pass
an observational or funding gate that requires those missing inputs.

## Experiment 1: angular activation

P_ij = sum luminosity_weight n_i n_j, U = trace P, D = 27 det(P)/U^3;
D=0 when U=0. Writing is proportional to U D, not to a class label.
Test isotropic, beam, planar and rotated controls; source subdivision and
luminosity rescaling; disk thickness 0.05, 0.1, 0.2 of disk scale; spherical
cluster profiles. Use deterministic source sampling at 4096/8192/16384 points
and fixed seeds. Evaluate at center and 0.5, 1 and 2 characteristic radii.
SPARC surface brightness and X-COP j proportional to n_e^2 sqrt(T) supply
independent *shape proxies*, not complete bolometric maps.

Numerical gates: exact angular controls, rotation and subdivision <1e-10;
0<=D<=1 within roundoff; median morphology-score change <0.05 between the two
highest resolutions. Mechanism criterion: median cluster/disk activation at
matched scaled radius >=3, robust to all three disk thicknesses, with overlap
reported. No pretense that this alone fixes the mass-speed relation.

## Experiment 2: maturing transport

df/dt + v n.grad f + nu df/ds =
q delta(s)/(4pi) - lambda0 exp(-s)(f-fbar) - f/tau.
A positive particle-energy representation solves the kinetic equation.
Each packet is supplied by photons; directions isotropize at the stated rate.
Time unit = 1 Mpc/c, distance unit = Mpc; v=1 therefore equals c and never
exceeds it. All scale factors are constant. Primary: lambda0=3, nu=0.3,
tau=3, source duration=5 tau. Controls: nu=0, no scattering, source off.
Sensitivity: lambda0 in {0.3,3,30}, nu in {0,0.3,1}, tau in {1,3,10}.
Use 32768/65536 packets with repeated fixed seeds and disclose sampling errors.

The Green response is convolved with a declared Coma ordinary-emissivity proxy
(gas-density squared primary; stellar-shaped and mixed sensitivities). The
optical potential responds to transported energy density, not to an inserted
gravitating packet mass. One response normalization is matched at 1 Mpc to
the old pressure-calibrated photon-written response, before Coma shear is
scored. It is an exposed pressure anchor, not an independently funded coupling.
The Coma source-distance normalization remains one nonnegative nuisance.
Primary parameters are never reselected from the sensitivity scan.

Gates: packet normalization and analytic stored+decayed budget <1e-10;
causal displacement <= v age; isotropic symmetry within five standard errors;
doubling packets changes shear shape by <15% on its normalized maximum scale.
Mechanism criterion: primary Coma shape chi2 <= baryons-only chi2, and the
effect at 5 Mpc exceeds the original finite-footprint effect. No absolute
cluster success without source geometry and luminosity normalization.

## Experiment 3: fluctuation-amplitude response

tau_a^2 a_tt + 2 zeta tau_a a_t - ell_a^2 lap a + a = xi_gamma;
Phi_extra = -Vstar^2 sqrt(time_average(a^2)).
Use a fixed finite spatial covariance for the forcing; it is proportional to
declared photon emissivity, not rotation-speed measurements. Test the
inhomogeneous linear stochastic field by two independent stationary covariance
solves (Lyapunov and normal-mode identities), finite-time samples, luminosity
scaling, source subdivision, and zero source. Exact stationary energy injection
and damping are computed from that covariance.

Run planar test orbits in a radial response driven from an initially empty
field, with response averaging times {0.1,1,10} reference orbital periods,
three seeds, 20 measured periods after a declared burn-in, and mean extra
support 10% at the launch radius. Refine timestep and spatial grid. These are
negligible-mass probes, not a self-gravitating galaxy.

Gates: covariance agreement, energy identity and source subdivision <1e-9;
amplitude scaling exponent 0.5 +/-0.01. Orbit quietness criterion: <5% radial
spread and <5% median-radius drift; angular-momentum error <1e-5. Report noise,
drift and sensitivity even on failure. The stochastic-field ledger is not the
combined matter-field action: finite-mass backreaction and real-galaxy profile
agreement remain distinct required gates, not assumed successes.

## Experiment 4: shock erasure

u_t + div J = q_gamma - (1/tau + eta sdot_irr/cv) u.
Erased energy becomes gas heat. Evolve a one-dimensional head-on collision of
ordinary ideal gas using conservative Euler fluxes, with two luminous stellar
populations passing through it. This is a controlled merger, not a measured
cluster reconstruction or a complete three-dimensional gravity simulation.
Photon emission is removed from a finite luminous fuel account. Field starts
empty. Fixed source tracks are stated external kinematic inputs, not emergent
galaxy trajectories; no invisible mass is included.

Use gamma_gas=5/3, CFL<=0.3, grids 256/512/1024, initially separated gas lumps
and stellar lights at x=+/-2 moving toward each other at speed 0.7; erasure
eta={0,1,5}. Evolve through and after collision. Irreversible entropy
production is estimated from the conservative gas update against advected
entropy. Report its numerical/diffusive contribution. Compute optical
curvature from the resulting field, not merely the location of stored energy.

Gates: gas+field+fuel+escaped/boundary energy <1e-8 relative; positivity without
unaccounted floors; source-off field stays zero. Refined optical-peak locations
agree within 0.15 distance units. Mechanism criterion: erasure moves optical
peaks closer to the luminous populations than the matched eta=0 run, while
the gas stays near the collision center. Otherwise record failure. No observed
merger or freely moving self-consistent cluster is claimed.

## Experiment 5: directional optical strain

Q is symmetric and traceless, with positive kinetic/gradient energy and a
positive quartic potential. Photon Hamiltonian is
H=c sqrt(p^T exp(4 Phi_dyn I/c^2 - 2Q) p).
The homogeneous energy test uses the reciprocal photon force -dH/dQ, rather
than a permanently prescribed pump, and counts damping as heat. Starts are
Q=Qdot=0; include anisotropic and isotropic photon populations and zero photons.
Check positive optical metric, finite-energy response, zero coupling, and
Hamiltonian ray integration under timestep/tolerance refinement.

For the observational screening, source radial STF strain from the anisotropy
of the ordinary stellar light, solve its linear weak-field spatial response,
and compute its line-of-sight tensor deflection. Use only the registry's
already adopted static-Euclidean distances and masses. Source profiles and
stellar-motion calculation contain ordinary stars, with stated population
mass and constant-anisotropy freedom. No halo or expanding geometry is active.
Shared lengths {0.3,1,3,10,30} kpc and one shared positive coupling; calibrate
on the first five named lenses, carry to J1630 without fitting it. This is a
role-respecting split of exposed data, not a fresh holdout. Record kinematics
independently: optical fitting cannot repair a failed dynamical baseline.

Gates: homogeneous field+photon+heat ledger <1e-8; isotropic and source-off
responses <1e-10; tensor positivity; refined ray deflection <1e-4 relative.
Screening criterion: every lens bend within 3%, motion chi2 per measured bin
<=3, and weak tensor magnitude <1e-3. Report unfulfilled global funding or
reciprocal spatial-backreaction requirements explicitly.

## Checkpoints and outcome policy

Commit and push this declaration to main before execution. Publish verified
implementation/results checkpoints after experiments 1-2, after 3-4, and after
5 plus consolidated review. Fetch before each push; never force-push.
Keep failed numerical attempts, disclose amendments before replacement runs,
and do not silently change gates. New result directories are immutable;
reruns use new names. Final report includes every experiment, controls,
resolution/sampling tests, parameter sensitivities, energy scope and unresolved
observational inputs. Passing a numerical gate is not proof of a new gravity law.
