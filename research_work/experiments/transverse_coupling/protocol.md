# TF-1: transverse coupling of the converted field
Declared 19 September 2026 before any calculation. Baseline `main` 9aff6b2, developed on
`claude/tests-clusters-lensing-xrhttm`. The owner's CWC-1 constraints are retained: static
coordinates, no dark matter, no expansion, no distance refit, fixed adopted observations, every
energy counted, every failure preserved. Formula provenance follows
[research_plan/formula-provenance.md](../../../research_plan/formula-provenance.md); the
attribution register is [provenance.md](provenance.md).

## The question, in plain language (the owner's communication rule)
CWC-1 showed that light's energy can be converted into a receiving field and that the field bends
light toward where it was generated. It failed three requirements: timing (the field changes light's
local speed as measured by material clocks), colour (equal-energy carriers of different wavelength
are stretched differently) and universality (the force on matter scales with internal action over
mass, not with mass alone). All three come from one feature: the field acts *along* the motion. It
changes light's speed and energy, and it pulls matter down a well.

TF-1 asks whether the same field can act *sideways* instead: turn light and stars toward the matter
without changing their speed or energy. Two sideways laws are tested, a swirl law with the structure
of magnetism and a steering law whose push does not depend on speed. The stage measures (i) which of
the campaign's failed requirements a sideways action passes exactly, (ii) what it costs in energy
conversion, since a sideways push does no work and therefore converts nothing, (iii) whether it can
focus light and support pressure-held systems at all, and (iv) whether, on the six SLACS lenses with
resolved stellar motions, a sideways action changes CL-2's verdict.

Why it matters: the objective is a gravity that keeps matter together, bends light closer and keeps
stars in orbit without unseen mass. A sideways action is the only kind that avoids the timing and
colour costs, because it does no work on light. If it fails the lenses too, the obstruction is the
field's shape and the programme returns to the shape (CL-2 stage 2). If it passes, the missing
ingredient is identified and its field theory becomes the next task.

What would stay unreliable if this were skipped: every energy-to-gravity conclusion would rest on a
coupling known to fail timing and colour; the sideways alternative would remain an untested idea; and
the speed-of-light mismatch of the swirl law would remain an assertion rather than a measurement.

## The two transverse candidates and the longitudinal channel
Units and fixtures are CWC-1's (c_ref = 1, length 1). Let P_perp(n) = 1 - n n^T for a unit vector n.

**(V) vector channel.** Established structure (Lorentz force; the gravitoelectromagnetic analogy,
see provenance.md); its use as a converted receiving field is a candidate, not a derivation.

    matter:  dv/dt = v x B,           B = curl A,  per unit inertial mass
    light:   dn/dt = n x B,           |n| = 1, per unit energy (achromatic by construction)

**(S) steering channel.** Proposed here; originality unverified.

    light:   dn/dt = g_perp P_perp(n) grad phi,   |n| = 1, frequency unchanged
    matter:  dv/dt = g_perp P_perp(v/|v|) grad phi, per unit inertial mass

It is the transverse part of CWC-1's index bending with the longitudinal part removed. To first order
in g phi the deflection of a ray by the index a = exp(-g phi) is g times the line integral of the
transverse gradient of phi (established: eikonal/Born approximation), which (S) reproduces with
g_perp = g. (S) has no Lagrangian in the family tested by G6 and no reciprocal recoil term; the
momentum it removes from light is not deposited anywhere. That is a declared limitation and a
promotion gate, not a hidden assumption. It cannot act on matter at rest (P_perp is undefined at
v = 0; the force is taken as zero there) and it does not turn radially moving matter.

**Longitudinal channel.** CWC-1's index and clock couplings with coefficient g_par:
a = exp(-g_par phi), clock rate exp(-b g_par phi), force b g_par I Omega grad phi. Conversion of
light energy into the field happens only through this channel, because a force perpendicular to the
velocity does no work (established mechanics). The transverse channel is exactly silent in the 1D
fixture, where no transverse direction exists.

**A derived limitation, recorded before the runs (derived here; originality unverified).** For a ray
Hamiltonian H = |p| h(x, p/|p|) that is homogeneous of degree one in p (achromatic), the ray speed is
|dx/dt|^2 = h^2 + |grad_{p-hat} h|^2. Exactly unit speed for every direction therefore forces h = 1
(free rays) or a degenerate H = p . e(x) whose rays are slaved to a fixed direction field. So no
achromatic Hamiltonian ray law bends light while keeping its speed exactly. The steering law escapes
by not being Hamiltonian; the vector law keeps a first-order speed change along A.

## Fixtures
- Primary generated field for the ray gates: CWC-1's archived 2D field
  `coupled_conversion/evidence/spatial2d-v1/primary-final-field.npz` (256^2, produced at commit
  13da431), loaded by SHA-256 and never modified. The archived index-ray angles in that directory's
  results are the comparison values.
- The CWC-1 1D fixture (L = 48, right-going packets at -14 and -10, width 0.8, wavelength 2, receiver
  at x = 6, four material oscillators, end time 26, n = 1024, b = 2) rerun through CWC-1's own
  `spatial.simulate` for the longitudinal scan. Archived 1024 cases (vacuum, colour vacua, coupling
  0.05/0.2/0.5 at clock exponents 0/1/2) are reused where they exist.
- Model fields for the geometry and orbit gates, explicitly fixtures and not astrophysical sources:
  a solenoidal vector field B = curl(A_phi phi-hat), A_phi = A0 R exp(-(R^2+z^2)/2w^2), w = 1; and a
  spherical steering field with g_T(r) = g0 r/(r+1)^2 (Hernquist-like force shape) together with the
  well of the same g_T as control.
- Lenses: the six SLACS lenses through CR-2's measurement interface as built in CL-2
  (`results/path-memory/cl2_sources.LensSystem`): published light-profile components, KCWI apertures
  and covariance, stellar masses per IMF, constant anisotropy beta in [-2, 0.45]. Primary geometry is
  the static Euclidean scenario PF1 (D = ln(1+z)/alpha0), consistent with the owner's no-expansion
  constraint; flat FLRW is a recorded sensitivity only. No halo model is fitted. The held-out lens
  remains unavailable (CL-2 amendment 1).

## Stage 1: exact properties (numerical gates, each with a negative control)
G1 work-free. Rays and probe particles integrated in the archived field: under (S) the direction is
   integrated as an angle, so the speed is one by construction and the frequency is not a variable;
   under (V) the speed drift over the crossing is < 1e-9 (DOP853, rtol 1e-11). Control: CWC-1's index
   rays, whose local speed a differs from one by up to g max|phi| > 1e-3 in the same field; rejected.
G2 composition independence. Two probes at one position and velocity with internal action over mass
   differing by two: accelerations equal to 1e-14 under (S) and (V). Control: CWC-1's material force
   at the same position, ratio two; rejected.
G3 colour and timing of the transverse channel. Exact: the channel enters neither the dispersion
   relation nor the clock term, so spectral stretch = event stretch = 1 and the fixed-ruler speed
   change is zero identically. Recorded as an identity, with the longitudinal scan (G4) measuring the
   only nonzero costs.
G5 steering rays on the archived field, six impacts +-0.5, +-1, +-2 from x = -8 to +8, g_perp = 0.2:
   (i) all six inward; (ii) relative difference from the archived index-ray angles < 2% of the largest
   bend, applicable because g max|phi| = 0.004 < 0.02; (iii) ray-step refinement changes < 1e-4 of the
   largest bend; (iv) zero field gives exactly straight rays; (v) the mirror residual is recorded
   beside the archived 1.78% (source asymmetry), not gated.
G6 Lagrangian family. For L = |v|^2/2 + psi(x)|v|^alpha, the Euler-Lagrange equations give
   a_par = (1-alpha)|v|^alpha (grad psi . v-hat)/(1 + alpha(alpha-1) psi |v|^(alpha-2)) and
   a_perp = |v|^alpha P_perp grad psi/(1 + alpha psi |v|^(alpha-2)) (derived here). Finite-difference
   Euler-Lagrange checks at 20 random states for alpha in {0, 0.5, 1, 2} agree to 1e-6. Reading:
   the speed is conserved only at alpha = 1, where the transverse force scales with |v|; a
   speed-independent transverse force is not in the family (alpha = 0 is a potential force, the
   control).
G7 velocity averages. (V): the mean of v x B over an isotropic Maxwellian population is zero within
   three Monte-Carlo standard errors at N = 1e6, and equals v_rot B for a co-rotating population to
   1e-3. (S): the mean support factor f(beta) = 1 - <v_r^2/v^2> for an anisotropic Gaussian with
   sigma_t^2/sigma_r^2 = 1 - beta equals the one-dimensional quadrature
   f = 1 - integral_0^inf (1+2s)^(-3/2) (1+2s(1-beta))^(-1) ds (established: Laplace identity for
   Gaussian moments) to 2e-3 against Monte Carlo at beta in {-2, -1, 0, 0.45}; f(0) = 2/3 to 1e-10.
G8 focusing geometry. Rays through the solenoidal model field from three viewing directions (along
   the axis, in the plane, at 45 degrees), impacts +-0.5, +-1, +-2: the focusing fraction, the
   deflection-weighted mean of the cosine between each deflection and the direction from the ray to
   the axis, is > 0.999 for (S) and for the well of the same field magnitude, and |fraction| < 0.1
   for (V). Recorded beside it: the light-to-matter transverse acceleration ratio at equal field for
   (V), c/v, at v = 200 and 300 km/s.
G9 orbits in the spherical steering field. A circular orbit keeps its radius to 1e-8 over ten periods
   and its period equals the well's to 1e-8; a radial orbit crosses unturned (direction and speed
   unchanged to 1e-12); an isotropic ensemble of 2000 probes released inside r = 2 is followed for 20
   dynamical times under (S) alone and under the well alone (control), recording the bound fraction
   and the anisotropy drift. Energy or speed conservation < 1e-8 per orbit is the numerical gate; the
   bound fractions are a reading, not a gate.

## Stage 2: the longitudinal scan (how much conversion the allowed costs buy)
Rerun the 1D fixture at n = 1024, b = 2, for g_par in {0.005, 0.01, 0.02, 0.1} at wavelengths 1, 2
and 4 at equal initial EM energy, plus the g_par = 0.2 primary as a reproduction run; reuse the
archived vacua and the archived g_par = 0.05 and 0.5 cases. Metrics as CWC-1: receiving gain,
spectral stretch, event stretch, their discrepancy, the equal-energy colour spread, the fixed-ruler
speed change.
Numerical gates: relative energy error < 1e-3 per run (CWC-1's); the reproduction run matches the
archived primary-1024 receiving gain and pulse metrics to 1e-6 relative; the fitted power laws
(gain and each cost against g_par over the scan) have residual scatter < 10% in the log.
Reading: the largest g_par at which timing (1%), colour (1%) and fixed-ruler speed (0.1%) all pass,
and the conversion it delivers. The comparison with the cosmological rate is arithmetic on the
observed rate H0/c (about 2.3e-4 per Mpc) and is labelled as such; no fixture-to-galaxy mapping is
made.

## Stage 3: the six lenses (exploratory, exposed data)
Two force families on each lens's deprojected stars, all amplitudes nonnegative:
- F1, the universal family of CL-2, widths >= 0.15 kpc for every lens (21 widths, 0.158 kpc to
  10 Mpc), one amplitude vector shared by all six lenses (the joint test) and, separately, free per
  lens (feasibility);
- F2, a free profile per lens: 16 thin shells at radii log-spaced from 0.2 to 50 kpc (a
  nonnegative enclosed-mass profile of no assumed form), the ceiling of what any shape can do.
Two couplings on each family:
- W, the well: light bends with the factor two (the L1 rule, CL-1), stars feel the full force;
- S, the steering law: light bends with factor one, stars feel f(beta) times the force through the
  first velocity moment of the collisionless Boltzmann equation (mean field; a steady state with this
  force is not established, which G9 measures).
(V) is not fitted: G7 gives zero mean support for a non-rotating population and G8 no focusing; the
c/v ratio is recorded as the reason.
Per lens and IMF: the Einstein row is imposed as a constraint (weighted 1e3 above its 3% error) and
the kinematics chi-square (in V, release covariance) minimised over beta by the CL-2 grid-and-bounded
search; the joint F1 solve alternates as CL-2 did, with iteration counts and objective monotonicity
recorded this time. Each nonnegative solve reports its optimality residual (< 1e-8).
Reproduction gates: the stars-only Einstein angles under PF1 reproduce CL-1's archived
theta_newtonian_L0 (Chabrier) for all six lenses to 1e-6; the stars-only chi-square and beta under
flat FLRW reproduce CL-2's archived CR-2 benchmarks to 1e-6.
Acceptance as CL-2: the lenses are described if max |Einstein residual| <= 3% and the total
kinematics chi-square <= 128 over 40 bins.

## Decision rule (fixed now)
(a) S passes G1-G9 and, under the joint universal family F1, describes the lenses where W does not:
    the transverse coupling is the missing ingredient; next, its reciprocal field theory.
(b) S passes G1-G9 but leaves the lens verdict unchanged under F1: the transverse coupling removes
    the timing, colour and universality obstructions of the scalar channel but not the shape
    obstruction. The F2 ceilings then read the light-to-star ratio: if the per-lens ceiling under S
    is worse than under W by more than a factor two in total chi-square, the ratio of S is
    disfavoured; better by a factor two, favoured; otherwise indistinguishable.
(c) S fails an exact gate or bends light outward: rejected.
(V) is rejected if G7 and G8 come out as declared; a pass of either keeps it.
Whatever the outcome, the three statuses stay separate: reproduction, numerical verification and
scientific reading. Outcome (a) is a description of exposed data, not a validated theory.

## Promotion gates (what would be needed before calling S or V a field interaction)
A reciprocal formulation: a field equation that receives the momentum the coupling takes from light
and matter, with the energy budget closed; the transverse channel's sourcing (V needs energy
currents with curl; S needs a retained scalar profile), formation and retention at the required
strengths; three-dimensional rod and clock completion; the galaxy and cluster tests with shared
coefficients. None is attempted here.

## Evidence, reporting and checkpoints
Commit this protocol and provenance before any calculation. Each stage runs into a new immutable
directory `evidence/<stage>-v1` with a manifest (commit, source hashes, consumed-input hashes,
versions, elapsed time) written before and after the run. Failed gates are recorded, never
re-thresholded; a corrective rerun needs a declared amendment and a new directory. The suite job
`tf1_checks.py` reruns the fast exact gates and anchors the archived results, registered before the
stage-8 job, which stays last. The report maps every gate to its evidence and keeps the three
statuses apart. Earlier CWC-1, PF5, CL-1 and CL-2 evidence is unchanged.
