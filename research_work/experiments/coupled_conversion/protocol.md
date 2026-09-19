# CWC-1: coupled conversion, clocks and wave response
Declared 19 September 2026 before execution. Baseline main: 0b9688d.
User authorization: comprehensive plan, provenance audit, then execution; retain
no dark matter, no expansion, fixed adopted observations, energy accounting,
and regular main checkpoints.

## Objective and decision
Test whether one explicitly coupled, initially empty receiving-wave model can
produce photon energy transfer, measured spectral and event timing changes,
matter response and optical bending together. This is a comprehensive first
coupled *model-family test*, not a promise to derive a complete gravity theory.
No measured galaxy or cluster profile is fitted by inserting unseen matter.
No result from a different model is imported as a success of this model.

The tested family extends the repository's matched-medium electromagnetic
action. That action and the general mathematics have known antecedents; see
provenance.md. No equation is certified historically novel. Our contribution
here is the explicit model specialization, implementation and falsifiable
cross-check campaign. A scalar receiving wave is not an ordinary graviton.

## Model, units, and conserved energy
Use static Cartesian coordinates, c_ref=1 and reference length=1. No scale
factor, cosmological distances, primordial population or fixed cosmic age.
One polarization A of a transverse electromagnetic field has canonical D.
The receiving scalar has coordinate phi and momentum Pi; initially phi=Pi=0.
Use a=exp(-g phi), n=1/a, corresponding to epsilon=mu=n.
For one-dimensional waves, or a two-dimensional TE wave A_z,

H = integral [ a (D^2+|grad A|^2)/2
               + Pi^2/2 + v_phi^2 |grad phi|^2/2
               + m_phi^2 phi^2/2 + lambda phi^4/4 ] dV
    + sum_j [ |P_j|^2/(2 M_j) + I_j Omega0 exp(-b g phi(X_j)) ].

All displayed energy terms are nonnegative. I_j>=0 is a conserved internal
oscillator action; its conjugate angle defines the modeled clock. M_j>0.
These are effective material oscillators, not a solved atomic/fermion theory.
Their baseline internal energy and its changes must be counted. Finite photon
energy is prepared in initial EM packets; there is no external pump after t=0.
No field profile, independent companion population or hidden source is inserted.

Clock rate Omega_j/Omega0=exp(-b g phi_j) and force
Pdot_j=b g I_j Omega_j grad phi_j follow from the same material term.
The corresponding source b g I_j Omega_j is included in the scalar equation.
The EM source is g u_EM; the scalar transfers energy reciprocally with EM and
material clocks. Moving material positions and momenta are evolved. Omitting
a reciprocal source is an explicitly rejected negative control.

b=0 is an assumed protected-clock completion, not derived atomic protection;
b=1 shares the optical rate; b=2 is the leading fixed-charge/mass Coulomb clock
scaling of the matched medium, used as an effective oscillator comparison.
It is not a complete atomic action. No completion is privileged after outcomes.
The exponential optical constitutive law is not Bekenstein's Z(phi)F^2 law,
and no claim is made that our scalar has spin2 or Gertsenshtein conversion.

No universal gravity equation has been derived from this Hamiltonian. Its
material force and ray bending are tested as candidate *extra* responses.
They cannot silently be identified with a general-relativistic metric or an
observed mass density. The ordinary matter oscillators are finite toy probes,
not galaxy mass reconstructions. A fitted coupling that varies by object is
forbidden.

## Stage 0: provenance and claim audit
- Inspect prior project formulations and primary publications on dynamic media,
  scalar/electromagnetic coupling, photon/graviton conversion and splitting.
- Record each component as established mathematics, existing project result,
  proposed specialization, or unresolved identification.
- Explicitly distinguish whole-mode photon/graviton mixing from partial
  energy transfer into a newly named receiving wave.
- Keep search queries, source URLs, dates, access limits and a claim register.
- Do not call a literature search proof of first invention or copy external code.

## Stage 1: derivation and independently checked small systems
1. Derive reciprocal EM, scalar, clock and matter equations from H.
2. Check Hamiltonian directional derivatives with finite differences.
3. Homogeneous exact-mode reduction: actions constant; EM energy A exp(-g phi).
   Evolve phi and finite clock energies from zero receiving field with
   g={0.05,0.2,0.5}, b={0,1,2}, m={0,0.4,1}, lambda=0.1.
4. Calculate ray arrivals, carrier ratios and clock-normalized timing from
   the same solution, including zero EM, zero coupling and zero clock energy.
5. Derive composition dependence of force and the fixed-ruler local-speed
   identity; distinguish numerical agreement from physical acceptability.

Numerical gates: homogeneous energy error<1e-9, independent clock/arrival
identity<1e-6, Hamiltonian gradient relative error<1e-5, nonnegative energies.
Scientific gates: detectable positive measured redshift>0.001; event/spectral
stretch agree within1%; fixed-ruler local optical speed differs from its
reference by<0.1%; force per inertial mass differs by<0.1% for two material
oscillators with a factor2 internal-action/mass ratio. A finite source-off
response must be reported; photon-exclusive generation requires zero response
to1e-10, not an ignored material energy source.

## Stage 2: spatial EM-to-wave conversion and pulse diagnostics
Implement the one-dimensional Hamiltonian on a periodic finite lattice.
Periodic boundaries close the energy ledger; analyze direct-arrival windows
before a packet circles the box. Edge-averaged positive a weights |grad A|^2.
Each split subflow includes its exact reciprocal momentum update. Use a
symmetric second-order composition; no damping, floors or uncounted filtering.

Primary: L=48, receiving speed1, m=0.4, lambda0.1, g0.2, b2.
Two initial right-going Gaussian packets, centers -14 and -10, width0.8,
wavelength2, total EM energy1; fixed receiver x=6, end time26.
Four mobile material oscillators initially at -1,-1,+1,+1, unit inertial
masses, internal energies0.0005,0.001,0.0005,0.001. Initial material momenta0.
These are prepared incoming pulses, not a reconstructed stellar emission history.

Use grids512/1024/2048 with dt<=0.08 dx (adapt end step only). Refine time
separately at fixed grid. Compare g={0.05,0.2,0.5}, b={0,1,2} at grid1024;
mass0/1 and receiving-speed0.5 sensitivities; wavelengths1/2/4 at equal
initial EM energy; total EM energy0.25/1/4 at the primary carrier. Preserve all.
Controls: zero EM, zero g, zero clock action, reversed pulse direction,
and an isolated free receiving-wave packet (only this explicitly labeled
transport control starts with receiving energy).

Save sector histories, field profiles, receiver time series, material recoil,
momentum diagnostic, clock rates, pulse arrival centroids, carrier phase/
frequency diagnostics, distortions and right/left-wave energy proxies.
Compare frequency and arrival changes with the same-grid uncoupled control,
not a differently fitted background. Receiver measurements are normalized by
the predicted material clock; no asserted universal clock slowdown.
Measure any loss of EM energy, energy return, and actual receiving energy
separately; a redshift alone is not a companion-production current.

Numerical gates: primary refined relative energy residual<1e-3 and improves
under separate timestep refinement; receiving-energy and direct-arrival
centroid changes on the two highest spatial grids<2% (energy normalized to
initial EM energy, arrival to one reference light-crossing unit); no-source/
no-coupling controls match their stated equations; discrete momentum error
normalized to initial total EM energy<2% on the refined primary. The lattice
has only discrete translation symmetry, so momentum must be checked rather
than claimed exactly conserved. Free receiving-wave energy error<1e-3.
Scientific tests: primary net photon-to-wave transfer>0.1%; clock-corrected
positive redshift and timing agreement as above; equal-energy carrier-color
fractional shift spread<1%; shape/sideband changes and intensity dependence
reported even on failure. Wave action is a diagnostic, not exact photon number.

## Stage 3: generated spatial response and optical deflection
Use the same Hamiltonian in 2D TE, L=32, four converging EM packets initially
at (+-6,0),(0,+-6), longitudinal width0.8, transverse width1.5, wavelength3,
total EM energy1. Same g,m,lambda,v,b and finite material energies; four
material oscillators at (+-1,0),(0,+-1). Zero receiving field. End time18,
grids128/192/256, dt<=0.08 dx/sqrt(2). Save snapshots near t=6,10,14,18.
Compare b0/b1/b2, and a source-off control. No pre-imposed central field.

From each generated final snapshot compute negligible-energy Hamiltonian rays
with H_ray=exp(-g phi)|p|, impact parameters+-0.5,+-1,+-2, entering x=-8 and
leaving x=+8. This is an explicitly *frozen-snapshot* optical diagnostic;
it does not certify propagation in the evolving field or real3D lensing.
Compare with independent weak-gradient line integration, ray-step refinement,
mirror symmetry and zero-field straight rays. Evaluate the material force
from the same phi and report its composition dependence.
Measure receiving energy remaining within radius3 after the light leaves.
This is finite-duration retention, not permanent capture or stability.

Numerical gates: refined energy residual<1e-3; generated deflection changes
between top grids<5% of the maximum bend scale; ray tolerance/step change<1e-4
relative to the maximum bend; weak-gradient comparison<2% when |g phi|<0.02;
mirror error<1% on the same scale. Scientific criterion: attractive, source-
funded nonzero bending and matter response with shared parameters and
composition-independent acceleration. Report failure or mixed sign explicitly.

## Stage 4: observation and gravity readiness audit
Execute a structured audit against archived measured targets: redshift/event
timing, stellar motions, Coma shear and cluster ordinary-light/gas profiles.
Record what the candidate actually predicts and what it cannot yet predict.
No arbitrary mapping of the dimensionless2D fixture to a galaxy, fitted
unseen density or added metric potential is allowed.
Required promotion gates: microscopic/gauge/spin identification; a consistent
matter/rod/clock completion; three-dimensional source/field/stress equations;
derived dimensional normalization; bolometric histories and formation;
joint motion/lensing prediction with fixed distances and shared coefficients.
When these are absent, mark not-ready/failed prerequisites and do not invent
an observational chi-square. This gate evaluation is an executed test of
readiness, not a claim that missing observational fits were run.
Supply diagnostics use actual finite delivered energy in the simulations.
No assumed cosmic age or universal finite luminosity history excludes the
broader idea.

## Evidence, reporting and checkpoints
Commit this protocol and provenance before calculations. Publish main after
the small-system/1D stage, after2D stage, and after final review; fetch before
every push, never force-push. Keep failed numerical attempts and amendments.
Write source/input hashes and runtime versions before a run; immutable output
directories. The final report maps every stage to evidence and separates
numerical, mechanism and observational status. No full historical suite whose
comparisons activate excluded premises is run. Earlier PF5 and RUT outcomes
remain unchanged. Passing this family does not establish historical novelty.
