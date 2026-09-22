# JR-6: source currents, orientation and three-dimensional companion structure

21 September 2026 (Pacific time). Exploratory toy construction, not a new real-galaxy fit or a complete gravity theory. Baseline inspected: 38573a58481bd27fa092fc13af1c2fe137171679. The preceding JR-5 gas-exchange evidence is unchanged. Code was run locally before being published; no GitHub Actions execution is claimed.

## Theory first

Matter and radiation may generate and exchange energy with a companion state. Its distribution, currents and history must eventually supply a connected response for matter and light while preserving ordinary internal planetary dynamics. This calculation explores relative motion rather than making the whole theory another static acceleration formula. It is one transport implementation within the broader memory/mass-current/multisector portfolio, not a requirement that all gravity obey diffusion or a single metric.

The user's ocean-current analogy motivates a map of moving components. Distinguish orbital motion, an object's own spin, and directional illumination. A round isotropic source has no face. A stellar bar's pattern speed need not equal the orbital speed of its stars. Mean velocity is also not enough: equal opposite streams can have zero mean current but nonzero velocity variance and kinetic energy.

## A concrete 3D construction

Solve in cylindrical coordinates R, phi, z:

    P_t = Dp lap(P) - kp P + km C - lambda_p P + J
    C_t + Omega_c(R) C_phi = Dc lap(C) + kp P - km C - lambda_c C

P and C are radiation-like and companion energy densities. Gas density is axisymmetric and identical in the primary cases. The two-lobed illumination pattern rotates at Omega_s, and the companion is hypothesized to be entrained by a gas-related angular current Omega_c(R). Entrainment is prescribed, not microscopically derived.

    rho_g = exp(-R/2) sech^2(z/0.35)
    h = rho_g/(rho_g+0.15)
    kp = 0.5 exp(2h), km = 0.5 exp(4h^2)
    J = J0(R,z)[1 + 0.6 R^2/(R^2+0.8^2) cos(2(phi-Omega_s t))]

J0 is a Gaussian of radial width 1.5, vertical width 0.35, normalized to total luminosity one. Dp=0.5, Dc=0.05, lambda_p=0.5, lambda_c=0.1 and Omega_s=0.6. R extends to 6 and z to +/-3 with reflecting boundaries. All values are toy units, not fitted astronomical constants.

Both exchange directions are present. Source power is withdrawn from a fuel account and the lambda losses go to an outgoing-energy account:

    d(Ep+Ec)/dt = Lsource - lambda_p Ep - lambda_c Ec.

This conserves channel energy including supply and escape. It is not a complete gravitational binding-energy or gas/stellar recoil, momentum and torque ledger. Diffusion is an effective transport approximation, not photon propagation at c. No direct magnetic or spin coupling of the companion is established.

The solver obtains the long-time periodic state through the m=0 and complex m=2 angular harmonics, retaining the R and z dependence. Full 3D fields are reconstructed. This is not a 2D mass sheet, but it is also not a primary empty-field evolution; one small-grid transient is independently evolved as a control.

## Executed campaign

Nine primary cases: aligned, opposed, stationary companion, differential angular flow, both signs reversed, two axisymmetric-source controls, and short/long companion-residence sensitivities. Four fine-grid cases compare aligned/opposed currents at 48x64 and 72x96 radial/vertical resolution, with 128/192 angular readout samples. One expanded-boundary case maintains comparable resolution. Each case solves two harmonics. No galaxy residuals were used to choose the parameters.

All force readouts come from a fixed softened three-dimensional potential of the resulting source. The face-on deflection uses the matching line-integrated potential. G=c=1 and Plummer softening=0.2; there is no separate light-deflection gain. Both companion-only and P+C source readouts are retained. These are conditional toy gravitational readouts, not measured galaxy forces or a derived microscopic energy-to-gravity relation.

## Result: currents reorganize the source without adding energy

Define the normalized two-lobed spatial moment

    qC = |integral C R^2 exp(2i phi) dV| / integral C R^2 dV.

Fine-grid results:

| Source/companion relative motion | qC | Companion-axis lag behind source | Maximum sideways force / mean inward force at R=2.5 |
|---|---:|---:|---:|
| Aligned, Omega_c=+0.6 | 0.0712105068 | 30.5519 degrees | 12.8962% |
| Opposed, Omega_c=-0.6 | 0.0159067682 | 55.8829 degrees | 4.4016% |

The aligned case has 4.47674 times the normalized directional moment. This is NOT 4.48 times more gravity or quantum coherence. The azimuthal half-range of inward-force variation is 1.9204% versus 1.6435%; the corresponding deflection half-range is 1.3162% versus 0.3682%. Different spatial readouts do not scale identically.

When BOTH energy channels are counted, the maximum sideways/inward ratios are 12.2104% aligned and 6.2729% opposed. The effect is therefore not produced solely by relabeling energy and ignoring the other channel.

Within each grid the primary aligned, opposed, stationary-companion and differential-flow cases have identical total P and C energies and identical azimuthally averaged radial force. At the finest grid Ep=1.7045404724 and Ec=1.4772976380. The supplied luminosity equals outgoing power one. Reversing both currents mirrors the field exactly; making the source axisymmetric removes the co/counter distinction.

### An exact limitation of this construction

The gas/rates are axisymmetric and the transport equations are linear. The m=0 and m=2 modes therefore do not mix. Changing angular currents can change the directional pattern and lag, but cannot change the mean source or mean radial force here. This toy alone does not repair the old overpredicted mean rotation curves. Nonaxisymmetric exchange, nonlinear saturation, radial/vertical redistribution or matter backreaction could connect orientation to the mean response; that connection still needs construction.

For the scalar-potential readout used here, a circulating energy current is compatible with curl(g)=0. A sideways force relative to the galactic center is not proof of a nonconservative vortex force. Velocity-dependent extra forces require their own coupling and energy/torque account.

## The useful reduced equation

For one ring,

    C_t + Omega_c C_phi = D/R^2 C_phiphi - C/tau + S(phi-Omega_s t).

A source harmonic m produces

    Cm = Sm / [1/tau + D m^2/R^2 + i m(Omega_c-Omega_s)].

Let tau_eff=1/(1/tau+D m^2/R^2) and X=m*abs(Omega_c-Omega_s)*tau_eff. Relative to a co-moving pattern, the amplitude is 1/sqrt(1+X^2), with lag magnitude arctan(X). Relative passage and persistence are linked in this equation; neither an arbitrary global clockwise sign nor independently fitted lags are needed. This is a consequence of the specified transport model, not a universal gravity law. Instantaneous local equilibration removes this particular memory lag while leaving static geometry effects possible.

## Other ingredients worth incorporating

| Component | Explicit information | Candidate mechanism; not yet established companion physics |
|---|---|---|
| Distributed stellar populations | Position, luminosity/spectrum, age and orbital velocity | Where power enters gas and how moving illumination is retained |
| Co/counter streams and shear | Separate stellar/gas velocities and velocity dispersion tensors | Residence, repeated encounters, winding and response-angle changes |
| Stellar spin and facing | Spin axis/rate, surface-brightness multipoles, magnetic/emission axes | Rotating anisotropic illumination and periodic response |
| Gas phases and topology | Thickness, clouds, holes, ionization, temperature and turbulence | Conversion/reversal, escape routes and compression; not one gas-fraction number |
| Bars, spirals and warped disks | Pattern speeds, radius-dependent tilt and gas-star offsets | Long-lived repeated exposure and off-plane structure |
| Dust and magnetic geometry | Frequency-dependent opacity, reradiation, field direction | Known light/plasma transport modifies sources; direct neutral-companion coupling requires derivation |
| Binaries, planets and moons | Orbit/spin axes, luminosity, intercepted area, tides and periodicity | Local receivers and periodic perturbers; galaxy-scale importance needs a derived weighting |
| Winds, jets and explosions | Power, momentum, opening angles and duty cycles | Cavities, asymmetric source injection and transient currents |
| Environment/history | Inflow, tides, encounters, old companion/current state | Offset structures, reorientation and boundary flux |

Gas may be treated as a fluid where appropriate; stars need particles or a phase-space distribution rather than being made viscous water. The moving gas, star, radiation and companion channels need not share one velocity. Every new exchange needs an explicit reciprocal energy/momentum rule.

### Stellar and planetary orientation must survive averaging

A simple axis-related pattern has Q=sum w_i(n_i*n_i^T-I/3)/sum w_i. Flipping spin along the same axis leaves this tensor unchanged. To distinguish spin sense, add an actual current, spin-vorticity relation, directional polarization or phase-sensitive response, with an explicit coupling.

For independent isotropic equal-weight axes, the RMS anisotropy relative to full alignment is exactly 1/sqrt(N). A separate 100-realization Monte Carlo check at N=10000 gives 0.9828%, versus the analytic 1%. For unequal weights use Neff=(sum w)^2/sum w^2. This is classical orientation averaging, not a claim that real stellar axes are independent or that their radiation is phase-coherent.

Fast spin or planetary forcing can average away if its frequency is too high for the field response. Correlated orientations or a derived resonance could preserve an effect, but neither follows just from having many objects. Weight planets by physical mass, available power, cross section or angular momentum, not by object count. For geometric interception of isotropic stellar light, the power fraction is Rp^2/(4a^2). Nongeometric strong coupling would require a mechanism, not an exception inserted into the fit.

## Numerical evidence and limits

- Primary/fine source-outgoing energy residual below 8.1e-15; linear residual below 5.5e-14; both channel densities positive at every azimuth.
- Fine-grid force changes at most 0.0342%, deflection at most 0.0181%.
- Enlarging the domain changes force by up to 0.792%, deflection by 0.110%. The GLOBAL normalized quadrupole changes 5.52% because distant material receives R^2 weight. It is not an infinite-domain invariant.
- An independent full R-phi-z finite-difference calculation uses 2560, 5120 and 10240 unknowns. Directional-field error decreases 6.05% -> 1.49% -> 0.371% with angular refinement; its angular mean agrees to 5.4e-16.
- A small-grid, independently propagated empty m=2 transient approaches the periodic state by toy time 80 to relative discrepancy 7.7e-16.
- 85 numerical/readback checks, zero failures. This is not 85 independent physical successes.

The next integrated construction is a moving distributed stellar source plus nonaxisymmetric, clumpy, tilted gas and a dynamically justified companion current. It should output full projected velocity fields and lensing structure. One-dimensional rotation curves do not determine 3D thickness, warp, spin or angular histories. No new maps, galaxy fits, clusters, stellar spins or planets were simulated here.

## Reproduce and provenance

From this directory, with a fresh output path:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python currents.py --output fresh-run
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit.py --results fresh-run

Executed locally on Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0. Full arrays, logs, the pre-execution local protocol, code and detailed report are in the originating conversation's JR-6 package. The local protocol was not claimed to have been published before execution. Original files are unchanged.

Exact executed source SHA-256:
- currents.py: 3d0d4649ecb5ad1e6b11f71707656e909b662fd84fe3035cd22a07630fa67455
- audit.py: bf402bc243ff5fb1054d14c85d042755235dd13f2f211bca94f601528bc1d304
- Full results JSON: 7709b0828e1620b46234d1ae1c457876baf3d292a603e7370107beb7d4a3b814
- Full audit JSON: c159d8219e1a586e809dfb3f0de1be29f0a9a864915c0f8ca11b720f653cb84f

## External primary-source context, not evidence of companion conversion

- Espinosa Lara and Rieutord (2011), Gravity darkening in rotating stars: https://arxiv.org/abs/1109.3038 . Rotation can change stellar emission with latitude.
- Capelo and Dotti, Shocks and angular momentum flips: https://arxiv.org/abs/1610.08507 . Gas shocks and stellar motion can behave differently in mergers.
- Beom et al. (2024), SDSS IV MaNGA, opposite gas angular momentum: https://arxiv.org/html/2410.06256v1 . Observed co/counter-rotating samples exist, with morphological and population differences to retain.
- Effenberger et al. (2012), anisotropic cosmic-ray diffusion: https://arxiv.org/html/1210.1423v1 . Directional charged-particle transport, not a direct neutral-chi magnetic interaction.
- Gaia Collaboration (2023), Mapping the asymmetric disc: https://www.aanda.org/articles/aa/full_html/2023/06/aa43797-22/aa43797-22.html .
- Walter et al. (2008), THINGS: https://arxiv.org/abs/0810.2125 . Resolved gas maps and kinematics; no new survey data ingested here.
