# Evolving propagation field (PF-1): one wave law for redshift, event stretching, photon number, brightness and energy exchange

Declared before execution, 13 September 2026. Baseline: `main` at 636a906. Opened at the project owner's request, holding observations fixed rather than the assumption that companions must become slow particles through collisions with ordinary matter. This is a proposed model; none of its postulates is established physics.

## Postulates

- **P1. Homogeneous index.** A homogeneous propagation index n(t) > 0, measured against fixed material clocks and rulers. In this proxy the matter sector does not couple to n. What n is operationally relative to, including any atomic coupling, stays an open requirement.
- **P2. Wave law.** Free waves obey d/dt(n dA/dt) - (c^2/n) Lap A = 0. Electromagnetic and gravitational waves share this law (common propagation). A photon-only variant is not adopted.
- **P3. The index is dynamical.** L_n = (M/2) ndot^2 - V(n), coupled to the wave Lagrangian L_w = (n/2) A_t^2 - (c^2/2n)|grad A|^2. Total energy is conserved by construction.

## Consequences to derive exactly and verify

Define tau = int dt/n; P2 is then the ordinary wave equation in tau. From this:

- 1+z = n(t_o)/n(t_e), and dt_o/dt_e = 1+z for fixed separation.
- Wave energy scales as 1/n, and photon number (energy/frequency) is conserved.
- Bolometric flux F = L/(4 pi D^2 (1+z)^2) in fixed Euclidean geometry.
- A linear index gives 1+z = exp(alpha D), with alpha = ndot/c.
- The wave sector loses energy at Q = (ndot/n) u_gamma, and the field gains it. The wave acts on n with dL_w/dn = E_w/n exactly.
- A thermal mode occupation stays thermal, with T proportional to 1/n.

## Numerical tests (tolerances declared now)

- **T1. Driven finite-difference simulation.** Integrate P2 in physical time, not in tau. A material-clock emitter drives Gaussian-envelope pulses at two times; an observer sits at separation D; the index is linear with z = 1.
  - The carrier-frequency ratio, the envelope-width ratio and the pulse-spacing ratio must each equal 1+z within 0.5%.
  - Wave energy times n must stay constant within 0.5% after emission.
  - Control with constant n: all three ratios equal 1 within 0.2%.
- **T2. Coupled Fourier-mode integration** with a dynamical index:
  - total energy (field plus waves) conserved to 1e-9 relative;
  - each mode's photon number conserved to 1e-9;
  - dE_w/dt = -(ndot/n) E_w, with the field gaining exactly what the waves lose;
  - the reduced mechanics M n'' = -V'(n) + E_tau/n^2 reproducing n(t) to 1e-8.
- **T3. Analytic observational consequences.**
  - Report the redshift-distance law and time-dilation exponent b = 1.
  - Report D_L = (1+z) ln(1+z)/alpha.
  - Tabulate distance-modulus differences against a flat FLRW comparator (Omega_m = 0.3) with the same low-redshift slope, at z = 0.1, 0.5, 1 and 1.5.
- **T4. Conditional exposed-data brightness test.** If archived supernova distance moduli and covariance are available in the repository, fit the single alpha together with the usual absolute-magnitude nuisance, and report residuals by redshift bin. Fit the FLRW comparator with the same data, covariance and nuisance freedom. These data are exposed, so this is not a blind test.

## Not claimed

- An operational clock completion, including why atomic transitions do not share n.
- A positive, eternal history for n: the linear form is an interval model only.
- Inhomogeneous n: geometric optics, lensing and local gravity.
- The origin of the microwave background: P2 can preserve a thermal spectrum but does not create one.
- Any connection to galaxy gravity, which is left to the RPG-1 branch.

## Assessment rule

PF-1 counts as an internally consistent propagation law only if T1–T3 pass. Its observational standing is read from T4 and the time-dilation comparison, with every limitation stated. No promotion over the archived conversion law or over expansion follows from exposed data alone.
