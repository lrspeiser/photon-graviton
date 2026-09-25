# Main audit and the energy-shift direction

Research review and bounded calculation, 25 September 2026.

## Provenance and limits

The GitHub main branch was read at commit `67075c2c59dbe3c7cf55af4eb1d9747bcf000ffb` (25 September 2026, 13:31:33 UTC). Reviewed source files include STEP-BACK-AUDIT.md, RULES.md, the hot-companion technical README sections 29.3-29.5, and code/wave_dark_v15.py and code/receivers_v15.py. This is not a claim to have independently re-audited every historical result or rerun the astronomical suite. No repository write was made.

The earlier separate `structured_medium_derivation` package supplies the effective exchange kernel and its range limitations. The new computation here diagonalizes fixed-excitation sectors of that kernel at finite occupation. It does not rerun the earlier wave/bath simulations or implement motion, new state preparation, radiation, screening, universal inertia or lensing. It is an effective weak-coupling calculation, not a new gravitational theory.

## Assessment

The energy-shift route is worth continuing. It removes the need to identify every static force with recoil from steady emission and permits attraction with no mean classical oscillator amplitude. However it requires correlations, does not yet give the required distance/mass law, and is not the same thing as an intensity-only receiver.

The audit usefully identifies the shared-clock assumption as model-specific. Its claims of general exclusions should be narrowed further:

- Force equals fed power divided by phase speed is a result for the studied travelling-wave recoil channel, not all conservative or wave-mediated forces.
- A finite kinetic-fuel depletion estimate assumes the adopted power law and source of that power. It is not a cost imposed on every conservative force.
- The specific full warm-force tests encode motion through internal mixing, but wave_dark_v15.py explicitly advances positions, recalculates couplings and redraws velocities. The missing result is a complete high-heat moving-source-and-receiver evolution with mechanical backreaction, not that no positions have ever moved.
- receivers_v15.py already includes passive and below-threshold amplifier controls responding as intensity. The exact proposed total-energy/flow receiver has not been constructed, but intensity-sensitive responses are not wholly untested.
- Project scope rules are not observational exclusion theorems. Existing failures are constraints on the adopted empirical law, even when input assumptions may contribute.

## New finite-occupation calculation

For eight identical two-level pieces with fixed seeded irregular positions, set

    H_n = omega_a n + sum_{i<j} J_ij (s_i^+ s_j^- + s_j^+ s_i^-),
    J_ij = -C exp(-kappa r_ij)/r_ij,
    C = 0.04/(4 pi), kappa = 1.3537825795626388, omega_a=1.2.

Each total excitation number n=0,...,8 is treated separately. Diagonalization selects the lowest energy within that sector; its dynamical preparation has NOT been demonstrated by this new check. Unlike the earlier one-excitation extension, sectors n=2,3,4,... are included here. No coefficient is refitted between sectors.

The force follows from the operator -partial H/partial X, evaluated in the selected eigenstate. The force agrees with an independent finite difference of that eigenenergy (Hellmann-Feynman). In a nonstationary state one must use the force operator; differentiating a state-dependent energy while forgetting changes of the populations can hide work/heat.

| Total excitations | RMS internal force relative to n=1 |
|---|---:|
| 0 | 0.0000 |
| 1 | 1.0000 |
| 2 | 1.7173 |
| 3 | 2.1488 |
| 4 | 2.2930 |
| 5 | 2.1488 |
| 6 | 1.7173 |
| 7 | 1.0000 |
| 8 | 0.0000 |

At every partially filled sector all 28 pair correlations entering the force are positive and the corresponding pair contributions are attractive. The full-cloud resultant is zero to floating-point accuracy. Every one-point <s_i^-> is zero: the state has a definite total excitation number. This removes a nonzero mean oscillating field, not pair correlations or all quantum coherence. It is not an intensity-only detector.

Zero at n=0 or n=8 means zero leading resonant exchange force in this effective Hamiltonian. It does NOT mean the full microscopic theory has no force: ground-state dispersion, higher-order contributions, radiative decay and other channels are outside this calculation. Full inversion blocks the elementary exchange because there is no unexcited site to accept it. This is a concrete reason not to assume that "more internal energy always means more attraction." The particle-hole symmetry comes from identical two-level exchange, not a general material law.

Validation in this run:

- maximum absolute force-versus-energy-derivative difference: 6.3931e-13 model force units;
- largest net internal force: 1.2644e-18;
- mean one-point oscillator amplitude: exactly zero in the sector construction.

This is an internal-force diagnostic at fixed mass/positions and changing excitation population. It is NOT a test of source mass scaling or external gravitational acceleration.

## Why sensing energy alone does not yet give the desired law

Let the local intensity outside an isolated source be I=K M/r^2, with K constant. Suppose the receiver contributes a universal local energy V=f(I), with no other independent field, gradient or history variable. Then

    F_r=-dV/dr=2 I f'(I)/r.

The desired magnitude F proportional to sqrt(M)/r would require 2 I f'(I) proportional to sqrt(M). At two positions/source masses with the same I but different M, the left side is identical and the right side is not. Therefore this restricted local-energy-only prescription cannot recover both exponents over independently varied M and r.

For f(I)=I^q, magnitude scales as M^q/r^(2q+1). Examples (sign is a separate microscopic question):

| Energy function | Mass exponent | Inverse-distance exponent |
|---|---:|---:|
| I | 1 | 3 |
| sqrt(I) | 1/2 | 2 |
| log(I) | 0 | 1 |

These algebraic exponents were checked with finite differences. This is NOT a no-go result for nonlocal, nonequilibrium, gradient-sensitive, collective or multichannel media.

For a direct momentum-transfer receiver instead, |F|=sigma_eff(I) I/v in a single travelling channel. Recovering sqrt(I) requires sigma_eff proportional to I^-1/2 over the operating range. This may be an emergent intermediate response of nonlinear matter, but it must be derived and its low-intensity crossover, power source and many-direction response must be calculated. A generic linear intensity detector is insufficient.

## Why retuning the previous bandgap alone is insufficient

For the fixed-correlation Yukawa force,

    p(r)=-d ln|F|/d ln r=2+(kappa r)^2/(1+kappa r)>=2.

For any positive, radius-independent weighted mixture,

    |F|=r^-2 integral w(kappa) exp(-kappa r)(1+kappa r) d kappa,
    d[r^2 |F|]/dr=-r integral w(kappa) kappa^2 exp(-kappa r) d kappa <=0.

A numerical mixture of 501 positive weights and 201 radii confirms this bound (minimum finite-interval exponent 2.01117 in that example). The analytic inequality is the actual result; the numerical example does not establish it generally. If correlations, occupation or the medium response vary with configuration, the fixed-weight assumptions no longer apply. That is a genuinely different next construction and needs its own dynamics and energy accounting.

## Next discriminating experiment

1. Retain the phase-locked source/receiver model and the linear structured-medium model as controls. Do not adopt the old empirical force in a new microscopic code.
2. Give all pieces a common local Hamiltonian, finite internal occupation and actual movable centres with kinetic energy. Compute the force as -Tr[rho partial H/partial X]. Track the work of any imposed trajectories or collision baths.
3. Add one physically motivated, energy-bounded medium nonlinearity at a time, for example a positive local anharmonic term beta X^4/4. Measure whether population changes only the transition frequency or changes the long-wavelength, configuration-dependent collective response. A generic quartic perturbation is not expected automatically to create a new asymptotic force exponent.
4. Test actual measured velocity dispersion from the evolved motion. Compare cold, freely moving, collisional and rigidly moving controls with identical initial internal fuel. Extend to the dimensionless cluster regime only after the mapping from velocity to response has been derived.
5. Measure U (energy density), J (energy flux), stress and two-point correlations separately. Use an exterior probe, two opposing fluxes and a central weak binary enclosed by a hot shell. Isotropy requires no net preferred force at the exact shell centre, but allows a changed response to a perturbation. Do not equate scalar energy with vector flux.
6. Independently vary physical source mass, resolution, source size, receiver distance and occupation. A new candidate must show a controlled interval approaching source mass exponent 1/2 and inverse-distance exponent 1, or propose a different observable law and test it honestly.
7. In parallel, remove the X-ray target-dependence in the cluster heat inputs and test the frozen velocity-dispersion lensing prediction. Keep no-hold screening and alternative static distance laws as registered comparisons under the same evidence standard, not automatic cures or silent changes to the baseline.

Stop/redirect if the nonlinearity only changes a Yukawa length, gives an occupation-dependent force inconsistent across material states, or produces the desired signal only through external work omitted from the ledger. A failure narrows that specific construction, not all nonlinear media.

## References used to bound broad claims

- Jaffe, Physical Review D 72, 021301 (2005), DOI 10.1103/PhysRevD.72.021301: Casimir interactions as forces of charges/currents; conservative attraction is not synonymous with continuous fuel consumption.
- Chen et al., Nature Photonics 5, 531-534 (2011), DOI 10.1038/nphoton.2011.153: backward scattering force through directional momentum redistribution, not a general requirement of a gain-bearing receiver.
- The supplied structured-medium README cites Hopfield, Notararigo/Passante/Rizzuto, and harmonic-chain reservoir literature; it explicitly distinguishes the new calculation from established ingredients.
- Any emergent nonlinear modified-Poisson law needs an equation-level comparison with existing theories, including Bekenstein and Milgrom (1984), DOI 10.1086/162570. Simply choosing the required nonlinear field energy would not be an independent microscopic derivation. The repository's RULES.md has not been changed.

## Reproduce

Python with NumPy and SciPy:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python checks.py

This recomputes results.json. There are no data downloads or astrophysical inputs. The directory is separate from the GitHub repository. Finite-difference errors depend weakly on platform/library precision; physical results here are model-unit diagnostics. These checks are not new astronomical evidence and not a probability that the theory is correct.
