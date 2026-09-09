# A dynamical wave action for uniform frequency stretching

## Purpose and what is new

The recovered propagation models already derive redshift and event stretching from a prescribed photon Hamiltonian, and already expose conflicts with atomic/cavity standards. This pass does not repeat those ray calculations as a new discovery. It fills a narrower gap: derive that Hamiltonian from electromagnetic waves coupled to an explicit receiving field, and test whether the common frequency shift and energy accounting survive beyond the slowly varying approximation.

This is an optional propagation-field comparison. The current photon-scattering candidate still has unresolved color and image-sharpness problems, so capture improvements alone cannot complete it. A changing propagation field offers a different mathematical route, but it must satisfy the same matter, clock, distance and gravity requirements. No field or varying-rod law is adopted here.

Time-varying electromagnetic media provide mathematical/laboratory analogues for frequency conversion, such as [Qu et al., Theory of electromagnetic wave frequency upconversion in dynamic media](https://arxiv.org/abs/1804.07358). That work does not establish this fictional constitutive law. The specific matched wave action and checks below are our conditional derivation.

## 1. An explicit positive classical subsystem

In a specified rest frame, set reference c, epsilon and mu to one. Let n>0 and take

    L = K/2 [(partial_t n)^2-v_n^2 |grad n|^2] - V(n)
        + n E^2/2 - B^2/(2n),       K>0.

Thus epsilon=mu=n: the wave speed is 1/n while the impedance sqrt(mu/epsilon) is constant. The electromagnetic terms are gauge invariant because they depend on field strengths. They are not simply Z(n) F^2; that older operator leaves the principal light cone unchanged. This constitutive model selects a frame. A relativistic matter/vector/clock completion is still required; writing a preferred frame does not establish a covariant gravity theory.

The positive electromagnetic energy density and propagation-field energy are

    u_EM = n E^2/2 + B^2/(2n)
    u_n = K/2 [(partial_t n)^2+v_n^2 |grad n|^2] + V(n).

For V>=0 both are nonnegative in this classical domain. Varying n gives

    K(n_tt-v_n^2 Laplacian n)+V'(n) = E^2/2+B^2/(2n^2) = u_EM/n.

Without charges, the local exchange is

    partial_t u_EM + div(E cross B/n) = -(n_t/n)u_EM
    partial_t u_n + div(-K v_n^2 n_t grad n) = +(n_t/n)u_EM.

Increasing n therefore reduces electromagnetic energy and increases the field's energy. No energy disappears. With matter present, electromagnetic work on charges and any additional matter-field couplings must also be included. The homogeneous n field is not thereby identified as ordinary gravitons or as a population already captured in gravity wells.

## 2. The photon-like mode Hamiltonian now follows from the wave action

Restrict to a spatially homogeneous n(t), a finite set of transverse modes in a periodic volume, and V=0. For a mode coordinate A_k,

    L_k = n A_dot_k^2/2 - k^2 A_k^2/(2n)
    p_k = n A_dot_k
    H = P_n^2/(2K) + sum_k (p_k^2+k^2 A_k^2)/(2n).

Define I_k=(p_k^2+k^2 A_k^2)/(2k). Hamilton's equations give

    dI_k/dt=0
    H_EM=sum_k k I_k/n
    n_dot=P_n/K; P_n_dot=sum_k k I_k/n^2.

These are exact classical finite-mode statements. They reproduce the receiving-field energy used in the earlier [field–ray model](../../../fixed_atom_time_field/fixed_atom_time_field.md), now derived from the electromagnetic action instead of prescribing omega=k/n for photons.

Using optical time eta(t)=integral dt/n, each mode obeys an ordinary oscillator equation. Explicitly,

    A_k(t)=A_k(0) cos(k eta)+p_k(0) sin(k eta)/k
    p_k(t)=p_k(0) cos(k eta)-k A_k(0) sin(k eta).

Every mode's phase derivative is omega_k=k/n. All frequencies therefore scale by the same factor, while spatial wave number remains fixed. Here stretching means a longer temporal period relative to the selected frame; it does not mean a longer coordinate wavelength in homogeneous space.

Unlike a generic varying medium, this exactly matched constitutive family does not require adiabatic evolution to conserve I_k. This is a special property of epsilon=mu=n, not a statement about arbitrary dispersive or mismatched media. In a rapidly varying background the phase derivative is well defined, but a finite-time detector Fourier spectrum need not be a narrow stationary line. Detector response and real emission/absorption still require modeling.

For a formal finite-mode quantization, the canonical oscillator number operators commute with this mode Hamiltonian because its n dependence multiplies the fixed oscillator Hamiltonians. This does not solve vacuum-energy renormalization, continuum cutoffs, field fluctuations, n-domain boundary conditions or the complete quantum matter theory. The reported computations are classical.

## 3. Numerical verification and backreaction

Four modes with wave numbers 0.3, 1, 3 and 7 were evolved simultaneously with n and its conjugate momentum. Three cases vary field inertia and initial rolling speed. The code checks:

- individual mode actions against their initial values;
- every wave coordinate and momentum against the exact optical-time solution;
- total electromagnetic plus field energy;
- full wave-coupled n against an independently integrated reduced field equation.

The maximum total-energy residual is below 2e-11 in the chosen dimensionless units; wave-coordinate errors are below 2e-10. In the rapid-change example, the initial maximum parameter |n_dot|/k is 16.7, so this test is deliberately outside the ordinary slow-change regime. Exact mode-action conservation still holds numerically.

Radiation backreaction matters. With field inertia K=10 and initial n_dot=0.02, the field reaches n approximately 3.983 at t=20; the electromagnetic energy lost matches the field energy gained. The background is not freely prescribed after adding radiation. An affine n and the associated exponential distance law remain negligible-backreaction limits, as the earlier field–ray report already established.

For fixed matter clocks the old ray mapping and event-stretch relation follow from this dispersion. That conditional statement cannot establish measured redshift until the same action specifies source and observer standards.

## 4. Keep the charged-matter consequence attached

This electromagnetic action also changes static permittivity. With minimally coupled charges, fixed charge magnitudes and fixed electron/nuclear masses, leading Coulomb physics gives

    atomic length proportional to n
    atomic optical frequency proportional to n^-2.

These are the a=1,b=0 member of the already recovered [atom–wave–cavity test](../../../minimal_clock_interaction/derivation.md), not a new precision atomic calculation. Relativistic, nuclear and material corrections remain absent.

If an atom emits its transition at epoch e and the photon is observed at epoch o, propagation gives nu_received=nu_atom,e n_e/n_o. Comparison with today's same-species transition then yields

    1+z_measured = nu_atom,o/nu_received = n_e/n_o.

For increasing n this is a blueshift relative to current atomic standards, even though the photon energy decreased relative to the chosen coordinate-time frequency. Three numerical ratios check that distinction. Thus the favorable propagation result cannot be retained together with an unexplained assertion that atoms are unchanged.

Varying masses, changing rods, environmental differences or a frequency-dependent matter interaction are separate possibilities requiring their own action and observations. The earlier exact cavity/atom identity remains relevant. A changing-rod compensation is not adopted as an operationally nonexpanding solution, and exact clock protection must not be substituted for an actual measured likelihood.

## What this establishes

There is now an explicit classical electromagnetic-field model with positive energy, a receiving n field, exact common mode-frequency scaling and backreaction. It supplies a wave-level foundation for one recovered ray hypothesis. It also retains that hypothesis's unresolved matter-standard problem rather than hiding it.

This is not a completed conversion-to-graviton model: n has not been shown to produce travelling companions that seek wells or form stable deposits. Preferred-frame completion, real dispersion, atomic spectra, source history, spatial field dynamics, local tests, supported gravity and lensing remain necessary. All 20 tasks and 32 observational requirements remain in scope.

Run `python -X utf8 research_work/results/matched-wave/check_matched_wave.py`. Saved evidence is in `matched-wave-results.json`; the full runner includes it. These mathematical checks are not astronomical validation.
