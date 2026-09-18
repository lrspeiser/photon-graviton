# RUT-1 stage 6: a supported population constructed in self-consistent equilibrium, and the two repairs that must precede any stability claim

Declared before any stage 6 computation. Its own protocol file; stages 4 and 5 stay frozen, and
`rut4_checks.py` and `rut5.py` pin what they ran.

## Why this stage exists

Stage 5 established that a **cold** ring is collectively unstable with or without memory, and that the lag
slows the instability rather than causing it. The owner's direction is to stop trying to keep that ring
intact and instead ask which orbital populations this response can support:

> Constructed in self-consistent equilibrium, then tested for stability.

Two implementation defects, both found in the owner's review of ceb0c86, must be repaired first, because
**establishing that one mode grows is far easier than establishing that none does**, and stage 6's whole
point is the second kind of statement.

## Part 1, declared and run here

### R1. An eigensolver that locates every root in a declared region

The stage 5 contour count certifies nothing and is itself wrong: on the no-attraction control, whose
determinant is exactly `s²(s² + Ω²)` with no root in the open right half-plane, it returns one unstable
root at every angular mode, because the double root at the origin sits on the contour edge with its phase
winding under-resolved. A count mismatch therefore cannot distinguish a missed root from a miscount.

The replacement is a **contour-integral eigensolver** for the nonlinear eigenproblem `T(s)v = 0`, of the
kind Beyn introduced (arXiv:1003.1580): probe moments

    A_p = (1/2 pi i) contour_integral s^p T(s)^(-1) V ds,        p = 0 .. 2K-1

with a random probe `V`, assemble the block Hankel matrices of those moments, take the SVD to find the
numerical rank, and read the eigenvalues off the resulting small linear problem. It needs no initial
guess and returns every eigenvalue inside the contour, which is what a completeness claim requires.

The contour must stay inside the analyticity region: `B_m(s)` has poles at `s = -1/tau + i k Omega`, so
the search region is declared as `Re(s) >= eps > 0`, and roots are reported with their distance from the
contour. Roots are deduplicated before counting. Continuation is kept as a cross-check, not replaced.

### R2. A seeded verification that perturbs the whole state

`seed_mode` perturbs positions and velocities only. The state of this model is four things — positions,
velocities, the excitation field `E` and the force-producing field `C` — so that seed is not an eigenmode
of the full state and necessarily launches other branches. That, not an intrinsic ill-conditioning, is why
G4's declared form failed.

The repair builds the mode's own history into the fields. For a mode growing as `exp(s t)`, the field
perturbations implied by its past source history are

    delta_E(x,0) = integral_0^inf (1/tau_form) exp(-u/tau_form) delta_S(x,-u) du
    delta_C(x,0) = integral_0^inf [tau_keep/(tau_keep - tau_form)]
                                  [exp(-u/tau_keep) - exp(-u/tau_form)] delta_S(x,-u) du

which this stage realises by **prescribing the mode's trajectory over a spin-up window and advancing the
shipped field update through it**, then releasing the bodies. That constructs `delta_E` and `delta_C`
through the same code the run uses, and it replaces the analytic continuous-ring priming with the actual
discrete base field. The declared check is the owner's: the measured initial mode amplitude scales with
the seed, while the growth rate and frequency converge on the prediction.

### R3. A supported population constructed in self-consistent equilibrium

An axisymmetric distribution function of the conserved quantities of the stationary potential,

    f0 = F(orbital_energy, angular_momentum),    orbital_energy = |v|^2/2 + Phi_total,0
    rho0(x) = integral f0 d^2v

solved together with the field it writes, for writing proportional to mass,

    S0 = alpha (K convolved with rho0),   E0 = S0,   C0 = tau_keep S0,
    Phi_total,0 = Phi_external - C0

by iterating the population and the field to a common fixed point. Finite mass, with explicit tapering.
Because the kernel is a Gaussian and the state is axisymmetric, the convolution reduces to the ring
integral stage 1 already carries, so the solve is one-dimensional in radius.

Declared family: annuli of three widths and several velocity dispersions, spanning cold to warm, at the
stage 4/5 writing label so the support is comparable.

**Existence and stability are separate questions and are reported separately.** Part 1 answers only the
first, plus the stationary control below.

## Gates, with thresholds fixed now

| gate | requirement |
|---|---|
| H1 eigensolver against known spectra | exact quartic of the zero-lag rung, the Kepler case (which must return **zero** roots), and the transcendental `exp(-s) = 1/2` whose roots are `ln 2 + 2 pi i k`: every root located to 1e-8 relative and the count exactly right |
| H2 eigensolver robustness | doubling the quadrature nodes and shifting the contour edge moves each located root by < 1e-6 and changes no count |
| H3 eigensolver against continuation | on the 32-writer two-stage ladder, every root continuation found is found by the contour method to 1e-6; any root found by only one method is reported |
| H4 full-state seeded verification | growth rate and frequency within 5% of the linear prediction, at two seed amplitudes differing by 10x, with the measured initial amplitude scaling with the seed to 10% |
| H5 equilibrium self-consistency | the population and field satisfy the coupled equations with relative residual < 1e-6, at finite total mass with declared tapering |
| H6 stationary control | evolving the constructed equilibrium with no added perturbation for 20 reference periods, its drift in mean radius, support and radial dispersion is within 3x that of a no-memory control built by the same sampling |

A gate that fails is reported as a failure. No threshold is relaxed after seeing a number, and where a
declared measurement turns out not to resolve what it was meant to resolve, that is reported as a failure
too, with the substitute named as a substitute.

## Part 2, declared here and NOT run

* The complete linear mode analysis around the constructed equilibria, for the instantaneous, one-stage
  and two-stage responses at matched stationary fields, with located and counted unstable modes, growth
  rates, pattern frequencies and numerical uncertainty.
* Targeted nonlinear tests chosen from that analysis: one stable candidate, one borderline case, one known
  unstable control, each initialised with full-state perturbations and followed through saturation where
  applicable.
* A formation test: whether a population that survives the equilibrium tests is dynamically **accessible**
  from an empty field, not merely constructible.
* The reciprocal energy and angular-momentum accounting carried **through** perturbation, growth and
  saturation rather than applied afterwards, with the distribution of angular-momentum changes and not
  only its mean, and with transfers to the damping reservoir and to the external centre separated.
* Complex mode coefficients of `S`, `E` and `C` saved per checkpoint, and full restart states including
  the fields.

The success criterion for part 2 is declared now so it cannot be chosen later: **positive support with an
orbital distribution whose continuing evolution is bounded, or slow enough to be consistent with the
regime being claimed, without unaccounted steering or energy exchange.** Not "nothing moves radially". If
only an extreme velocity dispersion suppresses growth, that is a constraint to report, not a solution.
