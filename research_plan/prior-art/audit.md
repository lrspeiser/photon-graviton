# Equation-level and mechanism-level prior-art audit

**19 September 2026; baseline `3a80fec8ae7d45cc6906c20cb9c537c3ac2c9698`.** This is a retrospective literature audit and mathematical comparison, not a preregistered physical experiment. Source IDs resolve to the links below and to the inspection records in [sources.json](sources.json). The separate [claim register](claims.json) controls the strength of the conclusions.

## 1. Scope and standard of evidence

We distinguish five outcomes: an exact mathematical identity; a known general architecture adapted with different ingredients; equality only in a stated limit; a physical analogy; and a candidate contribution whose originality remains unverified. Matching an equation is not enough to identify two full theories. Variables, units, source and force couplings, initial histories, boundary conditions, approximation regimes and observables must match as well.

The previous conversation correctly identified an exact interpolation-function antecedent and close memory-model architectures. This audit preserves those results and adds specific prior work on geometric instability, pressure-supported populations, phase-adaptive orbital stability, and a 2026 collective-cluster experiment. A finite search cannot prove that no equivalent construction exists.

The reviewed code anchors are [`fields.py`](../../research_work/results/path-memory/fields.py), [`formation.py`](../../research_work/results/path-memory/formation.py), [`reciprocal.py`](../../research_work/results/path-memory/reciprocal.py), [`ring_modes.py`](../../research_work/results/path-memory/ring_modes.py) and [`warm_modes.py`](../../research_work/results/path-memory/warm_modes.py). These links are relative to this document; the pinned baseline above defines the version being audited. Numerical outcomes are not revalidated here.

## 2. Exact antecedent: the PM interpolation function

**Provenance: exact published function, rederived and numerically implemented in this project.** Famaey and Binney print the following Bekenstein-toy function as equation 5 [FB2005]:

```text
mu_B(x) = [sqrt(1 + 4*x) - 1]/[sqrt(1 + 4*x) + 1]
```

`fields.py` uses:

```text
mu_PM(x) = 4*x/[1 + sqrt(1 + 4*x)]^2
```

**Established algebra:** multiply the numerator and denominator of the first expression by `sqrt(1+4*x)+1`. Since `(sqrt(1+4*x)-1)*(sqrt(1+4*x)+1)=4*x`, the two forms are identical. Rationalizing a known expression can improve numerical evaluation; it does not create a new constitutive law.

For the algebraic spherical relation, define `x=g/a_star`, `y=g_N/a_star`. The inverse is:

```text
y = x*mu_PM(x)
x = y + sqrt(y)
g = g_N + sqrt(a_star*g_N)
```

This is the same positive-acceleration relation used by PM-1 after its force-equivalent-mass substitution. Fitting its coefficient is a calibration, not independent discovery of an acceleration scale.

**Critical scope:** [FB2005] attributes the formula to a particular tractable Bekenstein construction and states spherical/acceleration-regime qualifications. The source's relativistic theory contains tensor, vector and scalar structure [B2004]. An algebraic identity does not equate RUT to TeVeS, establish the same lensing, or justify a universal disk relation. The exact antecedent is nevertheless sufficient to withdraw novelty for the function and its algebraic inverse.

## 3. Nonlinear Poisson, auxiliary fields and the danger of spherical equivalence

**Provenance: established field-equation architecture.** Bekenstein and Milgrom [BM1984, section II, equations 2b-3] derive a variational nonlinear-Poisson model:

```text
div[mu(|grad Phi|/a0)*grad Phi] = 4*pi*G*rho
mu(x) = F'(x^2)
```

Thus Completion I's use of this architecture, and computing a compatible field functional from a chosen `mu`, are not new principles. The small-argument evaluation and numerical verification can be useful implementation work.

**Derived here from a known function; no originality claim:** put `t=(sqrt(1+4*x)-1)/2`, so `x=t^2+t` and `mu=t/(1+t)`. Integration gives the project's `F=t^4+(2/3)*t^3`, up to an irrelevant additive constant. The requirement `d[F(x^2)]/dx=2*x*mu(x)` fixes it; agreement does not supply an independent physical prediction.

QUMOND is another established shared-field construction [M2010]:

```text
laplacian Phi = div[nu(|grad Phi_N|/a0)*grad Phi_N]
laplacian Phi_N = 4*pi*G*rho
```

AQUAL and QUMOND can share a chosen spherical algebraic relation and still differ for nonspherical sources. Likewise, Completion II's Newtonian field plus separate nonlinear auxiliary field is not proved identical to either merely because it yields the same spherical total acceleration. [ZF2006] is a further direct comparison for interpolation/free-function choices. Preserve all geometry and field-coupling distinctions.

## 4. The closest particle-field architecture predates RUT

The walking-droplet experiment that inspired this project demonstrates motion writing a persistent wavefield that changes subsequent motion [F2010]. The continuous trajectory-equation construction was developed explicitly in [O2013]. These are externally driven hydrodynamic systems, not measurements of gravitational memory.

**Structural antecedent, not exact physical equivalence:** [TDR2020, equations 1-2] couples inertial particles with drag to a relaxing field continuously sourced by the particles. Schematically:

```text
m*Xddot_i + D*Xdot_i = -m*g*grad(h)(X_i)
h_t + h/T_memory = source_amplitude*sum_i H(x; X_i)
```

The one-stage RUT architecture is:

```text
Xddot_i = -grad(Phi_ext)(X_i) + grad(C)(X_i)
C_t + C/tau_keep = sum_i q_i*K(x-X_i)
```

Renaming `C=-g*h` exposes the common source-relaxation-gradient structure, but it also changes how signs and source kernels must be mapped. The published driven wave kernel is not simply our everywhere attractive Gaussian. Their viscous particle drag and annular constraint cannot be silently deleted in an equivalence claim.

The genuinely important additional overlap is conceptual: [TDR2020] already separates geometric collective instability from memory-driven oscillatory instability and develops nonlinear analysis. Therefore RUT's discovery that a cold ring is unstable even without lag is an important internal control, but the geometric-versus-memory distinction itself is not new. Their experiment [TCB2020] and amplitude-equation extension [TDR2021] should be compared before claiming new collective mode selection or saturation.

[O2017] also finds that modifying the droplet's impact-phase response can matter to orbital stability. This is a useful modeling precedent: response dynamics can change an orbit's stability. It is not the same constitutive modification as our two-pole maturation rule.

## 5. New 2026 comparator: collective stability from unstable elements

Li and Valani's **preprint**, version 3 dated 21 July 2026 [LV2026], reports persistent stationary, rotating and translating clusters sustained by a shared wavefield. The isolated droplets fail through loss of bounce synchronization and coalescence. Their model retains external dual-frequency forcing, particle drag, an oscillatory Bessel-envelope kernel and short-range repulsion.

This is strong prior art against the broad claim “we are the first to show unstable individuals form stable moving collectives through a shared memory field.” It does **not** establish that our inertial cold-ring mode or our mass-sourced two-stage model is identical. The instability being stabilized and its energy supply are different. Treat the work as an especially close comparator, not as evidence that gravity has its mechanism; peer-reviewed publication was not verified in this audit.

## 6. Chemical trails, inertia and warm populations

Chavanis [C2010, equations 54-56] gives an inertial particle/chemical-field system of the form:

```text
Xdot_i = v_i
vdot_i = -xi*v_i + grad(c)(X_i) + noise
c_t = D_c*laplacian(c) - k*c + production*sum_i delta(x-X_i)
```

**Provenance: published model, adapted comparison.** Omitting noise and friction, setting diffusion to zero, replacing point sources with Gaussian footprints and adding a central potential yields a close structural relative of one-stage RUT. Those operations change the physical model; they are not a proof that the original paper already simulated our configurations. Conversely, calling the mediator gravitational rather than chemical does not create novelty in the architecture.

[C2010] and [CS2007] also provide kinetic/continuum descriptions. [STL2009] and [K2016] are direct antecedents for self-interaction with one's own chemical or persistent trail. [L2015] shows that delayed response and anisotropy can alter collective behavior even for chemorepulsive particles; the sign of a delay's effect is not universally stabilizing.

**Established principle, model-specific quantitative question:** [CS2008] studies pressure/sound-speed thresholds for chemotactic collapse. [H1970] reports cold stellar-disk instability, suppression of small-scale disturbances by velocity dispersion, and remaining long-scale structure. Warming an attractively coupled population is therefore not a new general stabilizing principle. The location and character of a threshold for a rotating population with our specific field response may still be a contribution, provided it is derived, verified and compared rather than merely renamed.

## 7. The two-stage cascade is an exact reformulation, not a new class of memory

**Provenance: standard linear-system algebra applied to a proposed coupling.** Our equations are:

```text
tau_form*E_t = S - E
C_t = E - C/tau_keep
```

Eliminating `E` gives:

```text
tau_form*C_tt + (1 + tau_form/tau_keep)*C_t + C/tau_keep = S
H(s) = tau_keep/[(1+s*tau_keep)*(1+s*tau_form)]
```

The causal temporal impulse response, for unequal positive times, is:

```text
h(u) = tau_keep/(tau_keep-tau_form)
       * [exp(-u/tau_keep) - exp(-u/tau_form)],  u >= 0
integral h(u) du = tau_keep
```

At equal times `tau`, `h(u)=(u/tau)*exp(-u/tau)`. Initial-condition terms must be retained when comparing nonzero histories: equality of transfer functions alone compares only the zero-state linear response.

Auxiliary-variable representations of exponential memory kernels are established techniques, including generalized Langevin embeddings [SGH2011]. That does not equate their stochastic bath physics to RUT. Possible novelty would lie in the specific coupling's collective predictions, not the conversion between two first-order equations, one second-order equation and a two-pole response.

## 8. Stationary RUT is an attractive pair interaction; finite memory is not

**Derived here within the declared reciprocal model; originality unverified as a full construction.** At stationary sources with `q_i=alpha*m_i`, the RUT field satisfies:

```text
C0(x) = alpha*tau_keep*sum_j m_j*K(x-X_j)
a_extra,i = alpha*tau_keep*sum_j m_j*grad K(X_i-X_j)
```

For the positive Gaussian `K`, this is attractive. Eliminating the equilibrated field from the reciprocal energy yields:

```text
U_effective = -(alpha*tau_keep/2)*sum_i sum_j m_i*m_j*K(X_i-X_j)
```

The diagonal terms are position-independent at fixed masses and width. The factor one-half follows only after including and eliminating the field energy; summing the bare interaction alone double counts it. A stationary attractive pair model and finite-history RUT have the same static limit, not the same trajectories or preparation dependence.

“Gaussian core model” is a literature-search lead, not an established full equivalence. In particular, an interaction's sign, repulsive regularization, thermodynamic environment and dynamical rule cannot be discarded because the word Gaussian matches. The original Stillinger source remains in the unretrieved-lead register rather than supporting an equation-level claim.

## 9. Reciprocal coupling and what the energy construction does not establish

The repository's reciprocal equations are a specific proposed effective construction:

```text
W*W = K; C = W*h; q_i = alpha*m_i

tau_form*h_tt + gamma*h_t + h/tau_keep = alpha*(W*rho)
gamma = 1 + tau_form/tau_keep
```

**Provenance: derived here within a model, complete originality unverified.** The shared source and force coupling give a matter-plus-field energy with a derived damping term. This is valuable consistency, but positive-kernel factorization, variational interaction coupling and elimination of hidden bath degrees of freedom are established mathematical strategies; [Z1973] is a primary system-bath antecedent, not a demonstrated exact equivalent of RUT.

A microscopic equilibrium bath generally imposes additional assumptions on noise and dissipation. No such microscopic mapping was established here. Damping identifies an energy-transfer channel, not automatically a physical medium that absorbs energy and momentum. A fixed central potential, instantaneous spatial kernel and finite-periodic-box approximation remain separate qualifications. Neither a correct balance nor an exact temporal response proves a relativistic or causal gravity theory.

## 10. Different meanings of gravitational memory

[HM2009], [RM2014] and [M2013] establish nonlocal-gravity frameworks in which historical constitutive response contributes to gravity or wave damping; rotation-curve applications predate this project. [M2011] instead studies history dependence in modified inertia. The latter changes the body's response rather than necessarily adding a common force potential.

[QW1997] supplies a radiation-reaction/self-force antecedent. [C1991] supplies gravitational-wave memory, a lasting wave-induced effect on relative test-body positions. The earlier [F2010review] pointer is retained as an author review, not a first-discovery attribution. These uses of memory are not interchangeable with a continuously written attractive orbital footprint.

The original vector-memory suggestion also resembles ordinary scalar/vector-potential minimal coupling. **Established Lagrangian identity, not an original force structure:**

```text
L/m = v^2/2 - Phi + v dot A
a = -grad(Phi) - partial_t(A) + v cross curl(A)
```

The work-free cross product alone is not an energy-accounted dynamical model; the time derivative and the field's source/action matter. [I2024] is retained only as a contextual rotation-dependent-gravity reference. Its clock-effect details are not adopted or reverified here.

Transported memory is another distinct constitutive choice: a material derivative includes motion through the field. [AP2024] is an inspected primary example with transported/deformed viscoelastic stress. No exact mapping from its tensor law to our suggested vector state was established. The inaccessible Oldroyd original is a follow-up lead, not evidence inspected at equation level.

## 11. Ring integrals, kinetic response and numerical methods

**Established mathematical identities:** the Gaussian ring average follows directly from [DLMF, 10.32.1]:

```text
average_theta exp[-(r^2+R^2-2*r*R*cos(theta))/(2*w^2)]
    = exp[-(r^2+R^2)/(2*w^2)]*I0(r*R/w^2)
    = exp[-(r-R)^2/(2*w^2)]*I0e(r*R/w^2)
```

The broadband wave-ring construction is another application of established Bessel identities to an assumed spectrum. Equal potential weight per logarithmic wavenumber can yield the on-ring expression `C*[J0(k_min*R)^2-J0(k_max*R)^2]/(2*R)`. That algebra does not derive the spectrum, its source normalization, its cutoffs or a universal gravitational law. The squared-Bessel expression is an on-ring case, not the force of an arbitrary fixed source ring at every probe radius. Its full physical construction remains originality unverified.

The wave-ring Bessel addition formulas, gradients of a prescribed potential, epicyclic linearization, Routh-Hurwitz tests, and coordinate Jacobians are established tools. In particular, converting `v_theta` to `L=R*v_theta` cancels the planar spatial `R` in the joint phase-space measure. Correcting that sampler was a correctness repair, not a new statistical-mechanics measure.

[P2024] is a direct prior-method comparison for distribution-function response matrices and complex stellar modes. [Beyn2012] supplies the contour eigensolver. [Liu2026] is a newer contour/region-partitioning method worth evaluating for difficult spectra. Adapting these methods to a specified interaction can yield original results; their methods should not be claimed as invented here. No formula from a marginally sampled contour establishes full spectral completeness without its own numerical domain/error checks.

## 12. Novelty assessment and constructive research direction

The search found an exact antecedent for the PM interpolation and strong antecedents for the broad architecture and several proposed qualitative insights. It did not establish a publication containing the complete RUT source law, reciprocal two-timescale field, orbital populations and the same quantitative predictions. This limited absence is **not** a novelty certificate.

Candidate contributions to investigate are: a general support-versus-torque relation under matched static attraction; predictive mode selection across zero-, one- and two-stage responses; a verified quiet region for self-consistent kinetic populations; and a history-dependent observable that distinguishes the model from its static closure. Each must state what earlier work already established and precisely what additional result remains.

The most productive next step is to reproduce the closest published controls and then change one ingredient at a time. Preserve the physics constraints of the comparator. A driven wave bath cannot quietly become an isolated gravitational reservoir; a chemotactic pressure threshold cannot be transplanted as a galactic stability threshold; the PM fit cannot be credited to RUT without a derivation connecting them.

Existing numerical failures and uncertainties remain untouched. This audit neither resolves the stage-8 ring-limit gates nor explains B13's rate discrepancy. Novelty, correctness and empirical success are independent. The fictional-universe approach is compatible with borrowing known mathematics, provided the borrowing is explicit and the physical proposal earns its own predictions.

**Permitted current claim:** “We develop and test a specified collective gravitational-memory hypothesis built from established mathematical structures. Several components have exact antecedents. Originality of the complete construction and its quantitative results remains under investigation.”

[FB2005]: https://arxiv.org/html/astro-ph/0506723v2
[B2004]: https://arxiv.org/abs/astro-ph/0403694
[BM1984]: https://articles.adsabs.harvard.edu/pdf/1984ApJ...286....7B
[ZF2006]: https://arxiv.org/abs/astro-ph/0512425
[M2010]: https://arxiv.org/html/0911.5464
[F2010]: https://arxiv.org/abs/1307.6051
[O2013]: https://doi.org/10.1017/jfm.2013.581
[O2017]: https://doi.org/10.1103/PhysRevFluids.2.053601
[TDR2020]: https://arxiv.org/html/2003.02220v1
[TCB2020]: https://doi.org/10.1103/PhysRevFluids.5.083601
[TDR2021]: https://arxiv.org/abs/2010.12655
[LV2026]: https://arxiv.org/html/2604.09506v3
[C2010]: https://arxiv.org/html/0804.4425v2
[CS2007]: https://arxiv.org/abs/0706.3974
[CS2008]: https://arxiv.org/abs/0708.3163
[STL2009]: https://arxiv.org/abs/0905.1316
[K2016]: https://arxiv.org/abs/1504.06814
[L2015]: https://arxiv.org/abs/1508.04673
[SGH2011]: https://arxiv.org/abs/1011.2848v2
[Z1973]: https://doi.org/10.1007/BF01008729
[AP2024]: https://arxiv.org/abs/2401.03981
[HM2009]: https://arxiv.org/abs/0812.1059
[RM2014]: https://arxiv.org/abs/1401.4819
[M2011]: https://arxiv.org/abs/1111.1611
[M2013]: https://arxiv.org/abs/1304.1769
[C1991]: https://doi.org/10.1103/PhysRevLett.67.1486
[QW1997]: https://arxiv.org/abs/gr-qc/9610053
[F2010review]: https://arxiv.org/abs/1003.3486
[I2024]: https://arxiv.org/abs/1407.5022v5
[H1970]: https://doi.org/10.1177/003754977001500307
[P2024]: https://arxiv.org/abs/2311.10630
[Beyn2012]: https://doi.org/10.1016/j.laa.2011.03.030
[Liu2026]: https://doi.org/10.1002/nla.70072
[DLMF]: https://dlmf.nist.gov/10.32
