# Candidate interactions and the missing derivations

This is preliminary T03 work. It supplies explicit comparison actions where possible and marks missing actions honestly. No candidate is adopted. The subsequent [universe contract](../../../research_plan/universe-contract.md) excludes assumed dark matter, expansion and Big-Bang premises; candidate C is inactive as an independent dark-matter substitute. Earlier open-branch statements below are subject to that restriction. The operational matter-frame choice in T02 remains open, and no complete microscopic photon–companion–deposit theory has yet been established.

## 1. Conservation must include momentum and interaction energy

An energy ledger is one projection of a fuller balance law. For a diffeomorphism-invariant action, the on-shell total stress tensor obeys:

```
nabla_mu T_total^(mu nu) = 0
T_total = T_photon + T_companion + T_deposit + T_matter
          + T_driver + T_receiving + T_interaction
```

This is a conditional statement about an action with its dynamical equations satisfied. A prescribed time-dependent background is not automatically a closed conserved system. Interaction energy cannot in general be uniquely assigned to one component; whichever convention is used must sum to the same total.

One may introduce exchanges Q_i^nu through nabla_mu T_i^(mu nu)=Q_i^nu, with sum Q_i^nu=0 after including the interaction/driver accounts. The earlier rates h u_gamma and Gamma u_c specify only energy exchanges in a chosen frame. Their spatial momentum exchanges remain unspecified.

For a fixed volume in a flat stationary frame, integrating the energy equation produces the earlier source-minus-boundary-flux ledger. For an isotropic comoving component in a homogeneous metric, the energy equation instead has the form dot rho+3H(rho+p)=Q, with c=1. Equivalently, d(rho V)/dt=QV-p dot V. The work term has an energy destination in the complete system. A homogeneous loss term cannot simply be added or removed without stating which geometry, clock and receiving sector define it.

In curved spacetime local covariant conservation does not by itself supply a globally conserved ordinary energy integral. A stationary spacetime with a timelike Killing vector supplies a corresponding conserved energy current when all relevant sources are included. A time-dependent fictional background requires its own full accounting. Thus the earlier fixed-volume equality is valid within its stated scope, not a proof of every proposed global history.

## 2. Static conversion and continuous photon redshift are different mechanisms

Consider a specified lossless two-mode propagation model at one frequency:

```
i d/ds (A_gamma, A_c)^T = H (A_gamma, A_c)^T
H = [[delta/2, g], [g, -delta/2]]
```

The real constants delta and g are inverse lengths. Hermiticity gives d(|A_gamma|²+|A_c|²)/ds=0. Starting in the photon mode:

```
P_gamma_to_c = g²/[g²+(delta/2)²]
               × sin²(s sqrt[g²+(delta/2)²])
```

This changes mode occupancy. In this stationary linear model each Fourier frequency evolves independently, so it does not shift that frequency. A frequency-dependent conversion probability could distort or selectively attenuate a spectrum, which must not be confused with a derived coherent redshift of all its lines. A time-dependent driver or an inelastic interaction can change frequency, but requires its own energy/momentum balance and spectral prediction. The numerical test verifies the stipulated matrix and its analytic probability; it does not derive g or delta from a field action.

Photon–graviton conversion in magnetic environments is a real comparison mechanism, with important background/gauge qualifications. It is not evidence that the fictional redshift-and-deposition chain already works. See [Hwang and Noh, On graviton-photon conversions in magnetic environments](https://arxiv.org/abs/2310.04150).

## 3. Capture can conserve energy while producing recoil

Take a massless incident companion of energy E and momentum E/c absorbed by an initially stationary freely moving target of rest mass M. Assume a single final excited target, with no other outgoing particles. Four-momentum conservation requires:

```
M_final² c⁴ = (M c²+E)² - E²
M_final = M sqrt(1+2x),      x = E/(M c²)
final speed/c = x/(1+x)
stored rest-energy gain/E = 2/[sqrt(1+2x)+1]
recoil kinetic energy/E = x/[1+x+sqrt(1+2x)]
```

The two energy fractions sum to one. At x=10^-6, 99.99995% of incident energy becomes target rest-energy increase; at x=1, about 73.2% does and 26.8% becomes recoil. Both fractions gravitate through the full stress tensor; the distinction matters when assuming a cold stationary deposit. A massive collective absorber can make recoil small. Equal opposed packets can have zero net momentum. An anchored absorber transfers momentum to its support, whose energy and stresses must be included.

This is not a requirement that capture wastes energy. It shows how to retain the energy while explicitly describing momentum and the deposited state's motion. If the target is ordinary matter, its increased mass already contains the captured energy: adding a separate deposit of E/c² on top of that would count it twice.

If the final target has a specified mass M_e, absorption is resonant in this isolated two-to-one example:

```
E_absorb = (M_e²-M²)c²/(2M)
```

The inverse decay into the original target and a massless companion is also kinematically allowed:

```
E_emit_in_excited_rest_frame = (M_e²-M²)c²/(2M_e) > 0
```

Kinematic allowance does not determine the lifetime. A matrix element, selection rules, environment and available final states are still required. An ordinary Hermitian transition interaction supplies the reverse vertex as well, so permanent storage cannot be inferred from writing down absorption alone. Candidate long-lived collective states, suppressed transitions, changes of phase and escape barriers remain open. Their storage times must be calculated, and “permanent” cannot be claimed from a merely long finite lifetime.

## 4. Bounded candidate action register

Use signature (-,+,+,+), c=hbar=1, and a physical matter metric g unless explicitly stated otherwise. A “graviton” below means an ordinary propagating tensor mode only in candidate A. Calling a new scalar or stored state a graviton does not establish that identity.

### A. Ordinary tensor radiation in Einstein–Maxwell theory

```
S_A = integral sqrt(-g) [M_P² R/2 - F_mu_nu F^mu_nu/4] d4x
      + S_matter[g,A,Psi]
```

The weak-field propagating content is two tensor polarizations and two photon polarizations, plus matter. Diffeomorphism and electromagnetic gauge invariance constrain the interactions. Expanding g=background+h/M_P gives a tensor coupling to the electromagnetic stress tensor, up to normalization/sign conventions set by the metric variation. External electromagnetic backgrounds can couple photon and tensor perturbations; the currents and stresses maintaining such backgrounds belong to the energy budget.

This is a low-energy comparison action; perturbative use requires energies and background curvatures within its effective regime, not a quantum-gravity completion. No independent D field exists in it. A radiation configuration or matter absorber must represent any deposit. A stored gravitational wave is not automatically cold dust.

Do not claim that massless radiation can never be gravitationally confined: [Andréasson, Fajman and Thaller](https://arxiv.org/abs/1511.01290) construct highly relativistic self-gravitating photon shells in the massless Einstein–Vlasov system. That existence result does not establish a dilute galaxy halo or this model's formation/stability requirements.

Open tasks: derive source conversion in the relevant backgrounds, test frequency and clock consequences, and construct a storage state with its full stress tensor. The homogeneous transparent version does not supply the intended nonexpanding redshift.

### B. A new scalar companion and explicit excited absorbers

A bounded illustrative effective action is:

```
L_B = M_P² R/2 - (1+g_gamma chi/Lambda) F²/4
      - (partial chi)²/2 - m_chi² chi²/2 - lambda_chi chi⁴/4
      + sum[a=g,e] bar(Psi_a)(i gamma^mu D_mu-M_a)Psi_a
      - [y chi bar(Psi_e)Psi_g + Hermitian conjugate]
      + L_other_matter
```

Here chi is a real scalar, not the ordinary tensor graviton; Psi_g and Psi_e are neutral massive Dirac absorber states with M_e>M_g. D_mu includes the spin connection. Their shared fermion number is conserved. Electromagnetic gauge invariance is retained. The displayed electromagnetic coupling is dimension five, and its cutoff Lambda must be specified before quantitative predictions. Require a positive electromagnetic kinetic coefficient and perturbative amplitudes below the cutoff; m_chi²>=0 and lambda_chi>=0 give a bounded uncoupled scalar potential. These limited conditions are not a complete interacting-background stability or radiative-stability proof.

Setting m_chi=0 is a candidate massless companion choice, not a protected prediction: radiative corrections and symmetries must be checked. The absorber transition demonstrates where recoil and inverse emission enter; it supplies no permanently stable excited state by itself. The chi F² vertex supplies interactions but does not derive a universal photon redshift or its rate. A background chi affects electromagnetic normalization and matter response, which must be calculated using the same action.

In particular, multiplying Maxwell's F² by a slowly varying positive scalar Z changes its principal equation to Z times the usual Maxwell principal operator. In geometric optics its leading characteristic remains the g-null cone; it is not, by itself, the earlier omega=ck/n(t) law. Gradients/background evolution, amplitudes and lower-order effects require separate treatment. This candidate has not escaped the earlier clock issue just by adding a field.

The candidate has a finite set of displayed couplings, not adjustable functions of radius. Its absorber density, level spectrum and source history are further physical inputs requiring independent specification. A complete effective expansion would also track other operators generated by its interactions.

### C. A stable massive stored-particle sector

An alternative bounded action is:

```
L_C = M_P² R/2 - [1+kappa |X|²/Lambda²] F²/4
      - |partial X|² - m_X² |X|² - lambda_X |X|⁴/2
      + L_matter
```

X is a neutral complex scalar with a global U(1) symmetry; it has two real field degrees of freedom and is distinct from the tensor sector. Work with m_X²>0, lambda_X>=0, a positive electromagnetic kinetic coefficient, and energies/field amplitudes within the chosen effective regime. These are vacuum/kinetic prerequisites, not a proof that a halo background is stable.

The interaction permits photon/scalar pair processes when their invariant energies meet the thresholds. It does not turn a single freely propagating photon into an arbitrary massive pair, nor supply the desired continuous redshift law. The conserved charge can protect an isolated lightest charged particle against decay into neutral states; particle–antiparticle annihilation remains possible. Thus charge conservation alone does not make an optically produced neutral collection permanent.

Massive slow particles offer a candidate storage identity. Production, cooling, capture, annihilation, spatial support and the source energy all remain to be calculated. In its present form this is a candidate radiation-fed extra matter sector, rather than a derivation that ordinary gravitons accumulate. It is kept as an explicit alternative, not silently substituted for the user's concept.

### D. A frequency-changing driver with a distinct companion law

No complete action is supplied yet for the intended combination of: nonexpanding material geometry, observable photon frequency loss, no secondary companion loss, protected local clocks, and a permanently stored halo state. Its register status is **incomplete action**, not a physically excluded concept.

The recovered preferred-frame electromagnetic/atomic toy and universal conformal action are useful comparison calculations. The former gives explicit competing atomic/cavity scalings; the latter gives 1+z=A_o/A_e in its homogeneous transparent matter frame and therefore no homogeneous redshift when that material scale is constant. Their existing derivations are retained rather than reinterpreted as a proof against every nonmetric, inelastic or environment-dependent theory.

For this candidate, the next action must specify the dynamical driver, the matter frame, finite coupling parameters, the mechanism of frequency transfer, and the companion's distinct propagation law. A prescribed n(t) and a separately imposed lambda_c=0 are not enough. Any new vector/preferred-frame or nonlocal operator requires its own degrees-of-freedom and stability analysis. Nothing has been inserted here solely to force the desired result.

## 5. Deliverables that would turn a candidate into a testable mechanism

1. Vary one chosen action for the background, matter, light, companion and stored state. Derive stress tensors and exchanged four-momentum from that action.
2. Calculate conversion amplitudes or a controlled classical limit, including reverse channels. Identify whether the result transfers whole quanta, gradually shifts photon energies, changes arrival rates, or combines these effects.
3. Derive a spectral/phase-space transport kernel. Its energy and momentum moments must reproduce the ledger without independently adjustable loss terms. Predict line shifts, broadening, image blur and transient arrival times.
4. Derive capture cross-sections, recoil, excitation/cooling and storage lifetime. A stationary cold deposit requires dynamical support and an equation of state, not just a scalar density history.
5. Evaluate the same interaction for atoms and cavities. If an environmental exemption is used, derive it from the interaction and boundary conditions.
6. Derive gravitational motion and lensing with all energy and stresses included. The fitted response multiplier eta remains a placeholder until these equations exist.
7. Recompute supply and halo evolution using those rates; preserve the archived benchmarks as comparisons. Only then can new observational fits test the unified candidate.

## Verification and provenance

`microphysics-checks.json` records five check groups: capture energy partition over 151 ratios, four-momentum mass shell, inverse-channel kinematics, opposed-packet momentum cancellation, and three static-mixing matrix cases. These are mathematical checks of stated examples, not validation of any action above. `check_microphysics.py` reproduces them.

Recovered sources inspected: `minimal_clock_interaction/derivation.md`, `conformal_action_derivation/derivation.md`, `interaction_stability_attempt/checks.json`. External primary sources were checked on 9 September 2026 and support only their explicitly attributed comparison claims. The scalar/absorber actions are proposed illustrative candidates in this work; they are not presented as established solutions or claims of novelty.

## Later capture/storage calculation

The [reversible receiver calculation](../capture-storage/derivation.md) develops the bandwidth/lifetime tradeoff, a three-state storage model, energy-release ledger, finite capacity and reverse bath transitions. It remains a conditional effective model; no actual material or ordinary-graviton capture strength is identified.
