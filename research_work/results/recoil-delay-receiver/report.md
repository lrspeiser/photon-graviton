# What a physical delay mechanism must account for

The previous forward massless-stream model cannot create a storage delay by rearranging positive-energy photons and companions alone, under the relativistic assumptions stated below. A pre-existing massive receiver supplies a constructive kinematic alternative: it can absorb a pulse, move while carrying its energy, and release it later. The receiver's displacement accounts for the radiation delay in the center-of-energy balance.

This demonstrates an allowed energy/momentum bookkeeping route for delay, not a microscopic interaction or a new theory of time. A constant residence time still does not stretch event spacing. A changing residence schedule can stretch events, but our calculation does not derive that schedule or connect it to the photon-frequency change.

## 1. The restriction on a photon-only input

**Known positive-energy relativistic kinematics, not a new project theorem.** Work in c=1 units. An initially isolated, exactly forward beam has total energy E and momentum P_x=E. For constituents with future-directed causal four-momenta,

    E_i >= |p_i| >= p_x,i
    sum_i (E_i - p_x,i) = E_total - P_x,total = 0.

Each term is nonnegative, so every occupied outgoing channel must satisfy E_i=p_x,i. Thus every such constituent is massless and exactly forward. A timelike storage carrier, an angled outgoing beam, or a stationary positive-energy deposit would require a positive contribution to E-P_x, absent from this input. Three linear-program controls maximize energy in slower states and obtain zero; a control without a forward null state is infeasible. The analytic inequality, rather than the finite velocity grid, supplies the general statement under these assumptions.

This assumes an isolated flat-spacetime system with standard conserved four-momentum and positive causal constituent energies. It is not a prohibition for all interacting fields or modified gravities. A pre-existing matter/field background changes the initial four-momentum; additional stresses, nonstandard energy conditions or a different spacetime law require their own explicit accounting. Calling the rate a time field does not by itself specify that accounting.

The energy–momentum and Lorentz-transformation framework is standard; see [Feynman's relativistic dynamics discussion](https://www.feynmanlectures.caltech.edu/I_16.html). Related established pulse/body center-of-energy reasoning is studied in the [relativistic Balazs thought-experiment analysis](https://www.sciencedirect.com/science/article/pii/S000349162030227X). The discrete absorption/release calculation below is a conditional diagnostic built with that mathematics, not a claim of a unique fundamental law.

## 2. An exact massive-receiver cycle

**Conditional derivation from known relativistic collision kinematics.** Let M denote the receiver's initial rest energy, not its mass in kilograms. It starts at rest. An incoming forward photon has energy E. After absorption, the excited receiver has total energy M+E, momentum E, rest energy mu, speed beta and Lorentz factor gamma:

    mu = sqrt(M^2 + 2 M E)
    beta = E/(M+E)
    gamma = (M+E)/mu.

Assume it stores the excitation for proper time tau. Its laboratory storage duration and displacement are

    T = gamma tau
    delta_x = beta T.

**Hypothetical emission channel.** At release, stipulate one forward photon with energy qE and forward massless companions carrying (1-q)E. The receiver returns to rest with rest energy M. The outgoing radiation carries total energy E and momentum E, so the receiver's absorbed momentum is returned to radiation. The value of q, the actual allowed transitions and their amplitudes are not derived here. Assigning photon frequency through E=h_P nu gives a spectral factor 1/q, but does not specify a coherent wave transformation.

For a detector ahead of the release point, the radiation arrives later than a freely propagating pulse by

    delay = (1-beta) T = (M/mu) tau
    M delta_x = E delay.                     [c=1]

With M still denoting rest energy, the last identity in ordinary units is M delta_x = E c delay. The receiver's final momentum and excitation return to their initial values, but its position does not. Resetting that displacement would require accounting for the corresponding additional process.

For M=E=1, tau=1 and q=0.8, the laboratory residence is 1.154701, displacement and arrival delay are each 0.577350, and the outputs carry 0.8 photon energy plus 0.2 companion energy. The excited internal energy is 0.732051 and its kinetic energy is 0.267949. The sum accounts for the entire incident energy.

Thirty-six cycles vary M/E, q and proper residence time. All displayed four-momentum balances pass, including 180 boosted-frame checks. Maximum boost residual is 3.7e-12 in the normalized units; the lab center-of-energy error is below 7.2e-15. Both the transient stored energy and recoil are included. These are free-particle/collision kinematic checks; a real controller, absorption process and emission interaction could require additional modeled sectors.

## 3. Delay still differs from event stretching

**Additional hypothetical clock schedule, not a derived cause:** let the receiver's proper residence time depend on its own age at absorption,

    w_n = 0.1 + kappa tau_receiver,n.

The receiver's clock is assumed to measure ordinary proper time. It is never reset; it advances during storage and the intervening rest intervals. We send four equal-energy pulses at source times 0,12,24,36, with a fixed source/detector separation and no overlapping storage episodes. The example's lengths, times and energies are normalized, not fitted astronomical quantities.

Writing a=M/sqrt(M^2+2ME), exact collision/clock accounting gives

    t_arr,n = t_emit,n + D + a w_n
    tau_receiver,n+1 - tau_receiver,n = Delta_t_emit + (1-a) w_n
    S_event,n = 1 + a kappa [1 + (1-a) w_n/Delta_t_emit].

These are conditional consequences of the stipulated residence schedule and standard kinematics. They are not unique laws of time. The receiver's displacement is included when computing the next interception and detector arrival.

| M/E | kappa | First event stretch | Last event stretch | Assigned spectral factor for q=0.8 |
|---:|---:|---:|---:|---:|
| 1 | 0 | 1.000000 | 1.000000 | 1.250000 |
| 1 | 0.2 | 1.124011 | 1.145862 | 1.250000 |
| 10 | 0.2 | 1.185358 | 1.191874 | 1.250000 |
| 10,000 | 0.2 | 1.199984 | 1.199992 | 1.250000 |

All twelve pulse-sequence cases agree with these recurrence identities. Constant residence gives a constant delay and unchanged event intervals. Changing residence stretches the intervals, but the factors depend on receiver state and prior pulses. Repeating the schedules with q=1 changes no arrival times, despite removing the assigned redshift. Thus the energy split and residence dynamics are independent in this construction.

Choosing q=1/S_event after computing the intervals would impose the desired correspondence, not derive it. A coherent whole-signal remapping also needs a phase equation that yields the frequency/time relation; our discrete collision prescription has none. The earlier stationary-versus-coherent quantum-receiver restrictions remain applicable and unresolved by this kinematic example.

## 4. What this changes in the research direction

The calculation makes the missing physical role concrete: a delay-capable receiver can carry excitation, momentum, displacement and a changing state. A proposed time/companion field could potentially fill that role, but its physical stress, initial state, propagation and coupling must be specified. This result does not establish a population of massive receivers in voids or adopt permanent void deposits.

The next candidate must join the receiving-state dynamics to the light's coherent phase and derive the relation between the energy fraction q and event stretching S from the same interaction. It must include controller/background energy and momentum, receiver preparation and repeated use. The present example demonstrates kinematic permission for delay, while leaving that central dynamical derivation open.

No astronomical data were refitted and no holdout outcomes were opened. The spectral/time tests, stable stellar-population likelihood, capture and supported well deposits, lensing, and deferred total photon-supply budget remain requirements of the full goal. Run `python research_work/results/recoil-delay-receiver/run.py` from the repository to reproduce the controls and numerical examples.
