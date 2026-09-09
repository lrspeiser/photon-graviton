# Unscreened time field: light, atoms, energy and radiation temperature

Research addendum — 8 September 2026. Fictional, spatially nonexpanding universe; real measurements used as benchmarks. These are restricted-model derivations, a numerical conservation experiment, and a reused FIRAS spectrum diagnostic. They do not establish a new physical theory or constitute a new blind astronomical test.

## What this step establishes

An explicit light–atom model can produce an observable redshift without screening. However, the simplest versions face a measurable tradeoff between atomic physics and the lengths of material rulers. A conservative field can receive the energy lost by radiation. Uniform, adiabatic frequency reduction can preserve an existing Planck distribution while reducing its temperature. None of these results independently determines the original coefficient or creates a thermal background from arbitrary light.

## 1. Specify what the field changes

Let a dimensionless field χ set an optical factor n = exp(χ), normalized to one now. In fixed background spatial coordinates, take the low-energy electromagnetic Lagrangian density to be

\[
\mathcal L_{EM}=\frac{\epsilon(\chi)}2 E^2-\frac{B^2}{2\mu(\chi)},\qquad
\epsilon=\epsilon_0 n^a,\quad\mu=\mu_0 n^{2-a}.
\]

Consequently light propagates with cγ = c₀/n. Positivity of ε and μ gives positive electromagnetic energy in this nondispersive approximation. Homogeneity conserves spatial wave number; slow temporal evolution gives ω = c₀|k|/n. This is frequency reduction at fixed coordinate wavelength, rather than literal elongation of the spatial wave in these coordinates. Astronomical “stretching” refers here to longer periods relative to a detector.

Use fixed charge e and ℏ, and an electron inertial mass mₑ = mₑ₀ nᵈ. The hydrogen approximation, with a heavy nucleus, gives

\[
H_{atom}=\frac{p^2}{2m_e(\chi)}-\frac{e^2}{4\pi\epsilon(\chi)r},\qquad
\nu_{atom}\propto n^{d-2a},\qquad a_B\propto n^{a-d}.
\]

Define p = 2a − d, so νatom ∝ n⁻ᵖ. This p is the atomic-clock exponent, not a gravitational-wave coupling. Reduced masses require corresponding nuclear mass laws; the example with changing masses assumes a constant proton/electron mass ratio.

**Relativistic completion assumption:** for the following fine-structure comparison, assume matter's relativistic limiting speed follows cγ, allowing the usual spectroscopic parameter α = e²/(4πεℏcγ) ∝ n^(1−a). This is an extra physical assumption. If matter and light have different limiting speeds, a separate calculation of relativistic atomic levels and Lorentz-violation observables is required; the α constraint below cannot simply be transplanted to that branch. The nonrelativistic hydrogen model alone does not settle that question.

This specification goes beyond calling something “time”: it says how traveling light and material standards respond. A common homogeneous lapse affecting everything identically can instead be a change of time coordinate and need not generate an observable redshift.

## 2. Derive the observed shift and duration stretching

Write S = nᵣ/nₑ for reception/emission optical factors. A photon emitted on the atomic line arrives with frequency νatom,e/S. Comparing it to the same local atomic transition gives

\[
\boxed{1+z_{obs}=S^{1-p}}.
\]

For a source and receiver fixed in background coordinates, differentiating the ray integral gives dtᵣ/dtₑ = S. Atomic clocks count dτ ∝ n⁻ᵖdt, so

\[
\boxed{\frac{d\tau_r}{d\tau_e}=S^{1-p}}.
\]

Thus frequency redshift and infinitesimal pulse-duration stretching agree in this model. This does not automatically derive every supernova evolution timescale: source dynamics and additional couplings must also be calculated. Peculiar velocities, gravitational potentials and finite-duration field evolution require additional terms.

| Candidate, for increasing n | Atomic-clock p | Observable 1+z | Advantage | Limitation |
|---|---:|---:|---|---|
| a=0, d=0: fixed Coulomb response and mass | 0 | S | Redshift and leading hydrogen clock/length remain intact | Under the shared-speed completion, α changes at the full field rate |
| a=1, d=0: equal electric/magnetic response, fixed mass | 2 | 1/S | α remains constant | Gives a measured blueshift for this field direction |
| a=1, d=2: equal responses, compensating masses | 0 | S | Retains redshift and constant α | Atomic lengths shrink as 1/n; mass generation and gravity remain unspecified |
| a=1, d=1: equal responses, fixed atomic length | 1 | 1 | Constant α and atomic size | The common shift cancels against the reference clock |

All candidates remain in the registry. A reversed field direction or additional interactions defines another candidate and requires consistent refitting.

## 3. Compare the coupling to a real clock measurement

Using the supplied κ = 0.000077315 per million light-years, the p=0 local rate is γ = κc₀ = 7.7315×10⁻¹¹ per year. Lange et al. report α̇/α = (1.0 ± 1.1)×10⁻¹⁸ per year [1]. Under the completion and standard sensitivity assumptions above,

\[
\frac{\dot\alpha}{\alpha}=(1-a)\gamma,
\qquad -1.50\times10^{-8}<1-a<4.08\times10^{-8}
\]

using a Gaussian 1.96σ interval from that published summary. The a=0 branch predicts 7.73×10⁻¹¹ per year and is strongly mismatched under these assumptions. The a=1 branch removes this particular variation by construction; that is not proof of compatibility with all experiments. The comparison uses published inferred parameters, not a reanalysis of clock time series.

For p≠0, maintaining the observed nearby slope requires (1−p)χ̇ today = κc₀ when today's clock is normalized to coordinate time. Accordingly the displayed numerical interval is conditional on the p=0 rate, not a general bound on every branch.

There is a compact tradeoff within this family. The atomic-length exponent ℓ = a−d satisfies

\[
\ell+(1-p)=1-a.
\]

Keeping α constant requires a=1. Retaining a nonzero redshift exponent then requires changing atomic length. For the p=0 case a_B ∝ 1/n: a fixed coordinate separation contains an increasing number of atomic rulers. Its present fractional increase is γ. This can resemble expansion operationally even though the background coordinates have fixed volume. Deciding what “constant-size universe” means therefore requires a physical metric and material measurement prescription. This is a limitation of this restricted family, not a theorem excluding every unscreened time theory.

## 4. Derive distance laws instead of imposing one

With reception at t=0, n(0)=1 and fixed separation R,

\[
R=\int_{t_e}^{0}\frac{c_0}{n(t)}dt.
\]

For a rolling field χ=γt, n=e^(γt), giving

\[
1+z_{obs}=\left(1+\frac{\gamma R}{c_0}\right)^{1-p}.
\]

For n=1+γt, the integral instead yields

\[
1+z_{obs}=\exp\left[(1-p)\frac{\gamma R}{c_0}\right].
\]

The second history reproduces the original exponential for p=0, but selecting that history is an assumption, not a first-principles derivation. It reaches n=0 at t=−1/γ, about 12.9 billion coordinate years ago. That marks a failure/boundary of this simple history, not an established universe age or a Big Bang. A continuation requires new dynamics. The simplest free rolling field produces the first law, and the already studied nearby data barely distinguish these possibilities.

An alternative derivation of the exponential assumes continuous, stationary multiplicative stretching along path segments: S(R₁+R₂)=S(R₁)S(R₂). That fixes functional form but neither its microscopic cause nor its numerical rate. An evolving field need not satisfy the stationarity premise across observation epochs.

## 5. Where the radiation energy goes

Give χ a positive kinetic term Kχ̇²/2 and potential V(χ). In a homogeneous fixed-volume box, adiabatically conserved photon occupations imply radiation energy Uγ = U₀e⁻χ. The reduced Hamiltonian is

\[
H=\frac{P_\chi^2}{2K}+V(\chi)+U_0e^{-\chi}+H_{matter}(\chi).
\]

There is no explicit time dependence. Its field equation contains the radiation source:

\[
K\ddot\chi=-V'(\chi)+U_\gamma-\partial_\chi H_{matter}.
\]

Radiation loses energy at −χ̇Uγ; the field/matter system receives the compensating amount. This is a precise version of the proposed energy exchange. It does not establish that time itself consumes energy. χ is a hypothesized physical degree of freedom.

The included code integrates the unscreened radiation-only case, K=1 and V=0, starting χ̇=0.2 in dimensionless units. Two initial radiation energies give:

| Initial radiation energy | Radiation energy change | Field energy change | Maximum relative total-energy error |
|---:|---:|---:|---:|
| 0.001 | −0.000635484 | +0.000635484 | 1.73×10⁻¹³ |
| 0.1 | −0.083783845 | +0.083783845 | 1.26×10⁻¹² |

These are synthetic conservation checks, not observed energy transfers. They also show that radiation backreaction accelerates this field: a prescribed constant χ̇ is not an exact solution when radiation is appreciable. A complete theory must specify matter rest energy, mass generation, gravity, the potential and initial conditions. The coefficient κ is still an empirical input.

## 6. Does stretching light change its temperature?

A single photon has energy E=hν; it does not have a unique thermodynamic temperature. For a thermal radiation bath the occupation is

\[
f(\nu)=\frac1{\exp[h\nu/(k_BT)]-1}.
\]

If every mode's frequency is divided by S and its occupation is preserved, the same distribution is Planckian at T/S. This requires a frequency-independent, sufficiently slow change with no absorption, emission, scattering distortion or particle production. Under the homogeneous model,

\[
\boxed{T_r=T_e/S},\qquad
\boxed{\frac{(k_BT/h\nu_{atom})_r}{(k_BT/h\nu_{atom})_e}=S^{p-1}=\frac1{1+z_{obs}}}.
\]

The second relation specifies what a temperature comparison against atomic energy standards means when those standards evolve. For p=0 it reduces to the familiar numerical Tᵣ=Tₑ/(1+z). Source temperature is not physically lowered merely because its emitted radiation arrives redshifted.

The code checks the occupation transformation over 1,000 dimensionless frequencies for S=2, 10 and 1100.7; maximum relative mismatch is approximately 1.4×10⁻¹⁴, numerical roundoff. This verifies the algebra, not nature.

Energy density needs additional care. With cγ=c₀/n and unchanged mode count in fixed coordinate volume,

\[
u_\gamma=\frac{\pi^2k_B^4}{15\hbar^3c_0^3}n^3T^4.
\]

Since T∝1/n, uγ∝1/n and photon number density stays constant in that volume. Substituting the ordinary fixed-c coefficient into u∝T⁴ would give an incorrect energy budget. At today's n=1 the usual formula is recovered. Evolving atomic rulers change the operational volume comparison as well. Other scalar coupling models can distort the distribution, so Planck preservation is a result for this specified adiabatic model, not for all time-based proposals [2].

## 7. What real CMB data establish here

The code reuses the 43-channel public FIRAS residual table [3], fitting temperature relative to its supplied 2.725 K reference and the supplied Galaxy template. The diagonal statistic is χ²=45.0187 for 43 channels and two fitted parameters. This is a limited spectrum diagnostic, not the complete covariance/calibration likelihood or a new absolute temperature determination. A published combined temperature estimate is 2.72548 ± 0.00057 K [4].

Any pair (source temperature, stretching) with Tsource/S equal to the fitted final temperature produces the same spectrum under our assumptions. For example, choosing source temperatures 10, 300, 3000 or 6000 K and adjusting R in Tsource exp(−κR) reproduces the same fitted final temperature to floating-point precision. **A good Planck fit therefore cannot determine the stretching history or identify the mechanism.** This degeneracy is the main data-facing result of this step.

As an arbitrary example, cooling an existing 3000 K thermal bath to 2.72548 K requires S=1100.72. The exponential law assigns about 90.6 billion light-years; the simple rolling-field law assigns about 14.2 trillion light-years, both for p=0 and fixed background geometric distance. These are large extrapolations of a nearby fit, not measured CMB distances or inferred ages. The 3000 K choice makes no assumption of a Big Bang or a known source.

Stretching does not explain why an initial bath was thermal, choose its temperature, erase arbitrary mixtures of different temperatures, or create the observed angular peaks, polarization or large-angle alignment. A spherical or constant-volume topology supplies no automatic thermalization mechanism. The earlier spectrum-mixing and angular-spectrum diagnostics remain relevant; no new alignment test was performed here.

## 8. Goals, desired outcomes and present status

| Goal | Outcome we would hope to obtain | What this step finds |
|---|---|---|
| Preserve observable redshift without screening | Light and atomic responses derive a nonzero shift | Achieved in restricted p=0 branches |
| Preserve laboratory atomic behavior | Dimensionless clock observables match measured stability | Fixed-Coulomb branch mismatches under shared-speed assumptions; equal-response branch avoids α drift only |
| Keep physical size constant | An invariant metric and stable rulers define a constant volume | Unresolved; α-preserving redshift branch changes atomic rulers |
| Conserve energy | Radiation loss appears in an explicit field/matter budget | Demonstrated numerically for a reduced radiation–field Hamiltonian |
| Derive the distance formula and rate | Dynamics selects a law and κ independently | Two different histories yield different laws; κ remains fitted |
| Explain radiation cooling | One field predicts both redshift and thermal evolution | Achieved conditionally for an already Planckian bath |
| Explain the CMB | Predict its origin, temperature, spectrum and angular structure | Temperature can be matched but is degenerate with initial conditions; origin and structure remain open |

The next discriminating theoretical calculation is a relativistic matter-and-gravity completion of the a=1,d=2 branch, including what operationally fixes the universe's size. A frequency-dependent light response is another retained branch, but would require causal dispersion, atomic spectra and wavelength-dependent redshift calculations rather than an assumed exemption for bound atoms. Neither branch has passed a Cassini residual analysis. No galaxy-rotation unification is claimed.

## Sources and reproducibility

1. Lange et al., *Improved limits for violations of local position invariance from atomic clock comparisons*, PRL 126, 011102 (2021). https://arxiv.org/abs/2010.06620
2. van de Bruck, Morrice & Vu, *Constraints on Disformal Couplings from the Properties of the Cosmic Microwave Background Radiation* (2013). Related theory, not validation of our candidate. https://arxiv.org/abs/1303.1773
3. NASA LAMBDA, FIRAS monopole spectrum and residual data. https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html
4. Fixsen, *The Temperature of the Cosmic Microwave Background* (2009). https://arxiv.org/abs/0911.1955
5. NIST, definitions and role of fundamental constants. https://physics.nist.gov/cuu/Constants/introduction.html

`calculate.py`, `firas.txt` and `results.json` accompany this report. Requires Python, NumPy and SciPy. The original rounded user coefficient is used explicitly; no refitting of that coefficient or additional galaxy data was done in this step.
