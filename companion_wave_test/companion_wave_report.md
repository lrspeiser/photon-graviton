# Companion waves: equations and first data tests

**Exploratory study, 9 September 2026.** This develops the proposed fictional universe using observed galaxy and cluster data. It does not establish a new interaction or a replacement for dark matter. The most promising feature in these tests is collection proportional to area combined with a broad, cored deposit distribution. The largest difficulty is the energy budget if deposited energy gravitates with ordinary strength.

There is no finite list of all possible theories. The following equations form a modular framework covering direct stored energy, amplified gravitational response, external collection, local production, finite residence, and redistribution. Each module introduces assumptions that must eventually be derived from one interaction.

## 1. What gravity lets us infer—and what it does not

For a circular orbit in the weak-field limit, define the inward acceleration magnitude

\[
g_{\rm obs}(r)=v_{\rm obs}^2(r)/r,\quad
g_b(r)=v_b^2(r)/r,\quad \Delta g=g_{\rm obs}-g_b.
\]

If the extra source is approximately spherical and obeys inverse-square gravity,

\[
\boxed{M_X(<r)=\frac{r^2\Delta g(r)}{G}
=\frac{r[v_{\rm obs}^2(r)-v_b^2(r)]}{G}.}
\]

The baryonic disk need not be spherical: its radial force is supplied separately by the SPARC mass model. Only the extra component is treated as spherical. M_X is an *equivalent gravitational source*, not a direct weighing of companion waves.

Suppose deposited energy density u_d has effective gravitational density ρ_X=ηu_d/c², with constant η>0. Then

\[
\boxed{E_d(<r)=\frac{c^2M_X(<r)}{\eta}},\qquad
\boxed{u_d(r)=\frac{c^2}{4\pi G\eta r^2}
\frac{d}{dr}\{r[v_{\rm obs}^2-v_b^2]\}.}
\]

η=1 is the ordinary rest-energy-equivalent budget benchmark. A bath of relativistic waves and its confining stresses require a relativistic calculation; this benchmark is not an exact general-relativistic treatment of trapped gravitational radiation. Factors of order unity from pressure cannot remove a shortage of thousands or billions.

η≫1 instead specifies a modified response, not a mechanism for creating extra energy. If η varies with position, the local density equation uses η(r), and E_d is an integral; E_d=c²M_X/η cannot use an arbitrary single average.

Gravity alone measures ηu_d in this prescription. It cannot separately determine η, incoming power, capture efficiency, or storage time. A negative inferred density would be incompatible with a strictly positive deposited source under these assumptions, but noisy differentiation and uncertain baryonic forces must be addressed before interpreting that as a physical failure.

For clusters the corresponding weak-field hydrostatic estimator is

\[
M_{\rm HSE}(<r)=-\frac{r^2}{G\rho_g(r)}\frac{dP_{\rm th}}{dr},
\qquad M_X=M_{\rm HSE}-M_{\rm gas}-M_\star.
\]

Nonthermal pressure changes this estimator. Hot gas is part of the baryon subtraction, not missing from it.

## 2. Photon energy transferred into companions

A phenomenological continuous loss rule is

\[
\frac{dE_\gamma}{ds}=-a(s,t,\nu)E_\gamma,
\qquad E_\gamma(s)=E_e\exp[-\int a\,ds].
\]

With identical, unchanged atomic reference frequencies and no additional endpoint shifts,

\[
1+z=E_e/E_o,\qquad E_{c,\rm generated}=E_e-E_o
=E_e\frac{z}{1+z}.
\]

The second equation requires that *all* of the lost photon energy enters the companion channel. A branching fraction b_c between zero and one can multiply it. The accompanying field, photons and matter must share one conserved total stress-energy tensor. Counting E_e plus the transferred energy would double-count energy.

The original candidate a=κ gives

\[
1+z=e^{\kappa s},\quad f_c(s)=1-e^{-\kappa s},\quad
\kappa=7.7315\times10^{-5}\ {\rm Mly}^{-1}.
\]

κ is retained from the supplied hypothesis; these gravity tests do not refit or validate it. Across one kiloparsec, the converted fraction is only about 2.52×10⁻⁷. For the energy tests below, s is explicitly set equal to the last measured rotation-curve radius; it is an illustrative local path, not a radiative-transfer average over the galaxy.

**A necessary distinction:** a static energy-loss rule does not by itself stretch arrival intervals. In a stationary propagation medium, successive identical flashes can suffer the same energy loss and still arrive with their original separation. The earlier time-history model remains a separate possible completion:

\[
\omega_\gamma=ck/n(t),\quad \dot n=\gamma n^p,
\quad 1+z=n_o/n_e,\quad dt_o/dt_e=n_o/n_e.
\]

Here k is conserved spatial wave number and the same ray equation determines color and duration. Its wave energy exchange rate is Q_γ=(ṅ/n)u_γ. Setting the tensor-wave dispersion to the same expression would be an additional shared-propagation requirement. Neither a companion interaction nor protection of atomic/cavity ratios has been derived by writing these equations.

The earlier improved history p≈0.264 is not silently substituted into a universal path-loss coefficient: a time history and a stationary distance law are different mechanisms. Their local, present-day energy-loss rates can agree without agreeing globally.

## 3. Energy-conserving transport and storage

In a static spatial background, a minimal coarse-grained bookkeeping system is

\[
\partial_tu_\gamma+\nabla\cdot\mathbf F_\gamma=j_\star-Q,
\]
\[
\partial_tu_c+\nabla\cdot\mathbf F_c
=Q-q_{\rm cap}+(1-f_h)u_d/\tau,
\]
\[
\partial_tu_d+\nabla\cdot\mathbf F_d=q_{\rm cap}-u_d/\tau,
\qquad \partial_tu_h=f_hu_d/\tau.
\]

The emitter reservoir loses j_★. Here u_h records thermalized energy, 0≤f_h≤1, τ is the residence time, and F_d can describe redistribution, for example −D∇u_d with D≥0. Adding these equations conserves energy including boundary fluxes and the emitter reservoir. This is bookkeeping, not yet a microscopic wave action; momentum, pressure, gravity and material response must also be supplied.

For absorption coefficient α_c≥0 and companion speed v_c,

\[
q_{\rm cap}=\alpha_c v_c u_c,\quad
\epsilon_{\rm ray}=1-e^{-\int\alpha_c ds}.
\]

For isotropic radiation of energy density u_c incident on a sphere of radius R, total incoming power is πR²v_cu_c. Thus

\[
\boxed{P_{\rm cap}=\epsilon\pi R^2v_cu_c},\quad 0\le\epsilon\le1,
\]

where ε is an angularly averaged capture fraction. This is the mathematical expression of the larger-collector idea. It assumes an external bath; gravitational focusing would require a separately calculated cross-section, not an arbitrary increase of ε above one.

For constant captured power, no initial stored energy and no net redistribution flux through the outer boundary,

\[
E_d(T)=P_{\rm cap}\tau(1-e^{-T/\tau}),
\qquad T_{\rm eff}=\tau(1-e^{-T/\tau})\le T.
\]

Therefore

\[
\boxed{M_X=\frac{\eta\epsilon\pi R^2v_cu_c T_{\rm eff}}{c^2}},
\qquad
\boxed{u_{c,\rm req}=\frac{c^2M_X}{\eta\epsilon\pi R^2v_cT_{\rm eff}}.}
\]

For time-varying production/capture, replace P_cap T_eff by

\[
E_d(T)=E_d(0)e^{-T/\tau}
+\int_0^T P_{\rm cap}(t)e^{-(T-t)/\tau}dt.
\]

Locally produced companions give an alternative supply:

\[
P_{c,\rm local}=b_c\int j_\star(\mathbf x)\langle f_c(\ell)\rangle d^3x,
\qquad M_X=\eta\epsilon P_{c,\rm local}T_{\rm eff}/c^2.
\]

An indefinitely old fictional universe can change the accumulation budget. That does not automatically supply infinite energy: finite τ limits storage, emitters have finite fuel unless replenished, and the external bath must itself satisfy a source/sink balance. Existing deposited energy from earlier epochs is another explicit initial condition.

## 4. Where must the energy be deposited?

If v_X²=v_obs²−v_b² is roughly constant in an outer region,

\[
M_X(<r)\propto r,\qquad \rho_X(r)\propto r^{-2}.
\]

A useful finite-central-density template is

\[
\rho_X(r)=\frac{\rho_0}{1+(r/r_c)^2},
\]
\[
M_X(<r)=4\pi\rho_0r_c^2[r-r_c\arctan(r/r_c)],
\quad
v_X^2(r)=4\pi G\rho_0r_c^2[1-(r_c/r)\arctan(r/r_c)].
\]

It must eventually truncate or steepen: its mass otherwise grows without bound. For an outer boundary R, one normalized deposition prescription producing this profile with constant residence time is

\[
q_{\rm cap}(r)=\frac{P_{\rm cap}}
{4\pi[R-r_c\arctan(R/r_c)](r^2+r_c^2)},\quad 0\le r\le R.
\]

Its volume integral is P_cap. This is an explicit *target profile*, not a derived capture law. If incoming energy is instead absorbed only in a skin, some redistribution mechanism must create this broad profile. Diffusion alone with a specified outer source does not automatically produce it.

A spherical shell entirely outside an orbit exerts no extra interior Newtonian force. Deposits need to lie within the radii where extra inward acceleration is required. Faster circular motion also does not necessarily mean infall: collection needs angular-momentum transport and energy dissipation.

## 5. New SPARC calculations

Source: [SPARC data release](https://astroweb.cwru.edu/SPARC/), [Lelli, McGaugh & Schombert 2016](https://arxiv.org/abs/1606.09251).

Starting with the 175-object catalog, the existing selection retains **149 galaxies and 3,150 rotation measurements**: quality flag 1 or 2, inclination ≥30°, positive disk scale length, at least five valid rows. Distances and inclinations are held at catalog values. Disk and bulge mass-to-light ratios are 0.5 and 0.7. Gas contributions retain the signed squared-velocity convention. The 26 exclusions and input hashes are recorded in results.json.

### Radial shape test

Each non-baryonic template gets one nonnegative amplitude per galaxy, fitted to the entire curve by minimizing mean squared velocity residual. The cored template fixes r_c to the observed stellar disk scale length; it does not fit a second halo radius. The shell occupies 0.9R_max to R_max with uniform density in that layer. Reported errors weight galaxies equally.

| Extra-source shape | Velocity RMSE, km/s |
|---|---:|
| Baryons only | 52.77 |
| Central point source | 44.86 |
| Outer 10%-radius shell | 48.64 |
| Uniform volume density | 26.55 |
| Density proportional to r⁻² | 20.89 |
| Cored r⁻² profile | **12.29** |

These are **descriptive fits**, not predictive validation or likelihood significances. All curves were used to fit their own amplitudes. They show which deposited-source shapes can resemble the data; they do not show that photon capture creates those shapes. The shell comparison is conditional on spherical extra gravity.

### Energy test

All 149 selected galaxies have positive inferred extra mass at the outermost retained point. Define a transparent luminosity conversion parameter

\[
L_{\rm bol}=B\,\ell_{3.6}L_\odot,
\]

where ℓ_3.6 is the numerical catalog luminosity in solar 3.6-μm units and L_⊙=3.828×10²⁶ W is the bolometric solar luminosity. **B=1 is a stated proxy, not measured bolometric photometry.** Results scale inversely with B and with the assumed illumination duration. A historical luminosity integral is required for an actual formation history.

With B=1, T=10 billion years, perfect retention and capture, and η=1:

| Quantity, median over galaxies | Result |
|---|---:|
| Required equivalent mass inside measured outer radius | 2.36×10¹⁰ solar masses |
| Time to supply required energy at proxy luminosity, even at 100% conversion | 5.77×10¹³ years |
| Required energy / all proxy starlight emitted over T | **5,775** |
| Companion fraction generated across one R_max | 3.18×10⁻⁶ |
| Required energy / local companion supply over T | **1.72×10⁹** |
| Required external companion energy density, v_c=c, capture=1 | **9.45×10⁻¹¹ J/m³** |

The external density is *required*, not observed. It is a possible target for an independent background-energy calculation; it must not be presented as a measured companion bath. Less than perfect capture or retention worsens the requirements. External illumination from many galaxies is not included in the local luminosity budget and must be calculated separately.

At B=10 the first energy shortage falls to about 577, and the local-production shortage to 1.72×10⁸. Thus a moderate luminosity-conversion change is insufficient. Changing disk mass-to-light ratio to 0.3 or 0.8, with bulge 0.2 higher, changes the median full-conversion shortage to 6,153 or 5,275. These are sensitivity tests, not fits to population-synthesis data.

Examples of the inferred ordinary-equivalent energy:

| Galaxy | Measured outer radius | Extra mass | E_required at η=1 |
|---|---:|---:|---:|
| DDO154 | 5.92 kpc | 2.49×10⁹ M_⊙ | 4.44×10⁵⁶ J |
| NGC2403 | 20.87 kpc | 7.43×10¹⁰ M_⊙ | 1.33×10⁵⁸ J |
| NGC3198 | 44.08 kpc | 1.85×10¹¹ M_⊙ | 3.30×10⁵⁸ J |

### Across-galaxy collection scaling

Fit one common normalization to outer extra masses in the prior 89-galaxy training subset. Compare M_X proportional to luminosity, R_max², or R_max³. The previously exposed 29-object validation and 31-object test subsets are reused; this is **exploratory reuse, not a new blind holdout**. All fits minimize squared log-mass residuals.

| Scaling | Training RMSE, dex | Validation RMSE, dex | Reused test RMSE, dex |
|---|---:|---:|---:|
| M_X ∝ luminosity | 0.468 | 0.532 | 0.530 |
| M_X ∝ R_max² | **0.269** | **0.215** | **0.262** |
| M_X ∝ R_max³ | 0.462 | 0.350 | 0.391 |

The area normalization is

\[
\boxed{M_X(<R_{\max})\approx1.390\times10^8M_\odot
(R_{\max}/{\rm kpc})^2.}
\]

A 0.262-dex RMS log residual corresponds to a multiplicative scale of about 1.83; it is not a measured 68% interval. This favors an area *description* among these three templates. R_max is an observing limit, not a measured capture boundary. M_X is itself constructed from radius and speed, and radius correlates with galaxy structure and luminosity. Consequently this is not independent evidence for incoming waves. The area scaling is across objects; it must not be substituted as M_X(<r)∝r² within each galaxy, which would contradict a flat extra rotation contribution.

## 6. Cluster transfer check

Sources: [Ettori et al., hydrostatic masses, Table 1](https://arxiv.org/abs/1805.00035), [Eckert et al., hydrostatic gas fractions, Table 2](https://arxiv.org/abs/1805.00034).

We transcribed R_500, M_HSE,500, and measured hydrostatic gas fractions for **11 clusters**. HydraA has no gas-fraction entry in the second table. A2029 has different mass values across the tables and is excluded rather than silently combined. Stellar mass fractions of 1%, 2% and 3% of total mass are explicit sensitivity assumptions; we have not measured individual cluster stellar masses here. The central estimate uses

\[
M_X=M_{500}(1-f_{\rm gas}-0.02).
\]

No correction that forces the cosmic baryon fraction is used. This avoids importing that particular cosmological constraint, but published distance conversion, hydrostatic equilibrium, and fitted mass-profile assumptions remain. R_500 is a conventional overdensity-defined radius, not a demonstrated absorption boundary in the fictional static universe.

**Using the galaxy area normalization unchanged**, required cluster excess exceeds prediction by a median factor **2.21**, spanning **1.92–2.63**. Changing the assumed stellar fraction from 1% to 3% changes the median factor from 2.24 to 2.19. The median extra mass is 4.69×10¹⁴ M_⊙, or 8.38×10⁶¹ J at η=1.

**Clarification: “underprediction” applies only to the deliberately fixed galaxy normalization, not to every possible companion-wave theory.** We did not independently measure the incoming companion flux. We inferred a combined normalization from galaxy gravity and asked whether it transfers unchanged. Allowing more stored energy or a stronger response can match the cluster requirement, but then a rule predicting that change is needed. Size alone therefore does not give one shared normalization in this check. The product η ε u_c T_eff would have to be about twice as large in clusters—or the radius definition and underlying gravitational inference would have to change. That is a requirement to explain, not a validated environmental law. We did not fit cluster radial profiles, lensing maps, a new static distance model, or a joint error covariance. These numbers are a transfer diagnostic rather than a formal exclusion significance.

## 7. Candidate families: advantages and remaining costs

| Family | Governing choice | What helps | What remains unresolved |
|---|---|---|---|
| Locally generated, ordinarily gravitating stored energy | η=1; P∝L f_c | Direct energy accounting | Severe energy shortage at specified duration; weak luminosity scaling |
| External collection with ordinary stored energy | P=επR²cu_c; η=1 | Best of three across-galaxy size/luminosity templates | Huge required, unmeasured bath; capture and confinement; cluster normalization |
| Very long storage/history | Large T_eff or initial E_d | Can increase available stored energy in an old static universe | Requires survival, replenished sources and global balance; not established by a fit |
| Amplified gravity response | η≫1 | Reduces actual energy needed for given extra force | η is a new coupling; needs stable dynamics, lensing and clock compatibility |
| Broad redistributed deposit | Cored density above; transport D | Best descriptive radial shape | Producing it from external illumination without arbitrary per-galaxy tuning |
| Thin outer-edge deposit | Capture in a skin | Natural possibility for an opaque collector | Poor radial fit; no interior shell force |

These choices can combine, but their parameters are degenerate. An unconstrained efficiency for every object could reproduce any positive inferred mass and would not test the hypothesis.

A modified-field alternative to literal stored mass could start with the weak-field equation

\[
\nabla\cdot[\mu(X,u_d)\nabla\Phi]=4\pi G\rho_b,
\quad X=|\nabla\Phi|/a_0.
\]

For a spherical system this becomes μg=g_b. A companion-dependent μ could amplify the response without assigning E_d=M_Xc². This is a candidate equation, not a derivation: one must specify μ globally and a field action with acceptable kinetic terms and perturbations. Merely defining μ=g_b/g_obs object by object is circular.

For instance, a weak-field variational starting point is a gravitational term proportional to −a_0²F(|∇Φ|²/a_0²,u_d)/(8πG), plus −ρ_bΦ; variation with respect to Φ gives the equation with μ=∂F/∂(|∇Φ|²/a_0²). This is only a static, nonrelativistic functional. It does not supply temporal kinetic terms, a causal tensor sector, or a conserved matter–light–gravity interaction, and its dependence on u_d must be included in the energy budget.

## 8. Tests still needed before calling this a theory

1. **Capture and residence:** predict α_c, τ and any redistribution coefficient from a common interaction, then solve the transport equation. Freely fitted radial energy profiles are not sufficient.
2. **Budget:** reconstruct luminosity histories and the external companion bath consistently with redshift energy loss. Compare a fixed common ηεu_cT_eff across galaxies, environments and clusters. Specify how the bath gravitates before capture and where energy goes after leakage.
3. **Lensing and motion together:** in a weak-field metric with potentials Φ and Ψ, slow matter follows Φ while deflection depends on the line-of-sight gradient of Φ+Ψ. Rotation curves alone do not determine Ψ. A matter-only extra force is not automatically extra lensing.
4. **Clocks and complete signals:** derive emission, ray propagation, detection, pulse duration and flux from the same interaction; retain the earlier atomic-ratio and cavity constraints. Photon energy loss alone is insufficient.
5. **Independent object tests:** predict amplitudes and radii without using those objects' rotation curves. Include distance/inclination and mass-to-light uncertainties, and cluster pressure and geometry uncertainties. Do not call the reused SPARC split blind.
6. **Dynamics:** account for energy and momentum, trapping, heating, angular momentum and stability. Familiar gravitational waves are not automatically absorbed into galactic rims. A companion channel with strong capture is new physics that must be specified.

**Recommended next candidate:** an externally supplied, broadly distributed companion field with a single capture/retention law. Carry an explicit response parameter η while testing η=1 first. Its attraction is the area scaling and broad radial shape; its unresolved cost is energy supply and coupling. This remains a nonluminous gravitational component even if it is a field rather than a population of dark particles.

## Reproducibility

The companion_wave_test archive contains raw SPARC inputs, the two primary cluster PDFs and extracted text, run.py, cluster_test.py, per-galaxy and per-cluster tables, and JSON summaries. Run run.py followed by cluster_test.py with Python, NumPy and SciPy installed. The cluster transcription is visible in cluster_test.py. Fixed catalog values and the proxy luminosity assumptions are intentional limitations, not uncertainties marginalized by the code. No new supernova, clock, CMB, lensing or photon–companion interaction fit was performed in this study.
