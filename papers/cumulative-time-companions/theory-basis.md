# Photon-Energy Transfer and Companion Deposition

## A phenomenological framework for redshift, additional gravity and reservoir emission in a nonexpanding universe

Working paper v1.0 | 13 September 2026 | Photon-Graviton research project

### Abstract

We formulate a hypothetical nonexpanding-universe model in which photons transfer energy into traveling companion excitations, a fraction of which is captured into gravitationally active reservoirs. An empirical one-third retention prescription supplies an effective deposited density, while optional settling, field-response and emission mechanisms are investigated separately. Existing calculations improve galaxy rotation predictions over ordinary matter alone, but the tested MOND relation performs better in aggregate. Cluster shear profiles can be fitted with free normalizations; transfer of the galaxy normalization and energy supply is not established. Six lens galaxies retain approximately 14% lens-angle discrepancies. A common phase-selector parameter set improves inner Milky Way circular-speed fits under two ordinary-matter baselines while leaving outer predictions unchanged. Conditional energy accounting permits the selected settling if released binding energy escapes. An optional nuclear-emission branch reinterprets black-hole-environment radiation as companion-powered. None of these findings derives the microscopic photon interaction, explains spectral redshift and event-duration stretching together, or establishes superiority to dark matter. We present explicit postulates, equation provenance, reproducible results and discriminating outcomes required for a stronger theory.

### 1. Motivation and scientific scope

The objective is to ask whether energy originating in radiation can account for some phenomena conventionally attributed to additional gravitating matter, while a related propagation process accounts for redshift without metric expansion. The nonexpanding universe is a stipulated setting of the model, not an observational conclusion of this paper. Observed spectra, stellar motions, lensing images and radiation remain constraints to explain. No particular universe age or total size is imposed.

The absence of an identified dark-matter particle motivates alternatives but is not evidence for this specific mechanism. For example, the LZ collaboration's December 2025 analysis reported no WIMP signal in its stated search range [1]. That result constrains particle properties rather than erasing gravitational evidence. Our companions are likewise unobserved. Calling them radiation-derived rather than primordial does not exempt them from the evidential burden placed on an additional gravitating component.

We distinguish three comparisons: ordinary-matter Newtonian dynamics as a control; specified MOND and dark-halo models as gravitational competitors; and expanding cosmology as a separate propagation and large-scale benchmark. A successful galaxy fit alone cannot establish superiority over all three. Published distances may be adopted as scenario inputs, but their calibration and any dependence on the competing model must be recorded. A distance inferred from the same redshift law under test cannot serve as independent validation of that law.

### 2. Core postulates and unresolved choices

- P1. Photons progressively transfer energy into a companion sector. A surviving photon can have lower energy; simple removal of whole photons is insufficient.
- P2. Traveling companions are taken to propagate at light speed and retain their energy until an interaction. Their microscopic identity, polarization and relation to gravitons remain unspecified.
- P3. Capture transfers energy and momentum into a bound reservoir. The reservoir adds gravity under the reference effective-density prescription. Ordinary matter continues to gravitate.
- P4. The calibrated retention fraction has an exponent fixed to one-third. It is an empirical organizing rule, not a derived fundamental constant.
- P5. A captured reservoir may redistribute internally. Movement to deeper binding releases energy that must be accounted for together with support, heat and outgoing flux.
- P6. Companion-powered emission in black-hole environments is an optional branch. Evaporation of stored rest energy and release of settling energy are distinct mechanisms and must not be counted together as the same available budget.

An environmental time field remains a possible cause of P1, but no unique environmental clock law is established. Low gravity, absolute potential, density, curvature and tidal field are different candidate environmental variables. Substituting one for another changes predictions and requires an explicit choice. A freely chosen time coordinate by itself is not an observable interaction.

### 3. Minimal equations and provenance

#### 3.1 Photon energy and redshift

Assume a local fractional energy-loss coefficient alpha per unit path length. With the standard photon energy-frequency relation, the deterministic reference law is

$$
\frac{dE_\gamma}{d\ell}=-\alpha E_\gamma,\qquad 1+z=\frac{E_{\rm emit}}{E_{\rm obs}}=\exp\left(\int\alpha\,d\ell\right).\tag{1}
$$

Provenance: the exponential solution is known mathematics and appears in other redshift models. Our proposed novelty would be an independently predictive physical equation for alpha, not equation (1). The archived illustrative constant is 0.0002488993 per Mpc; it does not establish broad-spectrum or independent distance-redshift agreement. At 100 million light-years it gives z approximately 0.00766 for the propagation component alone. Source/observer motions and ordinary endpoint gravitational effects must be treated separately.

An arrival-time map is also required:

$$
t_o=F(t_e,\mathrm{path}),\qquad \frac{\Delta t_o}{\Delta t_e}\simeq\frac{\partial F}{\partial t_e}.\tag{2}
$$

Provenance: a general timing identity. Equation (1) does not determine equation (2). A static delay common to both ends of an event does not stretch its duration. The observed supernova width-redshift relationship [2] must be reproduced without using it merely to assign a new arbitrary factor. Different responses of electromagnetic and gravitational signals must also meet the GW170817 timing constraint [3], including source emission-delay uncertainty.

#### 3.2 Transport and capture

For energy densities u and energy fluxes F, a minimal conversion/capture ledger is

$$
\partial_t u_\gamma+\nabla\cdot\mathbf{F}_\gamma=j_\star-Q_{\gamma c},\qquad \partial_t u_c+\nabla\cdot\mathbf{F}_c=Q_{\gamma c}-C.\tag{3}
$$

$$
\partial_t u_d+\nabla\cdot\mathbf{F}_d=C-S_d.\tag{4}
$$

Provenance: standard continuity equations; identifying the exchange terms with companions is proposed. S_d represents an actual energy transfer out of the deposited component and must enter another component or a boundary flux. If settling does gravitational work, kinetic and field-energy terms must be included; these three equations alone are not the full ledger. In a relativistic completion the total stress-energy exchanges must balance. This is not a claim that general relativity has a unique local gravitational energy density.

#### 3.3 Retention and the reference deposited profile

$$
X=\frac{L_{3.6}/(10^9L_\odot)}{(R_d/\mathrm{kpc})^2},\qquad \eta(X)=\frac{X^{1/3}}{1+X^{1/3}}.\tag{5}
$$

Provenance: the Hill-function form is established. Its radiation-density proxy, exponent choice and capture interpretation are project hypotheses. Here L refers to the adopted 3.6-micron luminosity proxy, not a measured total historical incident energy.

$$
a=q_aR_d,\qquad \kappa(r)=\frac{k_0}{[1+(r/a)^2]^2},\qquad J(r)=\frac{1}{2}\int_{-1}^{1}e^{-\tau(r,\mu)}d\mu.\tag{6}
$$

The optical depth tau is the integral of kappa along the incoming ray. The effective stored density is

$$
\rho_d(r)=\frac{2C_0\eta(X)J(r)}{[1+(r/a)^2]^2},\qquad M_d(<r)=4\pi\int_0^r\rho_d(x)x^2dx.\tag{7}
$$

Provenance: ray attenuation, angular averaging and spherical integration are known. The opacity shape and its stored-density interpretation are proposed. The current fitted constants are C0=4.72858924e7 solar masses/kpc cubed, k0=0.2039029004/kpc and qa=2.770766589. C0 is accumulated effective normalization, not a measured companion supply or universe age. The capture radius factor qa differs from the contraction factor s below.

$$
v_c^2(r)=v_b^2(r)+\frac{GM_d(<r)}{r}.\tag{8}
$$

Provenance: standard Newtonian spherical-source gravity with the proposed effective mass inserted. Equating a deposited energy density to mass density via energy divided by c squared is only a nonrelativistic reference closure. Traveling radiation, pressure and field stress require a more complete treatment. No adjustable energy-to-gravity amplification is derived here.

### 4. Optional physical branches retained for investigation

#### 4.1 Phase-dependent settling and its energy outlet

The current selector uses an ideal homogeneous Bose-gas criterion:

$$
T=\frac{m\sigma^2}{k_B},\quad T_c=\frac{2\pi\hbar^2}{mk_B}\left[\frac{\rho_d/m}{\zeta(3/2)}\right]^{2/3},\quad f=\max[0,1-(T/T_c)^{3/2}].\tag{9}
$$

Provenance: known ideal-gas statistics, proposed companion application. The adopted kinetic scale is G times enclosed monopole mass divided by 3r; it is not a measured temperature. Fraction f of each original shell is moved to sr and the rest stays. A shared inner-bin fit gives m=17.78 eV/c squared and s=0.8032. These values are effective fit parameters. Massive bound states would have to differ from the original light-speed traveling excitations; thermalization and particle-number behavior are unresolved. Published superfluid models motivate investigating phases [4] but do not derive this recipe.

With companion self-energy Ud, external-potential energy Ub and necessary global support energy K, the cooling requirement is

$$
E_{\rm mech}=K+U_b+U_d,\qquad Q_{\rm settle}=E_{\rm mech,i}-E_{\rm mech,f}.\tag{10}
$$

Provenance: standard energy accounting. Ud includes the self-energy factor one-half; Ub does not. K is assigned from the scalar force virial, a necessary global condition rather than a local equilibrium or stability proof. The tested shared model releases about 4.5e50 J, or about 1.2e-8 of deposited rest energy. A fully cooled solution requires this energy to leave the bound subsystem. All tested contractions pass the nonnegative-release condition, so conservation does not explain the fitted stopping radius.

#### 4.2 Refraction and relativistic completion

$$
\nabla\cdot[\epsilon(\rho)\nabla\Phi]=4\pi G(\rho_b+\rho_d).\tag{11}
$$

Provenance: density-dependent gravitational permittivity is an existing refracted-gravity approach [5]. Our companion-density driver and retained deposition profile are adaptations. Tested variants do not improve reserved outer Milky Way predictions, so refraction remains a secondary experimental branch. A nonrelativistic Phi does not alone determine light bending: a conventional weak-field metric involves a second potential Psi, and lensing probes their sum. No independent lensing multiplier is adopted.

#### 4.3 Microscopic conversion and environmental time

Photon-graviton mixing in an external magnetic field is an established theoretical process [6], but its ordinary coupling is weak and stationary whole-photon conversion does not redshift surviving photons. A proposed inelastic interaction or time-dependent field must calculate frequency change, photon-number change, momentum transfer, angular scattering and timing together. QED-like and scalar/axion effective-field descriptions remain candidates; string-inspired fields are a source of possible interactions, not numerical predictions without specified couplings [7]. None currently derives equations (1), (5) and (7) together.

For an illustrative Poisson process removing fraction delta per event, mean event count Lambda gives mean energy E0 exp(-Lambda delta) and relative energy variance exp(Lambda delta squared)-1. These are known Poisson identities applied to a hypothesis. At the archived alpha and a 100-million-light-year path, an illustrative added relative linewidth of 1e-5 requires per-event loss below 1.31e-8 and at least roughly 582,000 events. This is a design example, not an observed exclusion limit.

#### 4.4 Companion-powered nuclear emission

$$
L_\gamma=f_\gamma P_{\rm release},\qquad \dot E_{\rm res}=P_{\rm capture}-P_{\rm release}.\tag{12}
$$

Provenance: ordinary reservoir bookkeeping; attributing nuclear power to companions is our optional hypothesis. All released channels must be included. The radiated fraction f-gamma is not the retention fraction eta. Binding-energy release preserves constituent rest mass at Newtonian order; rest-energy evaporation changes the permanence postulate unless replenished. Observed black-hole-environment light is not established Hawking radiation, and no escape from inside a conventional horizon is assumed.

At a Sagittarius A* reference luminosity of 1e29 W or less [8], pure photon output corresponds to about 1.77e-11 solar masses/year of energy. The representative M87 core luminosity 2.7e35 W [9] corresponds to 4.77e-5 solar masses/year; its jet beaming, aperture and excluded mechanical power limit that comparison. These values reconstruct required throughput from observed brightness, rather than predicting brightness. Routing the whole Galactic settling budget through the faint Sgr A* reference channel would take about 1.43e14 years. This ratio imposes no universe age: it instead demands an explicit routing, efficiency and activity-history model.

### 5. Executed results and limits of the evidence

The project record contains both successes and failures. Parameter selection and previously inspected observations must not be described as blind tests. Reported RMS values are descriptive, not likelihood-based evidence for the entire theory.

| Test | Current result | What it establishes |
|---|---|---|
| SPARC rotation, 149 galaxies | Companion validation/test RMS 32.49/23.59 km/s; MOND 26.88/16.40 | Companions improve the ordinary-only control, but tested MOND is better |
| Six Coma shear bins | Companion/NFW omitted-bin sums 6.86/6.82 | Comparable freely fitted shapes; no transferred supply normalization |
| Six lens galaxies | About 14% lens-angle RMS | A material unresolved motion/lensing discrepancy |
| Shared Milky Way settling | Inner RMS 3.68/11.83 becomes 2.74/6.70 km/s | A shared conditional inner-fit improvement under two matter baselines |
| Same Milky Way outer bins | RMS remains 8.36/9.78 km/s | No outer improvement from the selected central rearrangement |
| Refraction sweep | No selected outer improvement | No reason to promote the tested refraction family |
| Nuclear release | Energy-throughput estimates only | A possible diagnostic channel, not a predicted spectrum or luminosity |

SPARC inputs [10] have 89 training, 29 validation and 31 test galaxies in the project split; all have been exposed during development. Companions have three shared fitted constants versus one for the tested MOND law. Restricted shared NFW mappings and per-galaxy NFW fits use different information budgets and cannot be pooled into one ranking. The Milky Way phase pilot uses a spherical gas monopole and two archived stellar baselines, with 20 inner bins for selection and 18 outer bins for comparison. It is not a comparison of two independent galaxies or individual bulge-star orbital predictions [11].

Earlier wave/particle support, saturation, inward migration and coupled exchange tracks remain in the research record. None supplies a demonstrated joint solution for stability, absolute supply and lensing. Numerical mass/energy checks and refinement protect against implementation errors; they do not verify the physical assumptions. The accompanying cross-scale-performance.md and linked executable reports contain the detailed audit trail.

### 6. Evidence that would materially strengthen the theory

The following outcomes are the main work program. Each explains what is needed and why failure matters. They are not claims already achieved.

1. A common mechanism with fewer adjustable choices. Specify an interaction and local state equations that yield conversion, capture and support with shared constants. Derive the one-third exponent or replace it with a better justified prediction. Otherwise independent fitting rules can conceal incompatible physics.

2. Absolute supply from independently constrained radiation. Integrate source luminosities, travel, capture, retention and release, keeping age/history as explicit inferred parameters rather than choosing unlimited history to rescue a deficit. Predict effective stored gravity without fitting C0 separately to the needed mass. Otherwise we have renamed a missing-mass profile rather than explained its origin.

3. Redshift and time stretching from the same propagation law. Use distances independent of redshift, control velocities, predict the same fractional spectral shift across tested bands, and reproduce event durations, sharp images and brightness. Include GW/EM timing. Otherwise matching a distance-redshift curve alone does not explain the observed propagation phenomena.

4. Joint motion and lensing prediction. Use one relativistic response and one deposition profile to predict circular speeds, stellar velocity distributions and strong/weak lensing in the same objects. Include distances, gas, mass-to-light and orbit uncertainty. Otherwise a rotation success cannot repair a lensing failure through an unrelated correction.

5. Genuine transfer to untouched systems. Freeze parameters and analysis choices, then predict additional galaxies and clusters spanning surface brightness, gas content, morphology and environment. Compare residuals, scatter and parameter burden against properly specified MOND and dark-halo models with identical inputs. Otherwise in-sample improvement may reflect flexibility or selection.

6. Formation, stopping and persistence. Solve evolving density, pressure or orbital support, angular momentum and cooling, including feedback from the companion field. Check Solar System and binary constraints in the appropriate limits. Otherwise a useful profile may never form, may collapse, or may disrupt systems already explained well.

7. A radiation prediction distinct from ordinary sources. From the same fixed settling law, predict where energy emerges, its spectrum and variability before using measured brightness. Test both diffuse emission and the optional nuclear channel, including nondetections and jet mechanical power. Otherwise any observed glow can be assigned to companions after the fact.

8. Broader nonexpanding-universe consistency. Address the measured microwave spectrum and angular correlations, large-scale clustering, abundances and thermal histories without importing an expanding-model solution as a premise. Retain observations even when their conventional explanation is set aside. Otherwise a gravitational phenomenology cannot claim to replace a cosmological framework.

### 7. What would justify a claim of advantage?

A restricted claim is the first realistic target: for example, a frozen companion law predicts joint motions and lensing in a new sample more accurately than specified competitors with comparable or lower effective flexibility. That would be a meaningful result even before every cosmic observable is explained. It should be reported with measurement covariance, nuisance-parameter uncertainty, model complexity and all failed targets included.

A stronger advantage would come from an independently confirmed connection between lost photon energy, accumulated gravity and released radiation. An environmental dependence or radiation feature predicted in advance, absent from the tested competing models, would be especially valuable. The signature must be quantitatively calculated; an unexplained residual is not automatically companion evidence.

Directly detecting an interaction or state associated with companions would be compelling, but is not the sole route to a viable theory. Conversely, dark-matter direct-detection limits do not automatically promote this theory. Both frameworks must meet the evidence appropriate to their claims. Since companions add an unseen gravitational component, the distinct scientific claim is its radiation-derived origin and calculable dynamics, not simply the absence of unseen matter.

### 8. Continuing work across tracks

The reference track remains exact-third deposition with transparent input and parameter accounting. The settling track seeks local support and cooling equations that determine the contraction rather than fitting it. The propagation track seeks an inelastic or environmental-time interaction that predicts both spectra and event timing. The gravity track seeks one relativistic lensing/motion response; unsuccessful refraction variants remain documented. The emission track now includes companion-powered nuclear radiation as an optional branch while retaining diffuse and invisible outgoing channels. The data-comparison track freezes rules before applying them to additional systems.

These tracks must eventually share one energy and momentum ledger. Until then they are labeled alternatives or partial closures; results from incompatible versions must not be assembled into a fictitious single successful theory. No background run or completed physical goal is implied by retaining a track in this work program.

### 9. Conclusion

The current basis is a radiation-conditioned extra-gravity hypothesis with an empirical retention law, several reproducible phenomenological successes, and substantial unresolved mechanisms. Selective settling can improve inner Milky Way fits without increasing the outer force; nuclear brightness can constrain a proposed release channel. Neither establishes photon-derived gravity. The strongest next advance would be a conserved, locally supported model with fixed parameters that predicts a joint gravitational and radiative observation in systems not used to construct it. This is a working hypothesis and results paper, not a demonstrated replacement for dark matter, MOND or expanding cosmology.

### References

[1] LZ Collaboration / Berkeley Lab (2025). LZ Sets a World's Best in the Hunt for Galactic Dark Matter. [Experiment result and search scope](https://newscenter.lbl.gov/2025/12/08/lz-sets-a-worlds-best-in-the-hunt-for-galactic-dark-matter/).

[2] DES Collaboration (2024). Slow supernovae show cosmological time dilation out to z approximately 1. [arXiv:2406.05050](https://arxiv.org/abs/2406.05050). Used as an observed timing constraint; expansion is not adopted as this paper's premise.

[3] LIGO/Virgo and partner collaborations (2017). Gravitational Waves and Gamma-Rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A. [Publication record](https://dcc.ligo.org/LIGO-P1700308/public).

[4] Berezhiani, Famaey and Khoury (2018). Phenomenological consequences of superfluid dark matter with baryon-phonon coupling. [arXiv:1711.05748](https://arxiv.org/abs/1711.05748).

[5] Sanna, Matsakos and Diaferio (2023). Covariant formulation of refracted gravity. [arXiv:2109.11217](https://arxiv.org/abs/2109.11217).

[6] Palessandro and Rothman (2023). A Simple Derivation of the Gertsenshtein Effect. [arXiv:2301.02072](https://arxiv.org/abs/2301.02072).

[7] Cicoli, Goodsell and Ringwald (2012). The type IIB string axiverse and its low-energy phenomenology. [arXiv:1206.0819](https://arxiv.org/abs/1206.0819).

[8] Event Horizon Telescope Collaboration (2022). First Sagittarius A* Results I. [doi:10.3847/2041-8213/ac6674](https://doi.org/10.3847/2041-8213/ac6674).

[9] Prieto et al. (2016). The central parsecs of M87: jet emission and an elusive accretion disc. [MNRAS 457, 3801](https://academic.oup.com/mnras/article/457/4/3801/2588956).

[10] Lelli, McGaugh and Schombert; SPARC data and related publications. [Project publications](https://astroweb.cwru.edu/SPARC/publications.html). Project-specific splits and fits are documented in the repository.

[11] Eilers et al. (2019). The Circular Velocity Curve of the Milky Way from 5 to 25 kpc. [arXiv:1810.09466](https://arxiv.org/abs/1810.09466).

### Reproducibility and research status

This manuscript summarizes the repository at main revision 35a7ea7 and the preceding analysis commits. Supporting records are in research_work/results/isotropic-galaxy-transfer/ and research_work/results/companion-extensions/. The latter includes refraction.py, phase.py, settling.py, their result files, numerical checks, protocols and nuclear-release.md. The companion paper cross-scale-performance.md links the complete cross-scale audit. Historical manuscript PDFs are preserved separately. This document introduces no new numerical fit, independent observational discovery, assigned external authorship or peer-review claim. Computational assistance was used in drafting and analysis; independent scientific review remains necessary.
