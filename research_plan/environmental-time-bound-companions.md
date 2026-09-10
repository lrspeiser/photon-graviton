# Environmental Time and Bound Companion Fields

**Historical branch:** the user subsequently selected [direct conversion without an added clock factor or screening](direct-conversion-bound-companions.md) as the active starting model. This document preserves the preceding hypotheses and their provenance. Its bound-mode concept carries forward; its time-rate prescription is not automatically part of the active branch.

Working research specification, 9 September 2026. This document consolidates the latest discussion; it is a hypothesis within a fictional nonexpanding universe, not an established branch of physics or a completed first-principles theory. The descriptive title makes no originality claim.

**Follow-up clarification and test:** the [screened-clock diagnostic](../research_work/results/void-screening/report.md) treats the environmental response in equations (3)-(4) as an additional factor multiplying a baseline clock law, rather than as the complete clock law. Its 180 illustrative cases show that accumulated local photon loss can be tiny without screening, while a slower-void clock factor can oppose ordinary gravitational gradients. A faster-void sign is tested as a separate alternative, not adopted. No clock amplitude or sign is established by redshift alone; the equations below preserve the earlier proposed branch for provenance.

## The concept in ordinary language

The universe has the observed stars, galaxies and astronomical signals of our own, but permits a different relationship between gravity, time and light. Deep gravitational wells suppress an environmental time effect. In weak gravitational environments, that effect enables light to transfer energy continuously into a companion field. The surviving photons become redder. Companions carry the transferred energy forward at the locally measured speed of light and do not collect in voids.

Inside sufficiently deep wells, an interaction can transfer traveling companion energy into long-lived bound configurations. These configurations are extended spatial patterns, loosely analogous to bound energy levels, rather than literal shells of spacetime. Incoming companions can populate those patterns. Their stored energy contributes to gravity and potentially to the motions of outer stars and the bending of background light by galaxies and clusters.

The central aim is to derive a single chain: environment -> conversion -> transport -> capture -> bound distribution -> gravitational observables. Giving each link a name does not yet derive the chain.

## Decisions and their current status

| Item | Status |
|---|---|
| Nonexpanding universe; no Big Bang or dark matter as active explanations | Adopted premise. Published distances may be stipulated, but their origin remains labeled. |
| Local clocks and local processes feel normal; local light speed is c | Required behavior. A complete matter/field coupling has not been selected. |
| A steady altered environment can continuously transfer photon energy | Adopted phenomenological choice. A static clock difference alone is not its derivation. |
| Total traveled path counts, including repeated ideal mirror reflections | Adopted idealization; reflection does not reset the transfer. A cavity wave calculation remains necessary. |
| Complete transfer of lost photon energy to companions | Adopted scalar energy-accounting assumption. Momentum, recoil and field work must also close in a microscopic theory. |
| Traveling companions move forward at c without further loss | Adopted propagation idealization. |
| Voids are transit regions, with no permanent capture reservoir | Adopted. Passing energy density is distinct from stored energy. |
| Neglect gravitational backreaction of traveling energy in the first calculation | Working approximation, to be checked later; energy is not counted twice. |
| Deep wells enable capture into bound configurations | Working hypothesis; threshold and interaction unspecified. |
| Extended bound levels or modes | Latest candidate, not a demonstrated spectrum or stability result. Wave-particle duality alone does not provide capture. |
| Extra gravity becomes relatively more noticeable outward | Preferred initial interpretation; enhanced capture specifically at edges remains an alternative to test. |
| Companion deposits have the distribution conventionally attributed to dark matter | Optional imposed benchmark, not a predicted map or an adopted dark-matter explanation. |
| No decay or escape of deposited energy | Idealized working limit. Its physical lifetime must be derived. Internal redistribution is allowed. |
| Total supply calculation | Deferred at the user's request; conservation remains part of the model. |
| Whether deposited gravity also regulates time/conversion/capture through W | Unresolved. Do not silently choose feedback or exclude it. |

## Formula provenance

Every formula below is classified as established mathematics, a proposed model rule, or a conditional consequence. Familiar equations do not become novel merely by changing the interpretation of a coefficient. No microscopic interaction or originality has been established here.

### 1. Environmental gravity-well depth

\[
W(\mathbf{x})=\Phi_{\rm ref}-\Phi_{\rm env}(\mathbf{x}). \tag{1}
\]

**Established potential-difference definition; proposed environmental use.** W is nonnegative within the selected wells and has units of energy per mass. A common reference and boundary prescription must be fixed; W is not simply local acceleration, which can vanish where opposing pulls cancel.

Initially, estimate the environmental potential from ordinary matter under an explicitly selected gravity law. Do not insert a dark-matter fit into that ordinary-matter map. Whether to include deposited companions in the environmental potential later is unresolved; the potential controlling motion need not be identified with this provisional environmental input without specifying the coupling.

For an isolated static Newtonian system with the reference at infinity:

\[
W=\tfrac12 v_{\rm escape}^{2}. \tag{2}
\]

**Established conditional relation.** It is not automatically valid for an entire nonisolated universe or arbitrary modified gravity. [Escape-speed derivation](https://galaxiesbook.org/chapters/I-02.-Elements-of-Classical-Mechanics_2-Escape-velocity.html).

### 2. Proposed time response and conversion coefficient

\[
q\equiv\frac{d\tau}{dt},\qquad
q(W)=1-\frac{\epsilon}{1+W/W_*},\qquad 0<\epsilon<1. \tag{3}
\]

**The clock-rate ratio is a standard definition; its dependence on W is a provisional hypothesis using a standard suppression function.** The reference time t and physical clock-comparison procedure must be specified. W_* is positive and has the same units as W. Deep wells approach q=1; an ideal void approaches q=1-epsilon.

\[
\alpha(W)=\kappa[1-q(W)]
=\frac{\alpha_{\rm void}}{1+W/W_*},
\qquad \alpha_{\rm void}=\kappa\epsilon. \tag{4}
\]

**Proposed interaction law, not a deduction from clock slowing.** Alpha and kappa have units of inverse length. A stationary environment is permitted to generate continuing conversion. Redshift generally constrains kappa times epsilon, not the clock alteration and interaction efficiency separately. No values of these parameters are locked in by this document.

This q is a physical hypothesis, not merely a relabeling of time coordinates. The earlier universal-clock Hamiltonian candidate instead generated cumulative matched-endpoint shifts from temporal change of q. Its conditional relation was

\[
\alpha_{\rm time}=-\frac{\partial_t q}{c q^2}. \tag{5}
\]

**Established Hamiltonian calculation within that separate candidate.** Equation (5) and equation (4) must not be imposed together without a compatible interaction. In that earlier coupling, spatial gradients also produced an additional force; slower void clocks were not a free change with no mechanical consequences. See [universal-clock calculation](../research_work/results/universal-clock-coupling/report.md). A common matter, photon and companion theory must settle this compatibility.

### 3. Photon energy, wavelength and cumulative redshift

\[
\frac{dE_\gamma}{ds}=-\alpha(W)E_\gamma,
\qquad
A=\int_{\rm traveled\ path}\alpha(W)\,ds. \tag{6}
\]

\[
E_{\gamma,\rm out}=E_{\gamma,\rm in}e^{-A},
\qquad
1+z_{\rm transfer}=e^A,
\qquad E_\gamma=h\nu=\frac{hc}{\lambda}. \tag{7}
\]

**Equation (6) is the proposed fractional-transfer rule; equation (7) uses its standard exponential solution and established photon relations, retained as assumptions in this universe.** It describes surviving photons losing energy, not just photon removal. The latter could dim a beam without reddening its surviving photons.

The spectral transfer factor is not automatically the entire observed redshift. Motion, endpoint clock comparisons and measurement conventions must be treated consistently; do not add competing causes twice.

For a constant environment, ideal stationary mirrors and no resetting at reflection:

\[
A=\alpha L_{\rm total}=\alpha cT. \tag{8}
\]

**Conditional consequence of the rule and constant local c.** T is local elapsed time in that stationary region. Total path length, not displacement, counts. A traveler can still feel normal while detecting a difference between stored and freshly emitted light. A full cavity theory must distinguish changing photon frequency from ordinary cavity energy loss and must respect its mode structure.

### 4. Traveling companions and spatial conservation

Before capture, for an isolated packet with no other exchanges:

\[
\frac{dE_c}{ds}=+\alpha E_\gamma,
\qquad E_\gamma+E_c=E_{\rm initial}. \tag{9}
\]

**Standard energy accounting under the proposed complete-transfer assumption.** Both components move forward at local c in the initial idealization. This statement does not prove their coherence, entanglement or identity as gravitons.

For an effective fixed-space transport calculation using one common energy and time convention:

\[
\begin{aligned}
\partial_tu_\gamma+\nabla\cdot\mathbf F_\gamma&=j_\star-Q,\\
\partial_tu_c+\nabla\cdot\mathbf F_c&=Q-C,\\
\partial_tu_d+\nabla\cdot\mathbf F_d&=C,
\end{aligned}
\qquad
Q=c\alpha u_\gamma,\quad C=\Gamma_{\rm cap}u_c. \tag{10}
\]

**Established continuity equations with proposed transfer and capture closures.** u denotes energy per volume, F energy flux, and j_star injected stellar photon power per volume. Gamma_cap has units of inverse time. F_d allows bound energy to redistribute; a static-store approximation sets it to zero. For multiple incoming directions, transport must retain angular information; the magnitude of the net flux is not c times the total energy density when beams cancel.

Adding the three equations gives:

\[
\frac{d}{dt}\int_V(u_\gamma+u_c+u_d)\,dV
=\int_V j_\star\,dV
-\oint_{\partial V}(\mathbf F_\gamma+\mathbf F_c+\mathbf F_d)\cdot d\mathbf A. \tag{11}
\]

**Derived conservation identity for this effective ledger.** Energy leaving a volume is not destroyed. Include source fuel for a closed total budget. Local and reference energies cannot be interchanged in a varying q field without deriving appropriate work and lapse factors. A full gravitational theory requires its own consistent conservation treatment. See [spatial transfer calculation](../research_work/results/transport/spatial-transfer.md).

### 5. Capture becomes available in deep wells

One possible smooth-onset prescription is:

\[
\Gamma_{\rm cap}(W)=
\begin{cases}
0,&W\leq W_c,\\
\Gamma_{\max}\left[1-e^{-(W-W_c)/\Delta W}\right],&W>W_c.
\end{cases} \tag{12}
\]

**Proposed illustrative rate, not a selected or derived microscopic law.** W_c is the threshold, Delta W its transition width, and Gamma_max a limiting rate. The expression implements zero capture in voids when their W lies below threshold. It does not itself explain binding, recoil, permanence or the final orbit of a captured companion.

The physical proposal is traveling excitation -> interaction with a well -> bound excitation plus the appropriate response of the surroundings. Wave-particle duality alone does not provide this transition. Calling an excitation a particle does not give it rest mass or stop its motion.

### 6. Bound levels as extended energy-storage patterns

For a simple description with approximately fixed, noninterfering or phase-averaged patterns:

\[
u_d(\mathbf{x})=\sum_n U_n f_n(\mathbf{x}),
\qquad f_n\geq0,\qquad\int f_n\,d^3x=1. \tag{13}
\]

**Standard normalized spatial decomposition, proposed as an effective model of stored companion energy.** U_n is the total energy in pattern n; f_n has units of inverse volume. This representation does not prove discrete quantum levels. Coherence, interactions and self-gravity can require additional terms and changing patterns.

Each allowed bound configuration would have a spatial extent, an excitation energy, and transition/capture rules derived from a field equation. Populating extended configurations could create extra gravity over a broad range of radii. It need not create a thin outer shell or make density increase outward. Nor do bosonic populations automatically fill like electron shells.

The closest literature analogy is hypothetical massive boson clouds around black holes, whose spectra can resemble atomic levels. Those models are not demonstrations of our photon-supplied, initially light-speed companion mechanism or literal quantized layers of spacetime. [LIGO research on gravitational-atom-like spectra](https://dcc-llo.ligo.org/public/0178/P2100343/008/O3_boson_cloud.pdf).

Long-lived extended patterns must be derived rather than presumed: rapid inward relaxation, leakage or decay would change the proposed halo. The repository already contains conditional [self-binding](../research_work/results/companion-self-binding/report.md), [bound-cloud exchange](../research_work/results/bound-cloud-exchange/report.md) and [pair-production balance](../research_work/results/pair-production-balance/report.md) studies. They constrain candidate mechanisms but do not establish this combined model.

### 7. Conditional gravitational response and the imposed-map benchmark

As a baseline only, suppose deposits behave gravitationally like nonrelativistic matter with small stress:

\[
\rho_d=\frac{u_d}{c^2},\qquad
\nabla^2\Phi_{\rm grav}=4\pi G(\rho_{\rm ordinary}+\rho_d). \tag{14}
\]

**Established weak-field/Newtonian mathematics under a proposed deposit equation of state and standard coupling.** This is not generally valid for arbitrary relativistic waves or arbitrary modified gravity. A companion-field model must determine pressure, momentum and gravitational response. [Poisson equation](https://farside.ph.utexas.edu/teaching/355/Surveyhtml/node50.html).

For spherical Newtonian circular orbits:

\[
v_{\rm circ}^2(r)=\frac{GM(<r)}{r}. \tag{15}
\]

**Established conditional relationship.** Approximately flat speeds require enclosed gravitating mass to increase approximately linearly with radius over the relevant range. They do not demand a density pileup at an edge. Lensing must follow from the same gravitational theory; pressure and the spatial metric can distinguish lensing from stellar dynamics. [Joint rotation/lensing analysis](https://academic.oup.com/mnras/article/372/1/136/971997).

The user's temporary distribution benchmark is:

\[
u_d(\mathbf{x},t_{\rm now})=c^2\rho_{\rm target}(\mathbf{x}). \tag{16}
\]

**Imposed target under the assumptions of equation (14), not a prediction.** Rho_target is a selected conventionally inferred dark-matter density map, reused as a fictional target with its inference assumptions labeled. Location alone is insufficient: normalization and stress/gravitational behavior matter. No particular target map has yet been chosen or loaded for this specification. Later, equation (13) and physical capture dynamics should predict a map to compare against this target rather than inherit it.

## What the combined model does and does not yet predict

It specifies a conditional cumulative redshift for any supplied W path and parameters; a transparent transfer ledger; no storage in voids; and a possible route from traveling energy to extended bound patterns. These are calculable assumptions, not a validated first-principles explanation.

In the tested stationary constant-speed packet model, photon energies change but the time spacing of separate emitted packets does not. Stretching a supernova's entire light curve therefore remains an independent requirement. A mode capture sketch also does not establish narrow spectral lines, sharp images, a thermal background, lensing strengths or a viable clock/force coupling.

No fit, calibrated clock slowing, capture threshold, bound-level spectrum, stable halo, or new observational success is claimed by this document. The earlier PDF remains a separate version; this specification does not silently update its contents.

## Next executable research steps

1. Specify an ordinary-matter W map, reference convention, and treatment of deposited-gravity feedback. Outcome: the same reproducible environmental input for every object, not a field inferred from each object's redshift residual.
2. Construct a candidate companion field with a freely propagating branch and bound configurations. Outcome: actual allowed spatial patterns and a lifetime calculation, not only a chosen halo shape.
3. Derive an interaction that links the proposed time environment, photon frequency/phase, traveling companions and capture. Outcome: conservation of energy and momentum, with clock, wave and recoil behavior from one model.
4. Predict which bound patterns are populated and whether they stay extended. Outcome: a radial distribution and its evolution to compare with the optional imposed target.
5. Calculate stellar motion and lensing from the same field response, then test redshift, brightness, event duration, spectral and image preservation, and joint light/gravitational-wave propagation. Outcome: shared-parameter predictions with reported failures.
6. Revisit the deferred energy supply and full formation history after the mechanism and spatial behavior are specified. Deferral is not evidence that this requirement has passed.

Supernovae and gravitational-wave standard sirens provide complementary observational tests; they do not directly measure W everywhere along a trajectory. Imported distance and waveform assumptions must be identified, especially when photon or gravitational-wave propagation is modified. [Standard sirens](https://ligo.org/science-summaries/GW170817Hubble/).
