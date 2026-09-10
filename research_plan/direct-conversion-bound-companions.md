# Active branch: Direct Photon Conversion and Bound Companion Fields

9 September 2026. Selected after the user's concern about the time-sign reversal and screening. This supersedes additional environmental clock factors and conversion screening as the active starting mechanism. Those earlier hypotheses and their tests remain available as historical alternatives, not ingredients silently retained in this branch.

## Physical proposal

Light gradually transfers energy into companions throughout its traveled path. There is no minimum path length at which conversion begins: short paths accumulate little loss. The initial model has one constant fractional-loss coefficient. Local clocks and the local speed of light retain ordinary behavior; standard gravitational and motion-related spectral shifts remain part of observation modeling. No extra q field is invoked to create the conversion.

Companions initially travel forward at c and retain their energy while crossing voids. They do not permanently collect there. An interaction in deeper gravity wells may transfer traveling energy into long-lived, spatially extended bound configurations. These configurations contribute gravity and are to be tested against galaxy rotation and lensing. The idea of bound modes is retained as a candidate; neither a mode spectrum nor a capture interaction has been established.

## What changes and what stays

| Decision | Active treatment |
|---|---|
| Nonexpanding universe with stipulated published distances | Retained. Catalog definitions and inference assumptions remain labeled. |
| Additional environmental time slowing or speeding | Not required in this branch; retired from the active conversion cause. |
| Screening of photon conversion in gravity wells | Not included in the initial model. Test actual local predictions rather than imposing it by default. |
| Long-distance photon energy loss | Retained as continuous fractional transfer, starting from zero path length. |
| Constant locally measured c and ordinary local clocks | Retained. |
| Energy transferred to forward companions traveling at c | Retained as an effective propagation assumption. |
| Further energy loss or permanent capture in voids | Absent by assumption. |
| Capture in deeper wells | Retained as an unresolved interaction. This is a storage condition, not suppression of photon conversion. |
| Bound levels/configurations | Retained as a candidate way to distribute stored energy. |
| Ordinary gravity | Initial baseline for matter and lensing; deposited-field stresses must be calculated. Modified response remains an alternative requiring derivation. |
| Conventional dark-matter map | May be an explicitly imposed distribution benchmark, never evidence that capture has produced that map. |
| Source-energy budget | Deferred for now, not declared adequate. Conservation is still required. |
| Comparison with expansion | Comparable explanatory and predictive performance is sufficient; superiority is not required. |

## Current equations and provenance

### Continuous conversion

\[
\frac{dE_\gamma}{ds}=-\alpha_0E_\gamma,
\qquad
\frac{dE_c}{ds}=+\alpha_0E_\gamma.
\]

**Proposed transfer postulate written as standard rate equations.** The second equation describes the isolated generated companion packet before capture and excludes additional exchanges. Alpha_0 is constant for this first branch. It must eventually come from a physical interaction; specifying it alone does not explain momentum, photon number preservation, phase or spectral linewidth.

\[
E_\gamma(D)=E_{\gamma,0}e^{-\alpha_0D},\qquad
1+z_{\rm transfer}=e^{\alpha_0D},\qquad
\frac{E_c}{E_{\gamma,0}}=1-e^{-\alpha_0D}.
\]

**Established exponential solution, using the usual photon energy-frequency relationship.** The formula is not novel and is of the familiar fractional-energy-loss/tired-light form. The proposed physical carrier, interaction, capture and joint predictions would need to supply any novel content; their originality is unverified. Alpha_0=0.0002488993286382367/Mpc is reused from an exposed-data fit, not derived or freshly validated.

### Spatial transport and storage

\[
\begin{aligned}
\partial_tu_\gamma+\nabla\cdot\mathbf F_\gamma&=j_\star-Q,\\
\partial_tu_c+\nabla\cdot\mathbf F_c&=Q-C,\\
\partial_tu_d+\nabla\cdot\mathbf F_d&=C,
\end{aligned}
\qquad Q=c\alpha_0u_\gamma.
\]

**Established continuity structure with proposed source terms, in a declared common energy/time convention.** Multiple incoming directions require angular transport. Internal capture terms cancel from the sum; source fuel, exported energy, recoil and gravitational/field work must be included where relevant. These scalar equations alone are not a complete relativistic conservation law.

\[
C=\Gamma_{\rm cap}(W,\text{local state})u_c,\qquad
W=\Phi_{\rm ref}-\Phi_{\rm env}.
\]

**A possible effective capture closure and standard potential-difference definition.** Gamma_cap is not chosen or derived. Require no capture in voids in the initial idealization. The environmental potential, its reference and possible feedback from stored energy remain to be specified. They no longer control photon conversion through alpha_0.

\[
u_d(\mathbf x)=\sum_n U_n f_n(\mathbf x),\qquad
\int f_n(\mathbf x)d^3x=1.
\]

**Standard decomposition, proposed as a phase-averaged/fixed-pattern approximation.** U_n is energy in pattern n. Derive allowed patterns, their occupations, stability, self-gravity and interactions before treating this as a predicted halo. Do not assume that being called a particle causes light-speed excitations to stop or bind.

### Conditional gravity

\[
\rho_d=u_d/c^2,\qquad
\nabla^2\Phi=4\pi G(\rho_{\rm ordinary}+\rho_d).
\]

**Established Newtonian weak-field baseline under a proposed effectively nonrelativistic, low-stress deposit state.** It is not automatically the right source law for relativistic waves. The deposited field's stress-energy and corresponding metric must predict both matter motion and light bending. Standard-gravity identities are conditional comparisons, not restrictions on every allowed modified-gravity branch.

### Whole-event timing and brightness

For the particular stationary fixed-speed closure:

\[
t_{\rm arrival}=t_{\rm emission}+D/c,
\qquad
\frac{\Delta t_{\rm arrival}}{\Delta t_{\rm emission}}=1.
\]

**Established travel-time calculation under the declared idealization.** Photon frequency loss does not automatically stretch a supernova's full light curve. This is a known deficiency of this minimal closure, not something to hide by importing a factor of 1+z without a mechanism.

For isotropic emission, photon-number retention, a static geometry and no extra event stretch:

\[
F_{\rm bol}=\frac{L_{\rm bol}}{4\pi D^2(1+z_{\rm transfer})},\qquad
D_L=D\sqrt{1+z_{\rm transfer}}.
\]

**Conditional consequence of the photon-loss law, inverse-square spreading and standard luminosity-distance definition.** These are bolometric expressions; bandpass, spectral changes, dust and surface-brightness-fluctuation calibration require their own forward model. The apparent distance inferred from brightness need not equal physical path length. Adopted catalog numbers remain unchanged; their interpretation must be explicit in each comparison.

## First integrated test and current status

The [direct-conversion test report](../research_work/results/direct-conversion/report.md) provides a 100-million-light-year prediction, five distance-selected group examples, all 164 predictions, and explicit expansion comparisons with matched distance conventions. It finds no demonstrated improvement over the expansion benchmark. The [screening diagnostic](../research_work/results/void-screening/report.md) explains why small local photon loss does not require a large clock modification.

The minimal model can produce an empirical redshift-distance trend while maintaining a scalar transfer ledger. It does not yet explain the microscopic loss, the observed event-duration behavior, capture, stable extended storage, lensing or the source budget. A useful surviving model must complete these links without per-object tuning.

## Goals and outcomes

1. **Redshift benchmark:** freeze this one-parameter starting rule; audit independent distance and redshift provenance, motion and grouped uncertainties; predict fresh objects. Outcome: actual measured versus predicted redshifts and uncertainties, not a fitted value for one selected galaxy.
2. **Interaction and wave propagation:** derive a number-retaining energy transfer with energy/momentum/field accounting, narrow lines and image preservation; specify its frame and symmetries. Outcome: a common wave/transport equation and calculable alpha rather than a name for the observed slope.
3. **Event timing and flux:** determine whether the same interaction can explain time-stretched transient profiles and their brightness, or document failure. Outcome: a joint prediction of spectrum, duration and flux. Do not reintroduce a sign-changing clock factor merely to force agreement.
4. **Capture and bound patterns:** derive conversion of freely traveling companions to supported long-lived extended states. Outcome: a predicted distribution with a stability/escape/transition calculation.
5. **Gravity:** predict galaxy rotation, cluster lensing and merging-cluster behavior using that distribution and the same gravitational response. Outcome: actual observables, not a halo profile imposed from the desired answer.
6. **Formation and energy supply:** return to the deferred budget and source histories. Outcome: a conserved, causal formation history with sufficient supply under the selected gravitational law, or a recorded failure.

Expansion remains a comparison theory with empirical tests, not an assumed explanation in this fictional universe. The user clarified that comparable performance is sufficient: superiority is not required. Establish comparability with predeclared tolerances, calibrated uncertainties and independent observations; a nonsignificant difference on an exposed sample does not prove equivalence. The [prior-art and parity note](../research_work/results/direct-conversion/prior-art-and-parity.md) records this clarification. Earlier work and all 32 observational requirements remain preserved in the [research roadmap](START-HERE.md).
