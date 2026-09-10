# Distinguishing mechanisms with timing, spectra and messenger arrivals

**Status: exploratory diagnostic, 2026-09-10.** No physical parameters fitted, fresh holdouts opened or unified theory validated. [Protocol](protocol.json), [calculation](run.py), [numerical results and input hashes](results.json). The familiar equations below are labeled by origin; the proposed physical interpretation has no established originality claim.

Similar-looking redshift formulas can represent very different physics. Ask whether a proposed mechanism also stretches an event, preserves the relative positions of spectral lines, and delays light relative to gravitational waves. Requiring one mechanism to answer all three is more discriminating than fitting another redshift curve.

## Candidate comparison

| Specified candidate | Event duration | Frequency behavior | Relative light/GW delay | Current interpretation |
|---|---|---|---|---|
| Achromatic photon-energy transfer with unchanged, fixed-speed transport | Unchanged | Common fractional shift by postulate | No additional flight delay | Poor exposed supernova-aging score; energy bookkeeping alone does not give event stretching |
| Photon-only inverse-affine time transport, uniform positive coupling, ordinary endpoint clocks | Multiplied by the spectral factor S | Achromatic in the assumed wave Hamiltonian | At least about 675,000 years for the stipulated GW170817 distance and fixed rate | Incompatible with the associated merger signals for this specialization and ordinary merger-scale emission offsets |
| The same inverse-affine transport for photons and gravitational waves | Multiplied by S | Achromatic under the same dispersion assumptions | Common flight delay cancels; observed source lag is multiplied by S | Not rejected by this relative-delay check; clock physics, GW waveform/amplitude and companion coupling remain unspecified |
| Spatially void-dependent inverse-affine coupling | Multiplied by exp(integrated coupling) | Achromatic if coupling is frequency independent | Depends on where conversion occurs, not just its total | Needs a measured environmental profile and a common rule; arbitrary profiles can evade a universal delay bound |
| Coherent quantum receiver map from earlier audits | Conditional Fourier-profile stretch | Conditional spectral rescaling with receiver-dependent broadening | No physical source-to-detector flight law supplied | Not yet scoreable as a complete propagation model |

These rows are physical specializations, not ideas attributed to historical scientists. A missing prediction is not a successful prediction. Giving photons and GWs the same transport law is a new candidate, not a retrospectively fitted repair counted as validation.

## A new conditional messenger calculation

**Proposed transport law**, using ordinary endpoint clock coordinates and static path length s:

\[
\frac{dt}{ds}=\frac1c+\kappa(s)(t-t_*),\qquad \kappa\geq0.
\]

Here t-star is a physical reference epoch, not a removable coordinate convention once the law is chosen. Nonnegative field age at emission is assumed. The photon speed in these coordinates is the inverse of the right side; a complete theory must still explain rods and clocks that make locally measured light speed c. This calculation does not supply that completion.

**Known integrating-factor mathematics, conditional on that postulate:** define

\[
K(s)=\int_0^s\kappa(u)\,du,\qquad S=e^{K(D)}.
\]

Neighboring emission times then have arrival separation multiplied by S. With the additional proposed nondispersive Hamiltonian omega=v(t,s)k, Hamilton's equations give d ln omega/ds=-kappa, so the same S stretches wavelengths. An energy-loss equation alone does not establish these Hamiltonian or clock assumptions. The exponential structure is established mathematics, not a new formula.

Relative to an unaffected GW traveling at c on the same path, **conditional deduction**:

\[
\Delta t_{\rm prop}=(S-1)(t_e-t_*)+
\frac1c\int_0^D\{e^{K(D)-K(s)}-1\}\,ds.
\]

For constant kappa=alpha and the most favorable age t_e=t-star:

\[
\Delta t_{\rm prop,min}=\frac{e^{\alpha D}-1}{\alpha c}-\frac D c.
\]

**Prior empirical calibration held fixed:** alpha=0.0002488993286382367 per Mpc. This is not derived from first principles. The measured host distance used here is 40.7 +/-1.4 random +/-1.9 systematic Mpc from surface-brightness fluctuations, not a Hubble-law distance. We stipulate it as a path distance in the fictional static model. Its brightness calibration must be revisited if the model changes the inference of distance. [Cantiello et al.](https://arxiv.org/abs/1801.06080)

The calculation gives S=1.0101816869 and a minimum extra delay of **674,646 years**. Distances 30 and 50 Mpc are also evaluated in results.json as illustrative sensitivity checks, not a distance posterior. This is a prediction for the transfer contribution, not a fit to the host's total measured redshift or its peculiar motion.

The gamma-ray burst arrived **1.74 +/-0.05 seconds after** GW170817. The measurement includes intrinsic source timing, so we do not equate 1.74 seconds with propagation delay or import the paper's conditional speed bound into our cosmology. [LIGO/Virgo/Fermi/INTEGRAL](https://arxiv.org/abs/1710.05834)

If gamma emission follows GW emission by delta_emit, the conditional relation is

\[
\Delta t_{\rm observed}=S\,\delta_{\rm emit}+\Delta t_{\rm prop}.
\]

Formally cancelling the minimum uniform delay requires gamma emission about **667,847 years before the merger**. This algebraic offset is not a physically acceptable emission model for the associated event; enforcing a nonnegative field age for that earlier emission only tightens the constraint. The conclusion depends on the stated law, uniform positive coupling, ordinary endpoint clocks and unaffected GWs. It is not an exclusion of every nonexpanding universe or every energy-transfer mechanism.

For two messengers following the same affine arrival map, the common propagation intercept cancels. Their lag is S times their intrinsic emission separation, allowing a merger-scale offset. A universal field must subsequently be checked against matter clocks and the GW waveform, not merely declared successful. If companions are ordinary gravitons and share a frequency-loss law too, their requested loss-free propagation needs separate justification; equal travel speed alone does not force equal energy exchange.

## Why void placement matters

**Conditional deduction, same law:** an ideal thin interaction at distance d from the source, with fixed total K, gives minimum extra delay (S-1)d/c. A uniform slab beginning at d0 with width w gives

\[
\Delta t_{\rm prop,min}=\frac{(S-1)d_0}{c}
+\frac{w}{c}\left[\frac{e^K-1}{K}-1\right].
\]

Concentrating all conversion arbitrarily close to the source makes the delay arbitrarily small. Therefore S and D alone do not establish a positive profile-independent lower bound. Conversely, distributing the interaction across ordinary intergalactic distances does not inherit this loophole automatically. Four synthetic slab placements with identical integrated conversion are retained in results.json; none represents a measured void.

As illustrative diagnostics only, assume an absolute intrinsic source lag no greater than 1, 10 or 1000 seconds and allow the quoted observed lag plus one quoted error. The corresponding maximum source distances of the thin interaction are approximately 0.55, 2.34 and 199 AU. These are conditional scale comparisons, not inferred confidence bounds or a suggested astrophysical mechanism. A void law must predict its spatial placement before seeing the target lag.

## Reused observations that already distinguish mechanisms

- All 35 previously exposed supernova spectral-aging rows give diagonal chi-square 150.569 for unchanged event duration, versus 26.949 for duration stretched by 1+z. Shared-template uncertainties, evolution and missing covariance prevent treating these as a complete likelihood. This reproduces earlier evidence, not a new independent result.
- The existing two independent methanol centroid groups test whether two frequencies acquire the same fractional shift. Previously explored order-unity power dependences predict differential velocities around 35,000 to 66,000 km/s in magnitude, whereas the published shared-gas offset is 0.66 +/-0.39 km/s. The p=-1 case gives approximately -45,285 km/s. The achromatic p=0 case predicts zero offset. Line profiles, gas association and calibration remain part of the interpretation; this is not an extreme-tail significance claim. The integrated loss is anchored to one absorber line because no independent path distance was used. [Methanol measurement](https://arxiv.org/abs/1412.7757)

We do not sum these heterogeneous diagnostics into a global score. The aging calculation conditions on observed redshift: it tests a timing relation, not the independently predicted redshift-distance relation.

## Verification and what follows

Independent numerical ODE propagation agrees with the analytic slab formula within 1.12e-15 in dimensionless arrival time. Finite emission intervals reproduce S, zero coupling removes delay, positive field age increases delay, and shared-messenger propagation cancels its common intercept. These verify the conditional calculation, not its physical assumptions.

The new four-stage research objective remains active. Stage 1 now has a concrete discrimination result; incomplete candidates still lack predictions. Stage 2 requires independently determined path environments. Stage 3 requires shared predictions for spectral shape, clocks and brightness. Stage 4 requires a frozen model and genuinely unopened outcomes. See the [forward test contract](../../../research_plan/joint-observation-test-contract.md). No gravity parameters, capture model, energy-supply status or reserved stellar outcomes changed.
