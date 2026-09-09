# Calibrating the inverse-well time effect with observed redshift

Local interpretation, latest user clarification: a traveler experiences normal clocks and local processes throughout the trip. The effective gearing below describes a comparison of signals across the journey; it is not a slowdown visible on the traveler's own watch. The eventual field/matter model must recover both local normality and the specified nonlocal observable stretch.

User direction, 9 September 2026: use the observed redshift to measure the impact of the proposed environmental time mechanism. Adopt this as an inverse problem: observations set the stretching target, while independently measured environments and a common physical law determine whether that target can be explained predictively.

## What observations determine directly

Established definition and mathematical transformation; not novel:

S_obs = 1 + z_obs, T_obs = ln(1 + z_obs).

S_obs is the received-to-emitted wavelength ratio under the stated spectral/reference-frame convention. T_obs expresses the same measurement on an additive scale, convenient for combining successive effects. It is dimensionless, not a number of seconds and not a local clock-rate measurement.

Established multiplicative factor bookkeeping, applied to this hypothesis:

1 + z_obs = (1 + z_time)(1 + z_motion)(1 + z_endpoint)(1 + z_measurement).

Therefore T_time = ln(1 + z_obs) - ln(1 + z_motion) - ln(1 + z_endpoint) - ln(1 + z_measurement).

This factorization requires consistent frame definitions and no double counting of endpoint/motion effects already represented by a propagation model or a catalog correction. Measurement factors describe an observation model, not a new physical energy channel. Independent motion constraints remain missing; no individual velocity is fitted to cancel a residual.

The companion CSV transforms all 164 already exposed observations without changing their CMB shifts or distances. Its all_shift_time columns are conditional targets assuming the entire observed shift comes from the proposed time mechanism. They are not corrected estimates or fresh predictions.

## Connecting observations to the inverse-well function

Optional joint parameterization proposed here; originality unverified. The saturating shape is generic mathematics, not a novelty claim:

I(W) = 1/[1 + (W/W_star)^2], alpha_time(s,t) = A I(W(s,t)), A >= 0, W_star > 0.

W measures physically defined nonnegative well depth along the ray. A is the maximum fractional stretch rate per unit distance. W_star is the well depth at half strength. This is a candidate phenomenological link between environment and stretch, not yet a derivation from a physical time field.

Established integration of fractional wavelength change, conditional on that proposed rate:

T_time = integral_path alpha_time ds = A integral_path I(W(s,t)) ds = A D_eff,

z_time = exp(A D_eff) - 1.

D_eff is the independently determined path length weighted by the environmental effect. It is not an additional geometric distance or expansion. Large redshift with large D_eff can calibrate A; variation among independently measured environments can constrain W_star. Use one common parameter set and training-only calibration, then predict untouched objects. The exponential solution itself is known mathematics.

If every path is deep in the void limit, D_eff is approximately the total distance and W_star cannot be determined from those data. If every path is deep-well dominated, only the product A W_star^2 is approximately identifiable for known W. Mixed paths spanning the transition are needed to distinguish amplitude and transition scale. These are parameter degeneracies derived from the proposed function, not observations establishing either regime.

No independent W(s,t) paths are currently supplied. Therefore no environmental parameter fit is performed. Defining W from each target's redshift residual would force agreement by construction and defeat the requested predictive test.

## Why this does not yet measure local clock slowing

### Mostly-void effective gearing, latest user direction

Assume a fraction F of the stipulated distance D is fully affected void, with negligible added effect on the rest. F=1 is the all-void approximation; F=0.9 and 0.75 in the table are sensitivity scenarios, not measured fractions. This two-region limit is an optional approximation to the smooth I(W) weighting.

Derived conditional calibration using established fractional-change mathematics; no novelty claim:

A_required = ln(1 + z_time)/(F D).

Effective signal-clock ratio, defined here as a reporting convention; not a local physical-clock law:

q_signal,eff = nu_obs/nu_emit = 1/(1 + z_time),

lambda_obs/lambda_emit = 1/q_signal,eff.

In this equation the frequencies and wavelengths refer to the isolated time-mechanism component with declared endpoint standards. For example, z_time=0.1 requires a total wavelength multiplier 1.1 and a frequency multiplier 0.90909. We can describe this as 9.09% effective signal slowing. It does not establish that clocks ran 9.09% slower throughout the void. The same total stretch can result from many different histories.

Under the constant void-rate approximation, the remaining signal-frequency ratio after affected path length ell is exp(-A_required ell); the wavelength multiplier is exp(A_required ell). These are derived consequences of the proposed constant fractional stretching law using established exponential mathematics. At ell=F D they reproduce the target by construction. They give the time-field derivation a specified pathwise target, not an independent test of it. Keep the distinction between this signal-evolution history and q(W), the physical local-clock ratio.

The illustrative F=0.9 scenario raises the inferred void stretching rate by 11.11% relative to F=1, without changing the observed total gearing. F=0.75 raises it by 33.33%. Redshift alone cannot select F, W_star or the physical local-clock law. Fit a shared rate across the calibration sample rather than allowing each object's A_required to become a free prediction parameter. These per-object inversions expose what the common model must explain and where it fails.

The earlier optional q(W) = 1 - epsilon I(W) clock prescription (proposed here; originality unverified) is separate from alpha_time = A I(W). A stretching rate per distance is not the same quantity as a fractional clock change. We have not derived a relationship between A and epsilon. Static clock differences alone do not establish lasting cumulative stretching between matched endpoints under the previously audited stationary propagation assumptions.

The next physical model must supply that relationship through time evolution or an explicit history-dependent interaction, including atomic clocks, local light speed, arrival-time spacing and wavelength standards. Until then, redshift calibrates integrated stretching, not a unique amount by which time itself slowed. The user has fixed the environmental direction; a physically realized sign and magnitude of local clock change remain unresolved.

Established photon energy relation and conditional algebra, not a novel law:

Delta E_companion/E_emit = 1 - 1/(1 + z_time) = z_time/(1 + z_time).

This follows if photon energy is proportional to frequency and all loss assigned to the time mechanism transfers into companions in the same energy convention. It does not determine the deposited gravitational response. Using z_obs in its place assumes every observed shift is attributable to this mechanism; the table labels that assumption explicitly.

## Research decision

Observed redshift is now the calibration target for the inverse-well premise. Keep the frozen constant-rate and unsuccessful smooth-rate comparisons. Acquire independent environmental path inputs, select a small physically motivated extension, calibrate only on exposed training data, and freeze its parameters and uncertainty before fresh evaluation. No new precision, clock-slowing measurement or fundamental explanation is claimed by converting observed redshifts into stretch targets.
