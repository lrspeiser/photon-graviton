# JR-3: does gravity-related photon exchange cause galaxy dimness?

Exploratory calculation, 21 September 2026 America/Los_Angeles. Baseline main 41ede5342a0ddabdff7dbea2a4f1951eceec65c6. The user proposes reversed causality: low observed light may result from conversion into a gravitational companion rather than imply less matter. Preserve the hypothesis and all observations.

## Separate physical questions

1. Hidden starlight: uniform band transmission T implies an intrinsic stellar source correction u=1/T at fixed intrinsic population M/L. Propagate that correction into ordinary gravity AND the frozen R10 source-linked companion. This is not a derived microscopic attenuation law, and it does not add a second copy of the R10 gravitational reservoir.
2. Repeat with equal attenuation of the gas flux as an explicitly achromatic sensitivity. It is not assumed that infrared and 21-cm photons share a rate.
3. Reversible exchange: solve a radiation/companion energy system with forward, reverse, photon escape and companion escape. Determine when a stationary stored state can coexist with lower outgoing bolometric luminosity.
4. Supply audit: if the existing R10 effective reservoir gravitates as ordinary stored energy, derive the missing radiative energy needed to supply it. Do not assume a galaxy or universe age, or substitute a 3.6-micron solar-unit luminosity for bolometric watts. Carry the unknown bolometric conversion B, capture efficiency eta and effective historical duration T explicitly.

## Frozen inputs and data scope

149 SPARC galaxies and six SLACS lenses, using the supplied JR1/JR2 packages. No new cluster forward model exists. Original radii, inclinations, gas maps, light shapes, velocities, covariance and lens angles stay fixed. All data have already been inspected in the project. No independent-confirmation claim. R10 is a phenomenological control, not the overall companion ontology.

## Disk calculations

Compare baseline, 25% and 50% net band dimming, independent per-galaxy optimal net dimming T in [0.1,1], and a signed gain/loss diagnostic u in [0.1,10]. Use both original quoted-error chi-square and fractional RMS; the primary shape/amplitude diagnostic minimizes fractional RMS and retains raw chi-square. No fitted noise floor. Free per-galaxy T is a feasibility/identifiability calculation, not a causal universal law or a measured transmission.

Fit two shared tau laws using the archived 89-galaxy training split only: constant log u=tau0, tau0 in [0,ln10]; and gravity-linked log u=tau0 H(u), with H equal to the modeled companion fraction of positive total squared circular speed at Re, evaluated self-consistently on the hidden-light-corrected source. Same tau0 in [0,ln10]. H never uses an observed rotation speed. It is not independent of the prior R10 fit. Solve all roots in a fixed 65-cell bracket and use the lowest branch; report multiple roots if encountered. Evaluate existing validation/comparison splits once after the one-parameter training fit; these are exposed reuses.

Record per-object error signs and parameter-bound hits, including all 54 prior outliers rather than only improved objects. No emission spectra are inferred.

## Lenses

At unchanged distances and angular profiles, solve the attenuation required for each ring by changing the intrinsic stellar normalization consistently. Compute the resulting stellar motions at frozen R10 orbital parameters and also with one orbital nuisance parameter refitted in [-1,.35]. Retain lenses that require net brightening rather than dimming. This is not a multiband population analysis, a full image fit or a common opacity-law derivation.

## Reversible energy model

    E_gamma_dot = L_in - k_escape E_gamma - k_plus E_gamma + k_minus E_chi
    E_chi_dot = k_plus E_gamma - (k_minus+k_chi_escape) E_chi
    L_out = k_escape E_gamma.

All rates are nonnegative. Stored plus radiated/escaped energy must equal injected energy. Test positive rates and compare analytic stationary solution, a direct linear solve and time evolution from zero stored energy. Distinguish bolometric total escape from a specific band, direction or observer aperture. The model is a rate-limit example, not a restriction on coherent/nonlocal exchange, formation suppression, or all possible gravitational mechanisms.

## Energy interpretation

Let Lcat be the dimensionless catalog 3.6-micron solar-unit luminosity and set L_bol,observed = B * Lcat * L_sun,bol by DEFINITION of unknown B. Then missing net luminosity is B*(u-1)*Lcat*L_sun,bol in an explicitly grey/constant-history illustration. Existing reservoir Mchi=A*rt/G requires

    B*eta*T >= Mchi*c^2/[(u-1)*Lcat*L_sun,bol].

Report also the finite enclosed mass at the last measured radius rather than only the extrapolated total. T is an effective historical integral, not an assumed actual age; no B or eta is measured. Returning the same energy cannot be counted as newly accumulated energy on every cycle. Modified response not proportional to E/c^2 needs its own energy-to-gravity relation and is not ruled out by this conditional inventory calculation.
