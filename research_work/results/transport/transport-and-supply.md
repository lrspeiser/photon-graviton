# Energy, escape and the next theory choices

Status: conditional calculations, not an established interaction or a fitted galaxy theory. Both a larger/older source history and a different gravitational response remain open. No new physical law has been adopted.

## What the calculation establishes

The proposed no-loss rule conserves energy within its stated three-reservoir model. To apply it to one galaxy, however, we must also count energy crossing the galaxy boundary. Keeping energy during flight does not by itself keep a companion inside that galaxy.

The extended solver passed nine checks: independent numerical integration, conservation through finite stellar-fuel exhaustion, analytic capture fractions, agreement with the earlier solver when escape is zero, and equal eventual deposits for different emission durations with the same fuel. The largest tested total-energy residual was 1.78 × 10^-15 in normalized energy units.

## Conditional derivation from stated assumptions

Use total energies in a fixed control volume, not densities in a changing volume. Let P be photon energy, C traveling companion energy, D deposited energy, X energy transferred into a specified but as-yet-unmodeled receiving sector, O cumulative energy exported across the boundary, and F available radiative fuel. All rates have units of inverse time; j is energy per time.

Assume constant nonnegative rates: h conversion, g capture, ep photon escape, ec companion escape, lc companion energy transfer, and ld deposit energy transfer. Then:

```
dP/dt = j - (h + ep) P
dC/dt = h P - (g + ec + lc) C
dD/dt = g C - ld D
dX/dt = lc C + ld D
dO/dt = ep P + ec C
dF/dt = -j
```

Adding the equations gives d(P+C+D+X+O+F)/dt = 0. The source switches off when the assigned fuel runs out. O records exported energy; it is not energy stored within the galaxy and is not available for local gravity merely because it appears in this accounting. Returning or incoming radiation requires extra boundary terms. X likewise needs a physical identity before this can be a closed physical theory.

Your suggested branch is lc = ld = 0. Setting ep = ec = 0 too recovers the closed three-reservoir equations, with the fuel account added. Escape is an independent choice from secondary energy loss.

For one initial unit of photon energy, no further injection, permanent deposits and positive h and g, integrate the photon equation to obtain ∫P dt = 1/(h+ep). Integrating the companion equation then gives ∫C dt = h/[(h+ep)(g+ec+lc)]. Therefore:

```
eventual deposited fraction = h/(h+ep) × g/(g+ec+lc)
```

If conversion or capture is zero, the deposited fraction is zero. This derivation uses fixed rates and a single mixed volume. It does not derive those rates, orbital trapping, a spatial halo, or microscopic particle behavior.

In the no-loss branch the second factor is g/(g+ec). Equal capture and escape rates retain half the converted energy locally; capture nine times faster than escape retains 90%. These are conditional examples, not measured rates. With a closed volume and positive conversion and capture, all initial photon energy eventually deposits. Finite time can yield much less.

## The supply–response tradeoff

Let Mreq denote the baseline effective extra mass inferred using the archived dynamical assumptions. Let eta denote gravitational response relative to the baseline energy-to-mass relation, Edep/c². It is a placeholder for a possible changed field law, not an energy source and not yet a theory consistent with lensing or conservation.

Let L0 be the catalog's present luminosity, T an illustrative emission duration, B the ratio of integrated available photon energy to L0 T, and f the fraction converted and deposited by the observation time. Then the necessary budget condition is:

```
eta × f × B × L0 × T >= Mreq × c²
eta_required = S × (10 billion years / T) / (B f)
S = archived full-conversion energy shortfall for that galaxy
```

The accompanying JSON evaluates 45 stipulated scenarios for each of the 26 matched galaxies. B and f are independent scenario inputs here; a physical transport calculation must determine them together. External sources may contribute to B only after their energy production, propagation, boundary inflow and allocation among receiving galaxies are accounted for. Energy cannot be assigned in full to multiple galaxies.

For B=f=1, the median required eta is:

| Assumed constant-luminosity duration | Median required response multiplier |
|---|---:|
| 10 billion years | 5,427.66 |
| 100 billion years | 542.77 |
| 1 trillion years | 54.28 |
| 10 trillion years | 5.43 |
| 100 trillion years | 0.54 |

Values below one indicate an energy surplus relative to this necessary budget, not a successful halo model. These numbers assume all the assigned photon energy can be captured by the observation time. The earlier nearby-catalog transport calculation has a much larger shortfall because it does not assume complete conversion and capture.

Merely making the universe older is insufficient if its total usable stellar fuel remains fixed. In the solver, stretching emission of the same fuel changes timing but not the eventual permanent deposit for fixed rates. An older-universe solution therefore needs a quantitatively specified integrated energy supply: source populations, fuel history, incoming radiation, or an explicitly modeled additional reservoir. The recovered fitted refractive-index history also cannot simply be extrapolated to these ages without revisiting its finite past endpoint.

## What remains for Codex to process

1. **Source-history branch:** specify age, stellar population and replenishment history; integrate emitted energy and fuel usage; track external inflows without counting them twice. Outcome: a justified B for each galaxy, with uncertainties.
2. **Capture branch:** compare free passage, scattering with energy retained, and a new binding interaction. Derive momentum exchange, residence times, spatial distribution and receiving-sector behavior. Outcome: f from transport, rather than a fitted retention assumption.
3. **Gravity-response branch:** formulate candidate field equations for eta and include all traveling and deposited energy. Outcome: predictions for dynamics and lensing, with energy–momentum conservation and stability checks. A multiplier alone is insufficient.
4. **Permanent versus finite storage:** preserve both possibilities until a physical storage mechanism is specified. Outcome: accumulation histories, escape/decay destinations, and halo evolution.
5. **Observational integration:** carry surviving conditional models into the full roadmap's optical, local-gravity, galaxy, cluster and cosmological tests. Passing an energy budget alone is not a full observational fit.

These are branches to investigate, not grounds for closing the fictional concept. The pending choice about prioritizing source history, gravitational response, or both remains with the user.

## Files

- `transport_ledger.py`: conditional rate solver.
- `check_transport.py`: reproducible verification and scenario generator; uses the existing local project paths.
- `transport-checks-and-tradeoffs.json`: nine verification results and 45 scenario maps over 26 galaxies.

This is a derivation from explicit rate-model assumptions. A derivation of those assumptions from an action or microscopic interaction is still outstanding.
