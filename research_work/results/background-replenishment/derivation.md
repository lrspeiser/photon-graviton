# Can continuing emission preserve the microwave background?

## Why we need this calculation

The previous pass found that lowering each photon's energy while keeping every photon changes the normalization of an initially thermal background. A continuing supply of light is a possible alternative history. Before invoking it, we must check whether ordinary positive emission can maintain the desired spectrum. Otherwise a fitted source could quietly require negative emission: removing photons while being described as a source.

This is a conditional inverse calculation. It finds what a source and sink would have to do; it does not derive actual stars, dust, absorption, thermalization or background formation. No expanding universe, Big Bang, dark matter population or thermal companion bath is assumed. A time-dependent or nonthermal history remains a separate problem.

## 1. Conservation in frequency space

Let N_nu be photon number per physical volume per frequency, j_nu the nonnegative emission rate in the same units per time, and kappa_nu a photon-removal rate. Assume a homogeneous fixed physical volume, fixed light speed and frequency drift

    d nu/dt = -H_c nu,       H_c > 0.

H_c is a conversion rate, not an expansion rate or Planck's constant. The frequency-space continuity equation is

    partial_t N_nu = H_c partial_nu(nu N_nu) + j_nu - kappa_nu N_nu.

The last term removes entire photons. It is distinct from partial energy conversion. Absorption transfers the removed energy to matter or another field; escape transports it elsewhere. Neither deletes energy globally. Re-emission, if present, belongs in j_nu.

For a stationary Planck target at fixed T, write

    x = h_P nu/(k_B T),
    N_nu = (8 pi/c^3) nu^2/(exp(x)-1),
    nu partial_nu ln N_nu = 2 - x/(1-exp(-x)).

The usual fixed-speed Planck occupation and mode count are reviewed in [Fitzpatrick's statistical-mechanics notes](https://farside.ph.utexas.edu/teaching/sm1/Thermalhtml/node103.html). We use this spectrum as a mathematical target, without importing a cosmological origin.

The required emission is therefore

    j_nu/N_nu = kappa_nu + H_c [x/(1-exp(-x)) - 3].

This equation is a requirement, not a new microscopic interaction law.

## 2. Emission alone cannot maintain this stationary target

With no removal, j_nu is negative below x = 2.821439372... and positive above it. Conversion moves photons from higher to lower frequencies; the low-frequency part needs a drain to avoid accumulating excess photons. Adding positive emission everywhere cannot supply that drain.

The same conclusion follows from total photon count: conversion preserves count, so a stationary finite photon count with no exits cannot accept a positive net photon source. This statement assumes frequency-boundary fluxes vanish, as they do for the Planck target. A condensate or unresolved zero-frequency reservoir would change that boundary assumption and needs its own accounting.

For a frequency-independent removal rate kappa, nonnegative emission at every positive frequency is possible if and only if

    kappa >= 2 H_c.

This follows because x/(1-exp(-x)) increases from 1 to infinity. More generally, frequency-dependent removal only needs

    kappa_nu >= H_c max(0, 3-x/(1-exp(-x))).

The gray bound is not a universal absorption bound. It does not apply unchanged to frequency-dependent drift, changing electromagnetic modes, a finite measured band, or a time-dependent target. Nor does mathematical nonnegativity establish that actual emitters can generate the required spectrum.

## 3. The source and the discarded photons have energy costs

Integrating the continuity equation with vanishing endpoint fluxes, for gray kappa,

    d N_total/dt = J_number - kappa N_total,
    d U_gamma/dt = P_source - (H_c+kappa) U_gamma.

A stationary bath thus requires

    J_number = kappa N_total,
    P_source = (H_c+kappa) U_gamma,
    P_conversion = H_c U_gamma,
    P_removed = kappa U_gamma.

At the smallest allowed gray rate, kappa=2H_c, the bath needs three units of input energy for every unit transferred by gradual conversion. The other two units leave through whole-photon removal. If that removal also feeds companions or deposits, it is an additional interaction channel and its energy and momentum must be counted explicitly.

The mean injected photon energy must be (1+H_c/kappa) times the mean bath photon energy. At kappa=2H_c this is 1.5 times the bath mean, or approximately 4.052 k_B T. This is an exact moment requirement for the prescribed source, not evidence that a stellar population has the required spectrum or power.

Maintaining the photon bath does not maintain the whole universe in a steady state. Fuel must decrease or some independently accounted energy source must provide P_source; receiving fields, matter or exported radiation take the outputs. Permanent deposits with nonzero capture continue to accumulate. Recycling the outputs into photons requires a physical process and cannot also retain the same energy permanently in deposits.

## 4. A same-temperature thermalizer is not enough by itself

A passive absorber/emitter at the target temperature with local detailed balance has j_nu=kappa_nu N_nu at that temperature. These two terms cancel on the target spectrum. The conversion drift then remains uncompensated. Strong thermalization can reduce deviations, but exact stationarity under nonzero drift requires an additional net source/heat flow or a different state. A realistic temperature distribution, opacity and heating budget must be solved together; simply naming dust or a thermalizer does not provide that solution.

## 5. Numerical checks and limits

The companion script independently checks the required source against a finite-difference transport derivative, integrates both source moments, and integrates the characteristic solution with the nonnegative prescribed sources at kappa/H_c=2,3,10. The unchanged target is recovered. This tests the inverse construction; no astronomical data are fitted and no physical source feedback or formation stability is demonstrated.

The next physical background model must specify its emission spectrum and energy supply, photon removal or escape, conversion rate and frequency dependence, material heating/re-emission, and companion/deposit destinations. These must predict the background spectrum and angular structure jointly with the other observations. The required source above cannot be treated as an observed or derived source merely because it reproduces the chosen target.
