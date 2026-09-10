# Companion speed must track the optical environment

**Setting the companion's fixed coordinate speed to one does not make it travel at light's local speed throughout this model.** The existing photon ray speed is c0/n_bar, while the field's principal wave speed is a separately fixed v. Their ratio is v*n_bar/c0. A shared change of clock or ruler units multiplies both speeds equally and cannot remove their mismatch.

This is a conditional consequence of the existing equations, not a claim about measured gravitons. It clarifies an unfulfilled user requirement: freely traveling companions were intended to move at the speed of light. Earlier half-speed scalar calculations are useful diagnostics but cannot be called that completed branch.

## Fixed-speed comparison

The same resolved smooth pulse is evolved at companion speeds 0.5, 0.75 and 1 in reference units with c0=1. Signal energies are 0.01 and 0.1. No redshift coefficient, smoothing, inertia, initial field or driving packet is retuned. Each case uses 32 and 64 field modes. The speed ratio is sampled where signal packets lie between the source and detector planes, not across every point and time of the box.

| Signal energy | Fixed companion speed | Sampled companion/photon speed ratio | Mean reference carrier stretch | Energy-width stretch |
|---:|---:|---:|---:|---:|
| 0.01 | 0.50 | 0.5019 - 0.5701 | 1.137366 | 1.137595 |
| 0.01 | 0.75 | 0.7529 - 0.8609 | 1.144222 | 1.144832 |
| 0.01 | 1.00 | 1.0038 - 1.1508 | 1.139778 | 1.140487 |
| 0.1 | 0.50 | 0.5020 - 0.5887 | 1.259238 | 1.258991 |
| 0.1 | 0.75 | 0.7530 - 0.9117 | 1.317028 | 1.319321 |
| 0.1 | 1.00 | 1.0041 - 1.2399 | 1.322069 | 1.328906 |

The fixed-speed-one field is faster than the photon rays by up to approximately 15.1 percent in the weaker case and 24.0 percent in the stronger case. These are dimensionless toy-model values, not astronomical measurements. A valid shared clock convention cannot make unequal speeds equal at the same location. The comparison also shows that changing the field speed alters redshift and timing, so the original results cannot simply be relabeled as light-speed companions.

The half-speed reference reproduces the archived smooth-pulse results. Energy and momentum checks pass at every speed, and mode refinement changes the reported carrier stretch, energy-width stretch and field-energy gain by less than 1e-6. Conservation alone therefore does not establish the intended propagation speed.

## A candidate change derived from a positive Hamiltonian

Provenance: a proposed modification constructed with established canonical variational calculus. It is not certified original, a derivation of the existence of time, or a validated theory of gravitons.

For a positive dimensionless local field n, consider replacing the fixed-speed gradient energy with

`H_field = integral [Pi^2/(2K) + K*c0^2/2 * |gradient ln(n)|^2] dx`.

Both contributions are nonnegative when K>0. Keeping n as the canonical coordinate gives, in one dimension and without photons,

`n_t = Pi/K`

`Pi_t = K*c0^2/n * partial_xx ln(n)`

or equivalently

`n_tt = c0^2/n^2 * n_xx - c0^2/n^3 * n_x^2`.

The second term is required by the Hamiltonian. Merely replacing a constant wave speed in the old differential equation and omitting that term would change the energy accounting. A coupled photon term contributes its own functional derivative as before; that full finite-pulse system has not yet been executed for this new field.

Linearizing about constant n=N>0 gives `delta_n_tt = (c0/N)^2 delta_n_xx`. Thus small field waves and photon rays share coordinate speed c0/N in a uniform background. If the same local clock/ruler convention applies to both, their local speeds are equal as well. This does not derive that convention or its material realization.

## Executed checks of the proposed field term

An independent finite-difference derivative of the discretized logarithmic-gradient energy agrees with its canonical force to approximately 1.8e-10 absolute in the chosen units. Small-amplitude waves were then evolved at background values N=1, 1.2 and 2, on 128 and 256 cells. Their measured phase speeds approach 1, 0.833333 and 0.5, respectively, with the expected grid dispersion. Measured frequencies agree with the discrete linear prediction to better than 2.4e-11 relative, and relative energy drift remains below 3.4e-11.

These tests establish positive-energy variational consistency and the linear speed property. They do not establish global nonlinear stability, finite photon conversion, transient stretching, astrophysical normalization, local clock compatibility, or gravity. The optical model averages n over a finite kernel, whereas this field's local characteristic speed depends on n itself. Exact speed matching in a strongly varying environment therefore also requires a consistent local or smoothed coupling; uniform-background agreement cannot hide that distinction.

## Equal speed and lossless travel are separate requirements

This change does not establish the stipulated no-loss companion transport. In a slowly evolving homogeneous background, the proposed small-wave dispersion is `omega_companion = c0*k/N(t)`, the same frequency dependence as the photons. At fixed wave number, increasing N lowers the wave frequency. In the usual adiabatic wave-action limit, its reference energy also scales with that frequency. This is an established WKB consequence applied conditionally, not a new observational result or a completed energy-exchange calculation.

Thus matching the speeds can reintroduce the very companion redshift that the intended lossless branch forbids. A complete model must specify how companions maintain their energy and where any compensating work comes from, while preserving the required photon shift and local detector standards. Calling the waves gravitons does not supply those equations.

The next full candidate must jointly satisfy speed, energy retention and observable redshift/timing. The old fixed-speed model and the newly proposed logarithmic-gradient model should remain distinct until that calculation passes. Neither is adopted here as the completed theory. No astronomical data or holdouts were opened, and the total photon-supply budget remains deferred.

## Reproduction

Run `run.py`, `matched_speed.py`, `report.py` and `verify.py`. `kernel.py` is a parameterized copy of the previous resolved-pulse solver, with the fixed-speed comparison checked against its archived half-speed output. The small-wave calculation is separate and saves its results in `matched-speed-check.json`. The former solver retains a fixed v; it must not be described as having implemented the new variable-speed field.

- [Resolved smooth pulse](../continuous-signal-feedback/report.md): baseline and packet-refinement checks.
- [Finite signal feedback](../finite-signal-feedback/report.md): energy exchange and endpoint clock controls.
- [Matter-clock closure](../matter-clock-closure/report.md): shared clock transformations and measured redshift limitations.
- [Source-timescale calculation](../source-timescale/report.md): the earlier subluminal approximation whose singular point-source formula must not be extrapolated to v=1.
