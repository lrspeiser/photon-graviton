# The energy ledger does not determine the observed light

This preliminary T04/T05 calculation compares three stipulated transport kernels. It does not adopt a microscopic interaction. All three can give the same mean photon loss and companion energy gain; their spectra and photon counts differ. The no-secondary-loss companion option remains explicit throughout.

## A common energy equation, three different outcomes

Let s be path length, h a nonnegative constant loss coefficient per length, and U the photon energy in an initially identical ensemble of photons with energy E0. For all three examples:

```
U(s)=U(0) exp(-h s)
companion energy gained=U(0)-U(s)
```

Here loss is entirely a transfer into companions. Subsequent companion capture does not alter this source-side accounting. No capture efficiency or gravity response is inferred from it.

| Stipulated mechanism | Photon count | Spectrum of surviving photons | Total photon energy fraction |
|---|---|---|---|
| Whole-photon conversion with removal rate h | Falls as exp(-hs) | A monochromatic input remains at E0 | exp(-hs) |
| Fractional transfer events | Preserved for event fractions below one | A distribution of shifted energies | exp(-hs) |
| Deterministic continuous energy transfer | Preserved | A monochromatic input remains narrow at E0 exp(-hs) | exp(-hs) |

Thus matching the photon-energy ODE cannot tell us whether we have explained redshift. Whole-photon conversion dims a line without shifting it in this example. This also limits what can be inferred from the earlier integrated source-budget tests: those tests constrain energy availability, not the identity of the optical mechanism.

## Fractional events: an exact energy-conserving example

Suppose events occur independently with constant rate kappa per length. At each event a photon gives a fraction epsilon of its current energy to a companion, while continuing with fraction q=1-epsilon. Take 0<epsilon<1. If N events occur:

```
N is Poisson with mean mu=kappa s
E_N=E0 q^N
sum of all companion transfers=E0(1-q^N)
```

The companion sum telescopes for every individual photon history, not only on average. This specifies energy only; directions, recoil, inverse processes and a physical event rate remain to be derived.

For a photon number distribution f(E,s), the corresponding energy-space master equation is:

```
partial_s f(E,s) = kappa [f(E/q,s)/q - f(E,s)]
```

Integrating over E preserves photon number. Its first energy moment obeys dU/ds=-kappa epsilon U, giving h=kappa epsilon. The second moment follows from the Poisson generating function:

```
<E>/E0 = exp(-mu epsilon)
<E²>/E0² = exp[mu(q²-1)]
Var(E)/<E>² = exp(mu epsilon²)-1
```

Define S=E0/<E>, a centroid energy-shift factor in this stipulated frame. It is not automatically an observed cosmological redshift until emission/detection clocks are specified. Then:

```
mu epsilon = ln S
relative energy RMS = sqrt(S^epsilon-1)
```

For this Poisson kernel, larger fractional events generate a broader line at the same mean shift. In the limit epsilon→0 with kappa epsilon=h fixed, the fractional width tends to zero and the process approaches deterministic drift. This provides an open direction toward preserving narrow spectra; it is not a derivation of a physical interaction or an infinite-event-rate completion.

As epsilon approaches one, photons with any event approach zero energy. At epsilon=1 the physical interpretation becomes removal of those photons; the remaining photons have their original energy. The number-preserving master equation above is not used at q=0.

## Illustrative spectral requirements

If the permitted added fractional RMS width is delta, this particular event model needs:

```
epsilon <= ln(1+delta²)/ln S
mu >= (ln S)²/ln(1+delta²)
```

At S=2, interpreting delta as an illustrative velocity width divided by c:

| Illustrative added RMS width | Maximum fractional energy transfer per event | Minimum mean event count |
|---|---:|---:|
| 1 km/s | approximately 1.61 × 10^-11 | approximately 4.32 × 10^10 |
| 10 km/s | approximately 1.61 × 10^-9 | approximately 4.32 × 10^8 |
| 100 km/s | approximately 1.61 × 10^-7 | approximately 4.32 × 10^6 |

These are scenario requirements, not measured exclusions. Actual spectra require intrinsic widths, instrumental response, calibration and the full possibly non-Gaussian line profile. The result applies to independent equal-fraction Poisson events. Correlations, coherent transfer, energy-dependent rates and different event distributions can change it; their kernels must be derived rather than borrowing this formula.

## Continuous drift and the time-dilation requirement

The deterministic energy law dE/ds=-hE gives:

```
partial_s f = partial_E(h E f)
f(E,s)=exp(hs) f_initial[E exp(hs)]
```

It preserves photon number and produces the same first energy moment as the other examples. This is a useful desired transport behavior, not a microphysical action. A phase-coherent realization must also explain how phase, frequency, group velocity, emission and detection relate. We have not established that realization by writing a drift equation.

If all emissions travel the same fixed path R at the same fixed speed v and emission/detection clocks use the same fixed time convention, then:

```
t_arrive=t_emit+R/v
d t_arrive/d t_emit=1
```

Neither stochastic energy transfer nor deterministic drift alone changes this mapping. Delays or scattering could broaden a pulse, but broadening is not automatically the required uniform dilation of the entire light curve. A time-dependent propagation law, clock relation or another explicitly derived transfer mechanism can alter the map. It must be evaluated jointly with the spectrum and the existing atom/cavity requirements.

This is a conditional result for fixed travel conditions, not a prohibition on a fictional time-dilation mechanism. The pending user question asks whether to compare a smooth wave/clock mechanism and many tiny events, or prioritize one.

## How this advances the roadmap

T04 now has exact energy-space benchmark kernels that a proposed microscopic calculation can be checked against. T05 has an explicit timing requirement that cannot be supplied by its photon-energy moment alone.

Next, a candidate should produce a joint kernel for output energy, direction and arrival time, conditional on the incoming photon and environment. Its energy and momentum moments must match the full conservation law; its photon-number, spectral-width and timing predictions must be retained rather than fitted independently. Whole-photon mixing remains a conversion comparison branch, and fractional/coherent mechanisms remain open redshift candidates.

T04/T05 remain incomplete: no rate has been derived from an adopted action; angular momentum, recoil, polarization, inverse processes, source/detector clock response and a joint observed likelihood remain unresolved. No new astronomical fit was performed.

## Verification

`frequency-transfer-checks.json` contains nine Poisson cases checked against explicit probability sums for the first and second energy moments, three whole-photon-removal cases, and 100 individual randomized fractional-transfer histories with explicit companion accounting. All three check groups passed. It also contains nine illustrative width requirements.

`check_frequency_transfer.py` reproduces the checks. These are mathematical checks of stipulated forward energy kernels, not verification of a microscopic theory or observed redshift law.
