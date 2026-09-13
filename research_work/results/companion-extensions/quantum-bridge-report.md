# Connecting redshift transfers to the bound-companion scale

13 September 2026. Cross-check of existing branches, not a particle detection or observational exclusion.

## Main finding

The settling branch's effective mass 17.7827941 eV/c^2 cannot simply be assigned to a particle created at every small redshift event. Under the existing illustrative Poisson linewidth model, the allowed individual transfers are many orders of magnitude smaller. Accumulation in an existing reservoir, soft collective excitations, or different traveling and bound descriptions remain possible, but need an explicit common interaction.

This distinguishes three quantities that previous branches do not yet connect: traveling quantum energy, bound-state creation gap, and the effective mass in the phenomenological phase selector. None has been measured as a companion property. The failed ideal thermal-phase support test means the fitted mass must not be treated as an independently established microscopic input.

## Existing work reused

The [capture ledger](settling-report.md) already includes reservoir recoil. [Bound-pair production](../bound-pair-production/report.md) and [production/loss balance](../pair-production-balance/report.md) already show reversible creation and inverse processes in conditional scalar models. This calculation does not repeat those results or claim permanent storage follows from them. It connects their unknown gaps to the [small-transfer spectral diagnostic](report.md).

## Formula provenance and assumptions

**Known Poisson probability identity applied to our hypothetical fixed-fraction loss process:** if N has Poisson mean Lambda and each event leaves fraction 1-delta of photon energy,

\[
E_N=E_0(1-\delta)^N,\qquad
\langle E_N\rangle=E_0e^{-\Lambda\delta},\qquad
\frac{\operatorname{Var}(E_N)}{\langle E_N\rangle^2}=e^{\Lambda\delta^2}-1.
\]

Set A=alpha*D=Lambda*delta using the archived fitted alpha, not a newly derived rate. For added relative RMS energy width w,

\[
\delta\leq\frac{\ln(1+w^2)}A,\qquad
\Lambda\geq\frac{A^2}{\ln(1+w^2)},\qquad
\Delta E_{\rm first}\leq E_0\frac{\ln(1+w^2)}A.
\]

These bounds apply to independent fixed-fraction Poisson events. They are not general limits on coherent frequency change, correlated events, deterministic transfer, arbitrary jump distributions or measured astronomical line widths. Later events from the same photon carry less energy because its remaining energy is lower.

**Known energy accounting:** if a separate creation process requires gap epsilon and draws only on accumulated transfers,

\[
N_{\rm contributions}\geq\left\lceil\epsilon/\Delta E_{\rm first,max}\right\rceil.
\]

This optimistic count assumes every contribution is maximal and can be delivered to the same store without loss. It is a count of energy contributions, not necessarily distinct photons or simultaneous events. Momentum, recoil, spatial collection, intermediate storage and reaction thresholds can make it insufficient. The whole-path quantity epsilon/[E_0(1-exp(-A))] is merely a photon-energy equivalent, not a deterministic minimum number of photons in a stochastic process.

## Computed scale comparison

Use D=30.660139 Mpc (approximately 100 million light-years), alpha=0.0002488993286/Mpc and illustrative w=10^-5. The mean-energy loss is 0.760224%, delta_max=1.31039e-8, and the minimum mean count is about 582,366 events per path. The table provisionally sets epsilon=17.7827941 eV solely to test identification with the phase parameter.

| Initial photon | Maximum first transfer (eV) | Minimum maximal contributions per hypothetical 17.8 eV gap |
|---|---:|---:|
| 1 GHz radio | 5.42e-14 | 3.28e14 |
| 100 GHz microwave | 5.42e-12 | 3.28e12 |
| 0.1 eV infrared | 1.31e-9 | 1.36e10 |
| 2 eV optical | 2.62e-8 | 6.79e8 |
| 10 eV ultraviolet | 1.31e-7 | 1.36e8 |
| 1 keV X-ray | 1.31e-5 | 1.36e6 |
| 1 MeV gamma | 0.0131 | 1,358 |

For optical photons the whole-path mean loss is about 0.0152 eV, or about 1,170 such path-energy equivalents per hypothetical single-state gap. The output also includes a two-constituent gap 35.5656 eV, and width allowances 10^-3 and 10^-6, without selecting one by observed agreement. None of these widths is an imported measurement.

## Can binding remove the mismatch?

A transition between levels of a preexisting store can have a very small gap even when its constituents have larger masses. That requires an actual level spectrum; it does not create the whole store from one small transfer. The seed store and its energy must remain in the ledger.

For hypothetical production from scratch with total constituent rest energy epsilon_rest, a small net gap could result from negative binding energy: epsilon_gap=epsilon_rest-B. If the gap is as small as the optical event budget, binding would cancel nearly all of the assumed constituent rest energy. One cannot then count the full 17.8 eV as newly deposited gravitational energy while accounting for only the tiny input: binding and interaction stress belong in the gravitational source. A nonrelativistic ideal-gas mass parameter cannot automatically be carried into that extremely bound regime.

## Traveling energy and bound effective mass need separate definitions

**Known special-relativistic dispersion**, not unique to this project:

\[
E^2=p^2c^2+m^2c^4,\qquad v/c=\sqrt{1-(mc^2/E)^2}.
\]

A massive free particle does not move exactly at c at finite energy. For a literal 17.8 eV rest mass, E=10 mc^2 gives a speed deficit about 0.5%; E=10^7 mc^2 gives approximately 5e-15. These are kinematic examples, not gravitational-wave observational bounds on companions. See the [Particle Data Group kinematics review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf).

Massless traveling waves and massive or effectively massive bound excitations are a possible two-stage description, but the transition must be derived. Calling the fitted mass an effective property of a medium avoids identifying it with a free particle's dispersion; it also requires deriving that medium's modes before reusing phase or support equations.

## What an achromatic receiving spectrum would need

**Known mean-loss bookkeeping:** a fixed gap epsilon with event rate Gamma(E) per path length gives alpha(E)=Gamma(E)*epsilon/E. Constant alpha would require Gamma(E) proportional to E, wherever the event is allowed. This adjusts the mean only; it does not repair thresholds, fluctuations, image changes or supernova durations. The earlier point-response pair interaction gave a different color dependence and remains unsuitable as the direct redshift law in that limit.

**Candidate interaction target, not a derived law:** a scale-invariant fractional jump kernel can be written dGamma/dDelta=(lambda/E)*f(Delta/E), where f is a normalized distribution on fractional loss 0<x<1. Its mean fractional-loss rate is lambda*integral x*f(x) dx, independent of E. An interaction must generate that spectrum and sufficiently small higher moments across the tested bands; writing it down does not demonstrate it. A fixed discrete gap does not have this scaling across arbitrary photon energies.

## Decision

Keep traveling-companion energy, bound creation gaps and the fitted phase mass explicitly separate. Do not use one-to-one production of 17.8 eV particles as the microscopic redshift mechanism. Preserve accumulated-store and collective-field branches, but require their mode spectrum, transfer/reverse rates and energy accounting before identifying them with the rotation prescription.

The next foundational calculation should test a proposed continuous receiving spectrum against both mean redshift and spectral broadening, retaining the earlier clock/event-duration requirement. This complements the environmental-capacity work; neither track alone completes the theory.

Reproduce with `python research_work/results/companion-extensions/quantum-bridge.py`. [Source](quantum-bridge.py) and [complete numerical record](quantum-bridge-results.json).
