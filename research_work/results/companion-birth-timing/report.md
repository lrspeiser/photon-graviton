# Packet birth settings determine whether the internal reservoir stretches events

**The energy-retaining reservoir does not by itself predict event stretching.** If every packet starts its internal variable at zero, successive packets take identical travel times and a supernova's duration is unchanged. If packet birth settings track a shared evolving reference, the same stipulated propagation rule stretches the event by the spectral factor. The previously exposed supernova aging measurements favor the latter timing behavior under the published source-template assumptions.

This is a conditional mechanism calculation plus a refit of already examined data. It is not fresh holdout evidence, a derived photon-conversion interaction, or a successful full theory. It connects the [internal reservoir construction](../internal-companion-reservoir/report.md) to the actual temporal observable, instead of treating energy conservation or two-frequency optical beats as an event-duration test.

## Formula provenance and assumptions

The ideal reservoir gave `ds/dt=1` and packet group speed `dx/dt=1/(1+a*s)` in dimensionless reference units with c0=1. Its propagation energy was `p/(1+a*s)`, with the complementary change stored internally. These remain hypothetical packet dynamics. Applying that energy to photon frequency is an additional optical interpretation, not a derived electromagnetic field interaction.

**New project parameterization of an initial-state choice, without a fundamental novelty claim:**

`s(t_birth)=kappa*t_birth`, hence `s(t)=kappa*t_birth+(t-t_birth)`.

Here kappa measures how birth settings change between emission events. It is not a new value of the ordinary clock rate. Kappa=0 resets the packet variable; kappa=1 synchronizes it with the selected common reference clock. This family assumes a fixed source-detector separation D and uses fixed ordinary endpoint clocks for its main comparison. The zero point and history of the background index are part of the stipulated environment.

**Conditional derivation using elementary integration:** writing flight duration as T and `n_birth=1+a*kappa*t_birth`,

`D = integral_0^T du/(n_birth+a*u) = ln[(n_birth+a*T)/n_birth]/a`.

Define `S=exp(a*D)`. The spectrum's reference energy ratio is 1/S for every packet. The arrival time is

`t_arrival = [1+kappa*(S-1)]*t_birth + (S-1)/a`.

Therefore, for two events at the same source and path,

`S_event = Delta_t_arrival/Delta_t_birth = 1+kappa*(S-1)`.

This is a conditional result of the stated birth rule and group velocity, not a newly invented exponential redshift law. It uses the same known accumulating-index kinematics as earlier project calculations; the new test identifies which initial states the internal-reservoir proposal must supply.

| Birth rule | Spectral stretch | Whole-event stretch | Interpretation |
|---|---:|---:|---|
| Reset, kappa=0 | S | 1 | Identical flight delay for every emitted packet |
| Partial tracking, kappa=0.5 | S | (1+S)/2 | Only part of the required arrival-interval effect |
| Shared phase, kappa=1 | S | S | Later packets experience longer flight times in this reference |

For an illustrative S=2, two source events separated by ten reference days arrive ten, fifteen or twenty days apart respectively. This is an example of the formula, not an observed pair of known emission/arrival dates.

The 54 numerical trajectories integrate the same canonical internal-reservoir equations, stop at the detector, and verify the analytic arrival times, energy transfer and 36 pairwise interval factors. Initial propagation energy is one, initial reservoir energy is separately one, and the propagation momentum is chosen accordingly at birth. This normalization does not derive emission or supply that seed. No source/detector motion, spatial gradients, waveform phase dynamics or ordinary-matter clocks are included in those trajectories.

## Comparison with the already exposed spectral-aging measurements

We reuse the unchanged 35-row Table 3 transcription from the [electromagnetic audit](../electromagnetic-audit/report.md). These measurements compare supernova spectral evolution with elapsed observer time using local spectral templates. They do not directly measure absolute flight times. The primary source is [Blondin et al. (2008)](https://arxiv.org/abs/0804.3595). Its expansion interpretation is not used as a premise of this fit.

Conditioning on measured z, and assuming the template supplies the appropriate intrinsic aging rate, the prediction becomes

`r_age = 1/[1+kappa*z]`.

**Status: conditional observational mapping, fitted by known weighted least squares.** There is no distance in this fit. It does not calibrate a, establish z(D), or show that the transfer mechanism caused the measured z.

| Model | kappa | Diagonal residual chi-square |
|---|---:|---:|
| Fixed birth reset | 0 | 150.5693 |
| Fixed shared phase | 1 | 26.9485 |
| Exploratory fit | 0.964237 | 26.8668 |

The formal delta-chi-square-one interval for fitted kappa is 0.844995 to 1.092880, within a search interval 0 to 4. The fixed cases have 35 residuals; the fitted model has one additional fitted parameter. The tiny improvement over kappa=1 is not a reason to replace exact synchronization with the fitted number. The two fixed scores reproduce old project results under a new mechanism interpretation; they are not independent confirmations to add to earlier evidence.

These scores assume diagonal published summary errors. Shared templates, population differences and selection are not fully represented; the interval is not a robust physical confidence bound on a universal clock. Allowing one global intrinsic-aging normalization as an explicitly exploratory sensitivity gives normalization 0.968768, kappa 0.845284, and chi-square 24.6122. This sensitivity is not a demonstrated correction or adopted source model.

More generally, the summaries measure `r_intrinsic(z)/[1+kappa*z]`. An arbitrary redshift-dependent intrinsic aging law can mimic a different kappa. Such a source law would require independent evidence, rather than per-object adjustments to rescue a preferred transport model. No covariance inflation, object removal or new holdout opening was used here.

## The clock cancellation remains

The main table treats ordinary endpoint clocks as fixed. For the kappa=1 branch, instead give them the common rate `q(t)=1/(1+a*t)` from the earlier matter-clock proposal. Then

`Delta_tau = integral q(t) dt = Delta_ln(1+a*t)/a`.

Since `1+a*t_arrival=S*(1+a*t_birth)`, finite source and arrival proper-time intervals are equal. Likewise local photon energy `E_reference/q` has the same value at emission and reception. Both the measured duration stretch and spectral shift disappear. Two numerical finite-interval checks reproduce this cancellation. It is not repaired by the presence of a conserving internal reservoir.

For birth-reset or partially synchronized packets, the internal index can differ for different packets at one spacetime event. Such a model cannot silently identify every packet's index with one universal local clock or metric; it needs a species/state-dependent interaction and a physical optical channel. In particular, matching a group-velocity and energy ansatz does not establish a consistent electromagnetic phase field.

## Next requirement and consequence

The shared initial phase is now a concrete requirement of this candidate, rather than an unspecified hope that a long journey stretches events. The next joint mechanism must prepare successive packet states from an evolving physical reference while explaining ordinary clock behavior and counting the energy of that reference. Simply resetting independent companion clocks at every conversion will not supply the needed duration effect in this model.

The [coherent receiver audit](../coherent-receiver-audit/report.md) already showed why a stationary receiver cannot provide arbitrary coherent fractional frequency remapping under its assumptions. This birth-state result is complementary: it concerns event arrivals, not that receiver's two-frequency visibility. Neither substitutes for a coupled, energy- and momentum-conserving photon-creation/transport/capture action.

The observational preference here supplies a target for the creation mechanism. It does not establish that the target exists physically, passes local standards, predicts the baryonic/companion field, or fits lensing. The full stellar orbital comparison, source-to-gravity coupling and other electromagnetic checks remain required. The total universe photon-supply calculation remains deferred.

## Reproduction

Run `python research_work/results/companion-birth-timing/run.py`. `results.json` records numerical trajectory checks, exact source and code hashes, fitted parameters, the source-normalization sensitivity, and the common-clock cancellation controls. `predictions.csv` retains every object, observed value, uncertainty and each prediction. No parameter from the older frozen distance, gravity or holdout pipelines is modified.
