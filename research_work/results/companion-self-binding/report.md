# Self-binding can hold a massive-wave proxy; sustained redshift remains separate

The user's proposed companion attraction leads to three distinct questions: can a bound configuration exist, can incoming companions populate it, and does it continue stretching light? This pass derives and tests part of each question without adopting a companion mass or new force.

The result is conditional. A nonrelativistic massive scalar proxy admits a collective gravitationally bound equilibrium supported by wave pressure. That is not an ordinary massless-graviton bond, a capture mechanism or a verified long-lived halo. Its stationary density gives no continuing stretch under the current optical law. A repeating optical response alternates stretching and compression instead of giving a universal positive duration stretch.

## Bound-state calculation

The [derivation](derivation.md) specifies the Schrödinger-Poisson Hamiltonian, a Gaussian comparison and the full spherical boundary-value problem. Two domain/tolerance choices solve a unit-mass nodeless stationary state, with independent normalization, energy and virial checks. The refined result is

| Quantity, in declared hbar=m=G=M=1 units | Result |
|---|---:|
| Kinetic/gradient energy | 0.05425640257 |
| Gravitational energy | -0.10851280517 |
| Binding energy | -0.05425640259 |
| Chemical potential, excluding rest energy | -0.16276920779 |
| Half-mass radius | 3.92510134 |

The independently integrated Gaussian has energy -0.05305165, above the numerical state as expected for a trial profile. Changing the domain from 30 to 40, doubling initial nodes and tightening the BVP tolerance changes binding energy by 3.42e-8 relatively and half-mass radius by 1.16e-7. The refined virial residual is below 2e-10. No full linear stability spectrum or dynamical formation simulation was performed.

The calculation assumes massive companions and a conserved nonrelativistic occupation. It does not identify the scalar with ordinary spin-2 gravitons, derive how photons create this occupation, or establish the preservation of its symmetry under an eventual photon coupling. There is no independent dark-matter population, expansion or Big-Bang input.

## Binding does not supply extra energy

The model fixes a relation between constituent mass, half-mass radius and circular speed. The [protocol](protocol.json) uses three illustrative target pairs to display that relation; they are not a fitted astronomical sample.

| Half-mass radius | Circular speed there | Required constituent rest energy | Magnitude of binding/rest-energy fraction |
|---:|---:|---:|---:|
| 1 kpc | 30 km/s | 8.95e-23 eV | 4.27e-9 |
| 10 kpc | 200 km/s | 1.34e-24 eV | 1.90e-7 |
| 100 kpc | 1,000 km/s | 2.69e-26 eV | 4.74e-6 |

A single constituent mass predicts a single product of half-mass radius and speed for this isolated ground-state family. The different masses in the table expose that restriction; they are not parameters silently fitted separately to actual galaxies. Other states, baryonic potentials or interactions require a new comparison. A finite single cloud also has an eventually declining exterior circular speed, not an automatically flat outer halo curve.

At the 10 kpc example, total rest energy is about 3.32e58 J. The binding correction is negative and only about 1.90e-7 of that total. Bonding cannot be counted as a large positive mass-energy bonus. Any fictional photon-supplied cloud would still have to obtain its energy from accounted production and capture.

The same tiny constituent mass gives a formal two-body gravitational Bohr radius about 2.43e121 m and binding magnitude about 4.91e-233 eV. The macroscopic state relies on a collective occupation, not compact atomic-style pair bonds. These are nonrelativistic massive-proxy formulas; no massless-graviton pair calculation is claimed.

For a slow incoming constituent, settling into the selected equilibrium requires releasing energy approximately equal to the magnitude of the chemical potential. At the 10 kpc example, that is about 5.69e-7 of its rest energy, plus any incoming kinetic energy. Twelve finite-addition ledgers close to better than 1.81e-12 relatively. They specify the required energy destination, not an interaction rate or an identified emission channel. Coherence, charge/occupation balance, ejecta, recoil and population history remain to be derived.

## Why a repeating bound wave is not automatically a steady redshift source

For a transparent periodic optical profile with period P and matched stationary clock standards, the arrival map obeys F(t+P)=F(t)+P. Its local stretch J=F' therefore averages to one over uniformly sampled emission phases. It cannot exceed one at every phase. A stationary density is the special case with fixed delay and J=1.

The optical test uses a prescribed cosine-squared spatial profile, offset 0.1, period 3, and modulation amplitudes 0, 0.05 and 0.1. It does not claim that this profile follows from the computed bound state. Twenty-four and forty-eight phase samples give 216 phase probes and six period-repeat controls. The finer sampling gives:

| Modulation amplitude | Sampled local stretch range | Mean local stretch | Mean received/emitted photon energy |
|---:|---:|---:|---:|
| 0 | 1 to 1 | 1 | 1 |
| 0.05 | 0.92866 to 1.07681 | 1 within 2e-14 | 1.00274192 |
| 0.1 | 0.86231 to 1.15965 | 1 within 5e-14 | 1.01100127 |

Some emission phases redshift, while others blueshift. The energy average is over equal-energy photons emitted uniformly in time; it is not the inverse of the mean stretch. Since energy ratio is 1/J, convexity requires its mean to be at least one in this prescribed model. Maintaining such a nontrivial oscillation while it gives energy to photons would require an accounted driver or change the field. This is not energy creation and not a closed photon/binding evolution.

A phase-sensitive coupling can also change the occupation conservation assumed in the massive proxy. Consequently, the stationary bound-state result and the optical diagnostic must not be combined into a claimed complete action. A growing population, nonperiodic field, different clock response, absorption or other source selection could change the premises and needs its own derivation. The calculation does not replace the user's cumulative-time proposal with ordinary gravitational endpoint time dilation.

## Verification and retained failure

The first optical pilot failed the unchanged 1e-7 frequency/event agreement gate, with a maximum residual about 2.68e-7 at the compact profile boundaries. The [failure record](pilot-verification.json), [original protocol](pilot-protocol.json) and [source snapshot](pilot-check.py.txt) are preserved. The diagnostic temporarily bypassed that gate solely to collect all residuals; its status remains failed.

The successful calculation keeps the same physical profiles and acceptance threshold. It tightens the imported probe solver to relative tolerance 1e-12, absolute tolerance 1e-14 and maximum step 0.01. Final frequency/event disagreement is below 8.45e-10; period-repeat arrival errors are below 1e-7, and phase-grid energy averages change by less than 1.8e-11. Static frequency/energy preservation agrees to roundoff. Bound-state normalization, chemical-potential, virial, Gaussian quadrature and domain checks also pass. These are numerical and mathematical checks, not observational acceptance.

The [saved results](companion-self-binding-results.json) retain both bound profiles, physical scaling and capture ledgers, all phase probes, convergence values and source hashes. The standard suite includes this diagnostic.

## Next step

Specify and derive a photon/companion interaction compatible with binding, energy and any occupation/charge constraints. Determine whether capture-driven, nonperiodic evolution can maintain the required optical time change without assuming a free pump or an independently fitted redshift rate. Ordinary gravitational self-confinement of massless waves, a separate new force, and a distributed cosmic source history remain distinct alternatives.

Stable halo support, actual capture throughput, available photon supply, matter clocks and a shared rotation/lensing response are not solved by this pass. The complete 20-task, 32-area goal stays active.

```sh
python research_work/results/companion-self-binding/check.py
```
