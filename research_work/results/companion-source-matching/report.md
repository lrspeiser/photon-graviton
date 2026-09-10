# Connecting the produced companion to the fitted galactic field

**The light-conversion calculation and the galactic wave fit are not yet one physical model.** Using the frozen lens-fit mass, ordinary photon-scale energy in a single companion quantum implies an extremely relativistic traveling state. The fitted galactic source instead assumes a slowly varying, nonrelativistic wave. The transformation between them must be supplied; calling both a companion does not perform it.

This check applies known kinematics to an existing fitted parameter. It is not a new fit, an astrophysical production rate, an estimate of total stellar photon supply, or a rejection of all stored-spacetime hypotheses. The calibration is read directly from `wave-lens-training/parameter-grid/selected-calibration.json`, with its hash recorded in `results.json`. Historical exposure labels in that frozen file describe the selection date; subsequent lens validation and test evaluations are recorded elsewhere and are not reopened here.

## 1. What the galaxy calculation actually assumes

**Known nonrelativistic field equations, applied hypothetically:** the stationary source used in the lens fit obeys

```
-(hbar^2/2m) laplacian psi + m Phi psi = mu psi
laplacian Phi_c = 4 pi G m |psi|^2.
```

Its frozen training-selected mass is m = 1e-24 eV/c^2. This is a calibrated candidate parameter, not a measured particle mass. The approximation requires momenta well below mc and energies close to mc^2 per excitation. The preceding [wave-source report](../self-consistent-wave/report.md) specifies this regime. Localized scalar configurations and their possible formation are established research topics, not a unique invention here; see the [boson-star review](https://doi.org/10.1007/s41114-023-00043-4). Existence of a stationary profile does not establish how photons populate it.

## 2. The resonant pulse example cannot be rescaled freely to that mass

**Known special-relativistic kinematics:** two vacuum photons of equal energy E meeting at full opening angle theta can produce one on-shell scalar of rest energy mc^2 only when

```
(mc^2)^2 = 2 E^2 (1-cos theta) = 4 E^2 sin^2(theta/2)
theta = 2 asin(mc^2/(2E))
gamma_scalar = 2E/(mc^2).
```

Quasi-parallel photon collisions are a known way to investigate low-mass scalar resonances; see [Homma, Hasebe and Kume](https://arxiv.org/abs/1405.4133). That experiment is not evidence for our mass or production process. This calculation provides necessary vacuum kinematics, not conversion probability, angular acceptance, cross-section, or a sky-averaged rate.

| Energy of each photon | Required full opening angle | Produced scalar Lorentz factor |
|---|---:|---:|
| 1 GHz, 4.136 micro-eV | 2.418e-19 rad | 8.271e18 |
| 1 meV | 1e-21 rad | 2e21 |
| 1 eV | 1e-24 rad | 2e24 |
| 1 keV | 1e-27 rad | 2e27 |
| 1 MeV | 1e-30 rad | 2e30 |

For two optical-scale 1 eV photons, the produced scalar would carry 2 eV and travel with 1-v/c about 1.25e-49. This is a traveling relativistic state, not the fitted slow galactic envelope. Merely having a tiny rest mass does not mean it carries little energy or can settle easily.

Conversely, producing this scalar at rest from two exactly opposing photons requires each photon to have energy 5e-25 eV, frequency 1.209e-10 Hz and wavelength 2.480e18 m. Those are formal vacuum values, not evidence that such freely propagating radiation exists in an astrophysical plasma.

The preceding dimensionless resonant experiment used scalar mass 2 and photon carrier 1. Mapping its scalar mass to 1e-24 eV/c^2 therefore also maps its photon carrier to 5e-25 eV. Mapping the carrier to optical 1 eV instead makes its scalar mass 2 eV/c^2. The same run cannot represent both optical photons and the frozen ultra-light galactic source simultaneously. Its mechanistic finding about spectra remains valid in its stated dimensionless setting.

## 3. Why slowing the same companion is not enough

**Known energy identity under a conditional same-species description:** an excitation moving at speed v has

```
E_cold = gamma_v mc^2,   gamma_v = 1/sqrt(1-v^2/c^2).
```

At the illustrative galactic speed 200 km/s, its kinetic/rest-energy ratio is 2.2253e-7. Thus E_cold is approximately 1.00000022e-24 eV for our frozen mass.

If a single incoming companion carries 1 eV and becomes one such slow excitation, only about 1e-24 of its incoming energy remains in that excitation. The remaining energy must go somewhere else. This is not energy destruction and does not prove that energy must escape the galaxy: an absorber, other field excitations, recoil or radiation could receive it. Those alternatives require explicit dynamics and a gravitational response.

To retain the whole 1 eV in slow excitations of this same species requires approximately

```
N_cold = E_in / (gamma_v mc^2) ≈ 1e24.
```

This is **conditional energy accounting using known physics**, not a claim that a real scalar has an exactly conserved particle number. In a classical field it expresses the large change in low-energy mode occupation required. The value varies by less than one part per million across the illustrative 100–300 km/s range. No galaxy age, stellar luminosity integral or photon-generation efficiency enters this calculation.

A lone on-shell scalar cannot simply decay in vacuum into many on-shell scalars of the same mass: its invariant rest energy m is less than the total rest energy Nm for N>1. A receiving environment or additional incident states can change the available energy-momentum and enable other processes. Likewise, under ordinary local vacuum dispersion an isolated photon cannot emit a massive scalar while remaining a photon: the final photon-plus-massive-scalar invariant mass is positive whereas the initial photon's is zero. These restricted statements do not apply unchanged to media, backgrounds, modified dispersion or multi-particle interactions.

## 4. What a completed branch would need

| Candidate continuation | What must be calculated | Consequence for existing fits |
|---|---|---|
| Same field, redistribution into many slow modes | Coupled production, scattering/relaxation, energy and momentum transfer, and population of supported states | Existing wave profile may remain relevant, but its normalization and formation cannot be imposed independently |
| Traveling companion absorbed by a different storage system | Transition probabilities, reverse emission, recoil, lifetime and stress-energy of the receiver | The fitted single-field mass is not automatically the receiver's parameter; recompute its gravity and lensing |
| Modified dispersion or stored geometric deformation | A complete field law covering clocks, propagation, conservation and gravitational response | Re-derive the nonrelativistic limit rather than silently retaining the current Schrodinger-Poisson source |

These are open branches, not three demonstrated solutions. The user's option of deposited well deformation remains allowed as a hypothesis. This calculation does not force the interpretation of deposited energy as individual cold particles; it identifies what is required **if we retain the specific wave model already used in our fits**.

The next joint mechanism must produce a phase-space or field-state source, not just a scalar total-energy injection term. That source must lead to a supported galactic profile and also predict the surviving light spectrum and event timing. No additional fit flexibility per galaxy is authorized by this result.

## Verification

`run.py` evaluates five production energy scales and fifteen storage energy/speed cases. It checks the photon-pair invariant independently at 110-decimal precision; maximum relative residual is below 1.89e-52. Stable formulas avoid rounding the required angles or speed deficits to zero. Reciprocal energy-partition checks pass. The source calibration remains unchanged and hashed; no held-out data or previous scores are modified.

Reproduce with `python research_work/results/companion-source-matching/run.py`.
