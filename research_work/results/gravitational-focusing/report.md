# CF-1 report: a slow companion bath is retained by clusters, overheats disks, and erodes everything at 3,000 km/s

**Result.** This is a conditional transport-and-capture diagnostic, run on potentials built from ordinary matter plus a counted seed.
- **Transport is verified.** Focusing raises the entry rate and the unbound density exactly as the analytic factors predict. With no capture, nothing is retained.
- **The declared interaction:** elastic, equal-mass scattering between incoming and bound companions. With it, net retention happens only when the bath is slower than the host's escape speeds.
  - At 300 km/s, every system retains. Coma retains 140–290 times faster than the Milky Way per unit seed and incident density, and its captured population stays bound.
  - The Milky Way's captured population overheats: each retained companion brings 26 times the seed's binding energy.
  - At 1,000 km/s, only the labeled full-well diagnostic retains.
  - At the trial speed of 3,000 km/s, every system erodes, including the full well.
- **Supply.** Retaining a mass equal to the baryons in 10 Gyr needs, at σ/m = 1 cm²/g, incident densities 10⁴–10⁶ times the cosmic mean matter density.

So one slow population under one law does produce the depth-keyed pattern: clusters retain, massive ellipticals retain marginally, and disks overheat. But only for a slow bath, and only with a supply that nothing yet explains.

Protocol: [protocol.md](protocol.md), declared in 012f6ca before execution. Code:
- [focus.py](focus.py): transport;
- [scatter.py](scatter.py): the interaction;
- [cf1.py](cf1.py): driver.

Results: [cf1-results.json](cf1-results.json). Suite job: [checks.py](checks.py).

## Validation (no capture), declared tolerances

| Check | Worst case over all systems and trial speeds | Tolerance |
|---|---|---|
| V1: entry area from integrated orbits against the exact effective-potential area, at the half-mass radius | 9.4×10⁻⁶ | 10⁻³ |
| V2: time-weighted density from orbits against √(1 + v_esc²/u²), at 0.1–10 half-mass radii | 2.1×10⁻⁸ | 10⁻³ |
| V3: retained fraction with capture disabled | 0 | 0 |
| Interaction: stationary-target capture and ejection against v_esc²/w² and u²/w² | 0.8% (Monte Carlo) | 1% |

- **The analytic entry formula** πR²(1 + v_esc²/u²) equals the exact area at every half-mass radius tested; no centrifugal barrier limits entry there.
- **Energy and momentum** are conserved to rounding in every scattering event.

## Potentials (baryons plus a 1% seed; escape speeds relative to R_b)

| System | Mass within R_b (Msun) | Half-mass radius | Escape speed: centre / half-mass radius (km/s) |
|---|---|---|---|
| Milky Way, baryon model I | 1.0×10¹¹ | 6.4 kpc | 689 / 310 |
| Six SLACS lens hosts, stars with Auger Chabrier masses | 1.6–4.1×10¹¹ | 6–15 kpc | 920–1,079 / 341–479 |
| Coma, low end of the gas and star brackets | 1.4×10¹⁴ | 1.54 Mpc | 1,099 / 798 |
| Coma, high end | 2.7×10¹⁴ | 1.54 Mpc | 1,506 / 1,093 |
| Inverse diagnostic, NFW 10¹⁵ Msun (c = 4), truncated at R₂₀₀ (labeled) | 1.0×10¹⁵ | 0.97 Mpc | 4,025 / 2,681 |

Escape speeds relative to infinity are at most 2% higher.

**A correction to the focusing estimates in the owner's discussion.** The Coma boost factors quoted there, about 2–3 at 3,000 km/s, used Coma's full inferred well (about 4,500 km/s). From baryons alone, Coma's central escape speed is 1,100–1,500 km/s.

## The declared interaction: growth and heating kernels

K is the net retention rate per unit seed mass and per unit (σ/m)ρ_∞, in km/s. Positive means the seed grows; negative means the bath erodes it. The heating column gives the energy each net retained companion brings into the bound population, over the seed's specific binding energy; above 1, capture unbinds more than it adds.

| System | K at u = 300 km/s | Heating per retained / binding at 300 | K at 1,000 | K at 3,000 |
|---|---|---|---|---|
| Milky Way | +10.1 | 26 | −908 | −2,969 |
| Lens hosts (six) | +124 to +424 | 0.49–1.94 | −783 to −873 | −2,928 to −2,957 |
| Coma, low | +1,432 | 0.07 | −477 | −2,823 |
| Coma, high | +2,953 | −0.02 (cools) | −22 | −2,668 |
| Full well (inverse, labeled) | +19,690 | −0.10 | +4,996 | −1,002 |

**Retention needs a bath slower than the host.**
- **Why.** For stationary targets the net retained per event is (v_esc² − u²)/(v_esc² + u²), and the targets' own motion makes ejection easier. So a bath faster than the host's escape speeds knocks out more companions than it captures.
- **The 3,000 km/s trial.** It is faster than every baryonic well here, and even the full well erodes at the seed's mass-weighted escape speed.

**Depth-keyed retention at 300 km/s.**
- **Contrast.** Coma's K is 142–293 times the Milky Way's, and 3.4–24 times the lens hosts'.
- **Focusing supplies only part of it.** At the half-mass radius the entry factor is 2.1 for the Milky Way and 14 for Coma's high end.
- **Retention supplies the rest.** In deeper wells a larger share of each event ends bound, as the owner's review anticipated: the fraction of kinetic energy that must be removed falls with depth.

**Energy decides whether a reservoir can form.**
- **The ledger.** Elastic scattering conserves the companions' energy, and each incomer arrives with +u²/2.
- **Milky Way.** It retains net companions at 300 km/s but gains 26 times their binding energy, so without an energy sink its captured population cannot stay bound.
- **Coma.** It keeps its captured population bound, and at the high end the population even cools, because ejected companions carry away more than incomers bring.
- **Lens hosts.** They are marginal.
- **What is missing.** CF-1 adds no sink. "What receives the energy" is still the open question for galaxies.

**Feedback.** Retention deepens the well, which raises retention further. Every system with K > 0 runs away once its exposure passes a threshold, so the outcome is a switch rather than a gradual trend.

## Supply (reported separately)

This is the exposure (σ/m)ρ_∞ for the retained mass to reach the baryonic mass in 10 Gyr, which serves only as a reference scale. It is converted to an incident density at σ/m = 1 cm²/g, a trial scale, in units of the cosmic mean matter density (H0 = 70, Ω_m = 0.3). The bath speed is 300 km/s.

| System | Exposure (kpc⁻¹) | Incident density / cosmic mean at 1 cm²/g |
|---|---|---|
| Milky Way | 1.5×10⁻² | 1.8×10⁶ |
| Lens hosts | 0.8–2.5×10⁻³ | 1.0–2.9×10⁵ |
| Coma, low / high | 2.6 / 1.3×10⁻⁴ | 3.0 / 1.5×10⁴ |
| Full well (inverse) | 1.9×10⁻⁵ | 2.3×10³ |

Even the most favorable baryonic case needs a bath 15,000 times denser than the cosmic mean at 1 cm²/g, or a cross-section larger in proportion. No such population or supply is identified. That is CC-2's question.

## Wave check

The de Broglie wavelength at 300 km/s is 30 kpc for the archived 1.34×10⁻²⁴ eV constituent and 0.4 kpc for 10⁻²² eV; at 3,000 km/s it is 3.0 and 0.04 kpc.
- The particle treatment holds for the heavier constituent everywhere, and for Coma at any tested mass.
- For the light constituent at retaining speeds, the wavelength exceeds the galaxies' half-mass radii (6–15 kpc). A wave treatment would be needed there, and CF-1's galaxy numbers do not apply to it.

## What this establishes and what it does not

**Established, conditionally:**
- Transport and focusing, verified.
- Under one declared interaction and one shared slow bath, retention is strongly depth-keyed. Clusters retain and stay bound, massive ellipticals are marginal, and a Milky Way–like disk overheats.
- A fast bath (1,000–3,000 km/s) erodes baryonic wells rather than filling them.

**Not established:**
- That such a bath exists or is supplied. The required density is extreme (CC-2).
- That clusters' observed excess comes from this. No lensing or mass profile was predicted, and Coma's baryonic inputs are brackets awaiting a matched model.
- Anything for galaxies, whose captured population cannot stay bound here without an energy sink. Their extra gravity would need another source or another interaction.
- Any fit. Nothing was adjusted to data.

**What would test it next:**
- a declared energy sink, or a second interaction, for galaxies;
- the equilibrium profile a retained cluster population would reach, to be compared with lensing under identical conventions;
- the supply mechanism (CC-2).

## Reproduce

```sh
python research_work/results/gravitational-focusing/cf1.py      # about 2.5 minutes; fails fast and logs progress
python research_work/results/gravitational-focusing/checks.py   # suite job, 13 s
```
