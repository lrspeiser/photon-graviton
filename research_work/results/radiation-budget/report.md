# RC-1 report: radiation cannot supply stage 2A's companions except through a hot, photon-removing history

The check covers stage 2A's 33 combinations within its regime. For each, it compares the companions' energy per unit volume with the radiation present today, and scores three declared ways of taking energy from the microwave background against FIRAS.
- **Energy.** At the required incident densities, the companions carry 1.1×10⁵–2.5×10⁸ times the energy of the microwave background (rest energy included, since their gravitating mass is what the benchmarks count). Against the microwave plus UV-to-mm background, the factor is 9.7×10⁴–2.1×10⁸. Present radiation holds the requirement in none of the combinations.
- **Laws (a) and (b) against FIRAS.** The diagnostic χ² may rise by at most 4 above no transfer.
  - Law (a), a frequency shift at fixed photon number, can then take at most 3.3×10⁻⁵ of the microwave background's energy.
  - Law (b), removing whole photons, can take 8.4×10⁻⁵.

  That is 10⁻¹³–7×10⁻¹⁰ of the requirement. FIRAS permits neither law in any combination.
- **Law (c), the Planck-preserving combination.** It is the only radiative route left. With the microwave background as the donor, it needs:
  - a past temperature of 50–340 K (69–420 K for the inner benchmark B2);
  - the removal of all but 5×10⁻⁷–1.6×10⁻⁴ of that background's photons;
  - 38–63 Gyr at the archived loss rate.
- **Its cost to every other beam.** Applied at the same rate to all photons, law (c) dims sources by 0.081 magnitudes at 100 Mpc, 0.81 at 1 Gpc and 2.4 at 3 Gpc. Keeping that under 1% at 1 Gpc needs a coupling at least 74 times weaker for starlight than for microwaves.
- **Controls.**
  - A thermal input stays Planck only under law (c). After a 1% transfer it deviates by 10⁻⁹, against 3.3% under law (a) and 1.1% under law (b).
  - A diluted 5,800 K blackbody stays equally far from any thermal spectrum under all three, so none of the laws thermalizes.

Protocol: [protocol.md](protocol.md), declared in cc02334 before execution. Code: [rc1.py](rc1.py), which is also a suite job. Results: [rc1-results.json](rc1-results.json).

**Why it was run.** It follows the owner's direction after 2c9a110. Before stage 2B selects an interaction, every energy transfer into companions must predict what it does to the photons: their energy, their number and their spectrum.

**What was already established.**
- **The fixed-volume spectrum** ([thermal conversion](../thermal-conversion/derivation.md)). A frequency shift at fixed photon number turns a Planck bath into q⁻³·B_ν(qT). Keeping a Planck shape needs photon removal at three times the shift rate.
- **Earlier supply comparisons.** Under other interpretations, radiation fell 10³–10⁶ short.

This check adds the supply the complete formation mechanism needs ([stage 2A](../companion-formation/report.md)), and asks which transfer law and history could deliver it.

## Validation (declared tolerances)

| Check | Tolerance | Result |
|---|---|---|
| W1, Planck number and energy densities | 10⁻¹⁰ | 2.2×10⁻¹⁶ |
| W2, law (a) keeps photon number and scales energy by q | 10⁻¹⁰ | 2.2×10⁻¹⁶ |
| W3, the thermal-conversion pass's FIRAS scores at q = 1, 0.9999, 0.999 | 0.01 | 45.02, 86.00, 4,339.94 reproduced |
| W4, law (c) keeps a Planck spectrum over ten e-folds, and splits energy 1:3 | 10⁻⁹ | 3.8×10⁻¹⁵; split 3.000 |

## The requirement against radiation

**Radiation today.**
- The microwave background at 2.72548 K holds 4.17×10⁻¹⁴ J/m³.
- The UV-to-mm background holds 1.9–7.1×10⁻¹⁵ J/m³ over its historical range of 45–170 nW m⁻² sr⁻¹.

**The requirement.** For each in-regime combination at B1, the table gives:
- the companions' energy over that radiation;
- the fraction of the requirement that laws (a) and (b) could supply within the FIRAS comparison;
- what law (c) would need.

| Combination | ρ_B1 (M_sun/kpc³) | energy over the microwave background | over microwave plus UV–mm (170) | fraction law (a) could supply | law (b) | law (c): past temperature, K | photons removed | duration at α₀, Gyr |
|---|---|---|---|---|---|---|---|---|
| Milky Way S 1, gravity on | 1.6×10⁵ | 2.4×10⁷ | 2×10⁷ | 1.4×10⁻¹² | 3.5×10⁻¹² | 191 | 1 − 2.9×10⁻⁶ | 55.7 |
| Milky Way S 0.1, gravity omitted | 1.7×10⁶ | 2.5×10⁸ | 2.1×10⁸ | 1.3×10⁻¹³ | 3.4×10⁻¹³ | 342 | 1 − 5×10⁻⁷ | 63.3 |
| Milky Way S 1, gravity omitted | 6×10⁵ | 8.7×10⁷ | 7.4×10⁷ | 3.8×10⁻¹³ | 9.7×10⁻¹³ | 263 | 1 − 1.1×10⁻⁶ | 59.9 |
| Milky Way S 10, gravity omitted | 2.4×10⁵ | 3.6×10⁷ | 3×10⁷ | 9.3×10⁻¹³ | 2.4×10⁻¹² | 211 | 1 − 2.2×10⁻⁶ | 57 |
| Milky Way M 0.1, gravity on | 2×10⁵ | 2.9×10⁷ | 2.5×10⁷ | 1.1×10⁻¹² | 2.9×10⁻¹² | 201 | 1 − 2.5×10⁻⁶ | 56.3 |
| Milky Way M 1, gravity on | 1.1×10⁵ | 1.6×10⁷ | 1.4×10⁷ | 2.1×10⁻¹² | 5.3×10⁻¹² | 172 | 1 − 4×10⁻⁶ | 54.3 |
| Milky Way M 0.1, gravity omitted | 1.2×10⁶ | 1.8×10⁸ | 1.5×10⁸ | 1.8×10⁻¹³ | 4.7×10⁻¹³ | 316 | 1 − 6.4×10⁻⁷ | 62.3 |
| Milky Way M 1, gravity omitted | 4.3×10⁵ | 6.2×10⁷ | 5.3×10⁷ | 5.3×10⁻¹³ | 1.4×10⁻¹² | 242 | 1 − 1.4×10⁻⁶ | 58.8 |
| Milky Way M 10, gravity on | 7.5×10⁴ | 1.1×10⁷ | 9.3×10⁶ | 3×10⁻¹² | 7.8×10⁻¹² | 156 | 1 − 5.3×10⁻⁶ | 53.1 |
| Milky Way M 10, gravity omitted | 1.6×10⁵ | 2.3×10⁷ | 2×10⁷ | 1.4×10⁻¹² | 3.7×10⁻¹² | 189 | 1 − 3×10⁻⁶ | 55.5 |
| J1630 S 0.1, gravity omitted | 1.2×10⁶ | 1.7×10⁸ | 1.5×10⁸ | 1.9×10⁻¹³ | 4.9×10⁻¹³ | 313 | 1 − 6.6×10⁻⁷ | 62.1 |
| J1630 S 1, gravity on | 1.5×10⁵ | 2.2×10⁷ | 1.9×10⁷ | 1.5×10⁻¹² | 3.8×10⁻¹² | 186 | 1 − 3.1×10⁻⁶ | 55.4 |
| J1630 S 1, gravity omitted | 3.9×10⁵ | 5.6×10⁷ | 4.8×10⁷ | 5.9×10⁻¹³ | 1.5×10⁻¹² | 236 | 1 − 1.5×10⁻⁶ | 58.5 |
| J1630 S 10, gravity omitted | 1.5×10⁵ | 2.1×10⁷ | 1.8×10⁷ | 1.6×10⁻¹² | 4×10⁻¹² | 185 | 1 − 3.2×10⁻⁶ | 55.3 |
| J1630 S 10, gravity on | 10⁵ | 1.5×10⁷ | 1.3×10⁷ | 2.2×10⁻¹² | 5.7×10⁻¹² | 169 | 1 − 4.2×10⁻⁶ | 54.1 |
| J1630 M 0.1, gravity on | 3×10⁵ | 4.4×10⁷ | 3.8×10⁷ | 7.5×10⁻¹³ | 1.9×10⁻¹² | 222 | 1 − 1.8×10⁻⁶ | 57.7 |
| J1630 M 0.1, gravity omitted | 9.1×10⁵ | 1.3×10⁸ | 1.1×10⁸ | 2.5×10⁻¹³ | 6.4×10⁻¹³ | 293 | 1 − 8.1×10⁻⁷ | 61.3 |
| J1630 M 1, gravity on | 1.4×10⁵ | 2×10⁷ | 1.7×10⁷ | 1.7×10⁻¹² | 4.2×10⁻¹² | 182 | 1 − 3.4×10⁻⁶ | 55.1 |
| J1630 M 1, gravity omitted | 2.9×10⁵ | 4.3×10⁷ | 3.7×10⁷ | 7.7×10⁻¹³ | 2×10⁻¹² | 221 | 1 − 1.9×10⁻⁶ | 57.6 |
| J1630 M 10, gravity on | 6.8×10⁴ | 10⁷ | 8.5×10⁶ | 3.3×10⁻¹² | 8.5×10⁻¹² | 153 | 1 − 5.6×10⁻⁶ | 52.8 |
| J1630 M 10, gravity omitted | 10⁵ | 1.5×10⁷ | 1.3×10⁷ | 2.2×10⁻¹² | 5.6×10⁻¹² | 170 | 1 − 4.1×10⁻⁶ | 54.2 |
| Coma, low bracket S 0.1, gravity omitted | 8,835 | 1.3×10⁶ | 1.1×10⁶ | 2.6×10⁻¹¹ | 6.6×10⁻¹¹ | 91.8 | 1 − 2.6×10⁻⁵ | 46.1 |
| Coma, low bracket S 1, gravity omitted | 2,821 | 4.1×10⁵ | 3.5×10⁵ | 8×10⁻¹¹ | 2.1×10⁻¹⁰ | 69 | 1 − 6.2×10⁻⁵ | 42.3 |
| Coma, low bracket S 10, gravity omitted | 939 | 1.4×10⁵ | 1.2×10⁵ | 2.4×10⁻¹⁰ | 6.2×10⁻¹⁰ | 52.4 | 1 − 1.4×10⁻⁴ | 38.7 |
| Coma, low bracket M 0.1, gravity omitted | 7,640 | 1.1×10⁶ | 9.5×10⁵ | 3×10⁻¹¹ | 7.6×10⁻¹¹ | 88.5 | 1 − 2.9×10⁻⁵ | 45.6 |
| Coma, low bracket M 1, gravity omitted | 2,440 | 3.6×10⁵ | 3×10⁵ | 9.3×10⁻¹¹ | 2.4×10⁻¹⁰ | 66.6 | 1 − 6.9×10⁻⁵ | 41.9 |
| Coma, low bracket M 10, gravity omitted | 795 | 1.2×10⁵ | 9.9×10⁴ | 2.8×10⁻¹⁰ | 7.3×10⁻¹⁰ | 50.3 | 1 − 1.6×10⁻⁴ | 38.2 |
| Coma, high bracket S 0.1, gravity omitted | 7,519 | 1.1×10⁶ | 9.4×10⁵ | 3×10⁻¹¹ | 7.7×10⁻¹¹ | 88.2 | 1 − 3×10⁻⁵ | 45.6 |
| Coma, high bracket S 1, gravity omitted | 2,415 | 3.5×10⁵ | 3×10⁵ | 9.4×10⁻¹¹ | 2.4×10⁻¹⁰ | 66.4 | 1 − 6.9×10⁻⁵ | 41.8 |
| Coma, high bracket S 10, gravity omitted | 795 | 1.2×10⁵ | 9.9×10⁴ | 2.9×10⁻¹⁰ | 7.3×10⁻¹⁰ | 50.3 | 1 − 1.6×10⁻⁴ | 38.2 |
| Coma, high bracket M 0.1, gravity omitted | 7,441 | 1.1×10⁶ | 9.3×10⁵ | 3×10⁻¹¹ | 7.8×10⁻¹¹ | 87.9 | 1 − 3×10⁻⁵ | 45.5 |
| Coma, high bracket M 1, gravity omitted | 2,373 | 3.5×10⁵ | 3×10⁵ | 9.5×10⁻¹¹ | 2.4×10⁻¹⁰ | 66.1 | 1 − 7×10⁻⁵ | 41.8 |
| Coma, high bracket M 10, gravity omitted | 779 | 1.1×10⁵ | 9.7×10⁴ | 2.9×10⁻¹⁰ | 7.4×10⁻¹⁰ | 50 | 1 − 1.6×10⁻⁴ | 38.1 |

The requirement is an environmental density around each system. If it held everywhere, the companions' mean density would also exceed the cosmic mean matter density by the same factors as in stage 2A (19–42,000).

## Laws (a) and (b) against FIRAS

Both laws leave a spectrum that FIRAS can tell from a blackbody.
- **Law (a)** keeps photon number, so the spectrum is q⁻³ times a Planck spectrum at the lower colour temperature.
- **Law (b)** removes photons at every frequency, so the spectrum is a scaled Planck spectrum.

With the colour temperature and the Galaxy coefficient refitted, the diagnostic χ² (45.02 with no transfer) rises by 4 at:
- an energy fraction of 3.3×10⁻⁵ for law (a), colour temperature 2.72493 K;
- 8.4×10⁻⁵ for law (b), colour temperature 2.72506 K.

These are diagonal-error diagnostics without the FIRAS calibration model, not official limits. But the requirement exceeds them by 9 to 13 orders of magnitude, far beyond any plausible change in the likelihood.

## Law (c): the history it would need, and its cost

Law (c) keeps a Planck spectrum while its temperature falls, so FIRAS cannot see it today. What the check can report is what it would take:
- **The past.** The microwave background would have had to hold E_req more energy: a past temperature of 50–340 K for B1, and 69–420 K for B2.
- **The photons.** All but 5×10⁻⁷–1.6×10⁻⁴ of its photons would have been removed.
- **The duration.** 38–63 Gyr at the archived rate, whose temperature e-folding time is 13.1 Gyr. No age is imposed; the duration is reported, not judged.
- **The cost to other beams.** If the same removal acts on every beam, sources dim by the extra photon-number factor e^(−3α₀D): 0.081 magnitudes at 100 Mpc, 0.81 at 1 Gpc and 2.4 at 3 Gpc. The existing brightness fits do not contain that. Keeping it below 1% at 1 Gpc needs the optical removal rate to be at most 0.013 of the microwave rate.

## Controls

| Input | Before | Law (a), 1% | Law (b), 1% | Law (c), 1% |
|---|---|---|---|---|
| Thermal, 2.72548 K: largest deviation from the best-fitting unit-amplitude Planck spectrum | 5×10⁻⁹ | 0.033 | 0.011 | 10⁻⁹ |
| Starlight proxy, a 5,800 K blackbody diluted 10⁻¹³: the same measure | 1.2×10¹³ | 1.2×10¹³ | 1.3×10¹³ | 1.2×10¹³ |

Only law (c) keeps a thermal spectrum thermal, and none of the three makes a nonthermal spectrum thermal. The plan's thermalization test therefore still needs an explicit absorber or emitter, which none of these laws is.

## What this does and does not show

**Shown.**
- **Present radiation falls short.** The radiation present today cannot supply stage 2A's requirement: it is 10⁵–10⁸ times short.
- **FIRAS rules out laws (a) and (b)** at that scale.
- **What law (c) would need:** a much hotter past, the removal of nearly all the original photons, tens of billions of years, and a coupling that spares starlight.

**Not shown.**
- **A thermalization mechanism**, or the origin of the microwave background.
- **Official FIRAS limits.**
- **A source.** Nor its kinematics: one photon cannot become one slow massive companion while conserving momentum.
- **Any change to the redshift or brightness fits.** The dimming is reported as a cost.

## Consequences for stage 2B

- **A source drawing on present radiation cannot meet stage 2A's requirement.** Stage 2B must do one of two things.
  - **A non-radiative source** with its own energy ledger, such as the owner's example of a field decaying into companion pairs.
  - **Law (c) with its history.** It must supply:
    - a microwave background that was once 50–340 K;
    - a photon-removal channel into companions that spares starlight by a factor of at least 74;
    - an origin for the hot initial state, which the recording-transition branch would then have to explain.
- **The requirement stands either way.** It is stage 2A's: 19–42,000 times the cosmic mean density in companions around the systems tested.

## Reproduce

```sh
python research_work/results/radiation-budget/rc1.py   # a few seconds; suite job
```

It compares its output with the archived results exactly and overwrites them only with `--canonical`.
