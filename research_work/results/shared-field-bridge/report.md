# BRIDGE-1B: the field that reddens light can make the companions, but it must already carry their energy

One interaction, calculated in a homogeneous closed volume: radiation gives energy to the propagation field n of the matched-wave action, the field produces companion excitations through m_C²(n) = m0² + g²(n − n\*)², and those excitations react back on the field. There is one accounting of the photon's loss. This is the project owner's BRIDGE-1B proposal, implemented independently from their description; their prototype package was not readable from this session.

- **The source is no longer a set of knobs.** Production happens in one non-adiabatic crossing of the mass minimum, with the spectrum exp[−π(k² + m0²)/(g|ṅ\*|)]. Afterwards the momenta are frozen and the mass grows as g(n − n\*), so
  - the mass production rate is q = N g ṅ, tied to the same ṅ the redshift measures;
  - the companion density follows the redshift history, ρ_C(t) ∝ n(t) − n\*;
  - the birth speed is not chosen: v_rms(t) = k_rms/m_C(t).
- **Coldness is not the difficulty.** At the reference point the products' r.m.s. speed is 2.3×10⁻¹⁹ km/s one Gyr after the crossing. A companion travels 89 m before it is slower than 3 km/s, and 427 m in the whole 9.8 Gyr since. The prototype's 0.054c is a property of its code units (its run ends at m/k ≈ 20), not of the mechanism.
- **The energy is the difficulty, and it is the field's own.** Between the crossing and today the radiation delivers 1.5×10⁻⁶ of the companions' rest energy. The photons can pay in full only if the history begins at rest at a turnaround, where the radiation must hold (R_E + 1) times the companions' energy: at R_E = 10 that is n_s = 4.3×10⁻⁸, a 6.3×10⁷ K background 12.7 Gyr ago, whose cooled relic would be today's microwave background.
- **Producing companions bends the brightness the wrong way, and the field must be heavy to stop it.** The drain decelerates n, which makes distant sources brighter, and PF-1 is already too bright at high redshift. Keeping the Pantheon+ χ² within 1 of PF-1's coasting fit needs the field's present kinetic energy to be **at least 145 times the companions' rest energy**: 3.0×10⁸ times the microwave background's energy, about 1.4×10⁴ times the critical density.
- **A driving potential moves it the right way.** If a linear potential supplies 4–6 times the companions' drain, the history accelerates instead, and the fit improves to χ² = 840–841 against PF-1's 871.55 (flat FLRW: 836.51 at one parameter). That is a one-parameter fit to exposed data, and the potential's energy is stored energy, not light.
- **The optical gate fails, now with a number.** Within the minimal coupling family, fixed rulers and a fixed fine-structure constant force the hyperfine-to-optical ratio to drift at exactly 2p·ṅ/n, where p is the fraction of the coordinate shift that survives as measured redshift. The Yb⁺/Cs bound then gives **p ≤ 6.6×10⁻⁷**: a homogeneous propagation field cannot supply the measured redshift in this family, by six orders of magnitude. None of the five tested completions passes the four optical gates.
- **A prediction, if it were true.** The companions' density tracks the redshift history: ρ_C(z)/ρ_C,0 = 0.56 at z = 0.5, 0.34 at z = 1, and zero at the crossing redshift, 3.2 for the reference point.

Protocol: [protocol.md](protocol.md), declared before execution. Code: [modes.py](modes.py) (the semiclassical integrator), [bridge.py](bridge.py) (driver) and [checks.py](checks.py) (suite job). Results: [bridge-results.json](bridge-results.json).

**Why it was run.** The owner's proposal of 15 September: make the shared interaction the research target, derive the companion source from the dynamics that changes light, and stop treating the source rate and the birth speed as independent knobs. Their prototype showed the energy chain closing but its products hot, a source-power ratio of about 2.75 million, and a damping approximation that moves the brightness prediction the wrong way. This calculation takes the first two of their three calculations: the production spectrum with feedback, and the measured shift against atoms. The third, the RC-2a interface, is specified here and deferred.

## The model

In a homogeneous closed volume, with ħ = c = 1 in the frame the matched-wave action selects:

- **Radiation.** u_γ = A/n, so u̇_γ = −(ṅ/n)u_γ. This is the only photon loss in the calculation.
- **The field.** K n̈ = A/n² − V′(n) − g²(n − n\*)·S, with V the renormalized potential.
- **The companions.** A free scalar with m_C²(n) = m0² + g²(n − n\*)², carried mode by mode in the adiabatic basis, so the occupation is |β_k|² and the energy density is Σ w_k ω_k|β_k|². The zero-point sum and its force are subtracted; what renormalization absorbs is reported below, not hidden.
- **Exactly one ledger.** E = Kṅ²/2 + V + A/n + ρ_C is conserved by these equations; the numerical residual is the first validation.

## Part A: the code-unit trials

| Trial | What it tests | Result | Tolerance | |
|---|---|---|---|---|
| V1 | the crossing spectrum against exp[−π(k²+m0²)/(g\|ṅ\|)] | 3.1×10⁻⁴ absolute, 0.22% relative | 10⁻³, 2% | pass |
| V2 | energy conservation in every trial | ≤ 10⁻¹¹ relative | 10⁻⁸ | pass |
| V3 | the back-reaction's stopping point against K ṅ²/(2Ng) | 31.168 against 31.006, 0.52% | 5% | pass |
| V4 | frozen momenta and v_rms = k_rms/m_C after the crossing | occupations 7×10⁻⁴, speed law 0.13% | 10⁻³, 1% | pass |
| V5 | the owner's four controls | all four behave | — | pass |
| V6 | the arrival-interval stretch against the carrier redshift | 1.3×10⁻⁹ after first-order extrapolation, convergence ratio 2.000 | 10⁻⁷ | pass |
| V7 | freezing the occupations once every mode is adiabatic | stopping point moves 2×10⁻⁴ | 1% | pass |

**The controls.** With no radiation and the field at rest, nothing is produced. With the coupling removed, the field takes exactly what the radiation loses and no companions appear. With the index frozen, there is neither redshift nor production. **With the field moving and no radiation at all, companions are produced and paid for by the field's own kinetic energy** — the control that matters for the energy question below.

**The reference trial** (radiation the only initial energy, field at rest, code units): the radiation loses 13.322, the companions gain 12.208 and the field's kinetic energy gains 1.114, with a total-energy residual below 10⁻¹¹. The owner's own reference trial (17.344, 6.549, 10.795) used parameters this session could not read; the structure, not the numbers, is what these two share.

## Part B: the physical mapping

Fixed: PF-1's archived rate (ṅ/n = cα = 74.618 km/s/Mpc), the microwave background (4.17×10⁻¹⁴ J/m³), and 2B-F1's reference source as the abundance scale (q\* = 1,437.8 M☉ kpc⁻³ Gyr⁻¹ over 10 Gyr: 2.1×10⁶ times the background's energy, about 93 times the critical density). RC-2a stage 1 showed that rate was selected by a 300 kpc zone, so it is a scale, not a target.

Free: the index at the crossing n\*, and R_E, the field's present kinetic energy in units of the companions' rest energy. Everything else follows.

| n\* | crossing (Gyr ago) | R_E | field energy (ρ_crit) | k\* (eV) | m_C today (eV) | v_rms at +1 Gyr (km/s) | travel since crossing (m) | Δμ(z=1) | Pantheon χ² | turnaround T (K) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.9 | 1.3 | 10 | 930 | 2.8×10⁻⁸ | 4.7×10¹⁶ | 1.5×10⁻¹⁹ | 276 | −0.096 | 896.65 | 6.3×10⁷ |
| 0.5 | 6.4 | 10 | 930 | 2.0×10⁻⁸ | 1.2×10¹⁷ | 2.1×10⁻¹⁹ | 390 | −0.058 | 896.10 | 6.3×10⁷ |
| 0.237 | 9.8 | 10 | 930 | 1.9×10⁻⁸ | 1.6×10¹⁷ | 2.3×10⁻¹⁹ | 427 | −0.039 | 887.15 | 6.3×10⁷ |
| 0.237 | 9.9 | 30 | 2,791 | 1.9×10⁻⁸ | 1.6×10¹⁷ | 2.3×10⁻¹⁹ | 430 | −0.013 | 876.51 | 1.8×10⁸ |
| 0.237 | 10.0 | 100 | 9,305 | 1.9×10⁻⁸ | 1.6×10¹⁷ | 2.3×10⁻¹⁹ | 431 | −0.004 | 873.01 | 5.8×10⁸ |
| 0.237 | 10.0 | 300 | 27,914 | 1.9×10⁻⁸ | 1.6×10¹⁷ | 2.3×10⁻¹⁹ | 431 | −0.001 | 872.04 | 1.7×10⁹ |

The full six-by-six scan is in the results file. The pattern does not depend on n\*: the microphysics moves by tens of per cent across the whole range, and the brightness cost is set by R_E alone.

**What the companion would be.** At the reference point (n\* = 0.237, R_E = 10): produced with momenta of about 1.9×10⁻⁸ eV, a number density of 3.4×10⁻⁶ m⁻³, a mass today of 1.6×10¹⁷ eV, a de Broglie wavelength of 96 m, and a coupling g = 2.1×10¹⁷ eV. The field's stiffness K = 2R_E ρ_C,0/h² corresponds to a decay constant √K of 75 reduced Planck masses, and 285 at the brightness threshold below.

**The energy ledger.**
- Between the crossing and today the radiation delivers u_γ,0(1/n\* − 1) = 1.3×10⁻¹³ J/m³, which is 1.5×10⁻⁶ of the companions' 8.7×10⁻⁸ J/m³. In present power the same statement is the owner's: 1.01×10⁻³¹ W/m³ against 2.77×10⁻²⁵ W/m³, a ratio of 2.75 million. Both numbers are reproduced here.
- A V = 0 history that begins at rest turns around at n_s = u_γ,0/[(R_E+1)ρ_C,0], and there the radiation holds everything the field and the companions later have. At R_E = 10 that is 1+z_s = 2.3×10⁷ and T = 6.3×10⁷ K, 12.7 Gyr ago; at the brightness threshold R_E = 145 it is 3.1×10⁸ and 8.4×10⁸ K. **So "photon-powered" is not impossible — it is a hot beginning**, and the microwave background today is that bath cooled by the same factor. Whether such a history is admissible is exactly the question RC-1 left open, and it is not settled here.
- Without that beginning, the field's energy is stored energy that the model does not explain, and the photon's contribution is a part in a million. The owner's rule applies: calling it the same field that affects light does not explain its supply.

**The brightness cost.** The companions' drain gives the field a constant deceleration, so n changed faster in the past, distances at a given redshift shrink, and distant sources get brighter. PF-1's residuals already say its sources are too bright. The cost falls as 1/R_E:

| R_E | field energy | Δμ(z = 1) | Pantheon+ χ² (960 rows, one offset) |
|---:|---:|---:|---:|
| 1 | 93 ρ_crit | −0.322 | 1086.06 |
| 10 | 930 ρ_crit | −0.039 | 887.15 |
| 30 | 2,791 ρ_crit | −0.013 | 876.51 |
| 100 | 9,305 ρ_crit | −0.004 | 873.01 |
| **145** | **1.4×10⁴ ρ_crit** | −0.003 | **872.55 = PF-1 + 1** |
| ∞ (coasting) | — | 0 | 871.55 |

**The owner's damping approximation, scored.** Their K n̈ + KΓṅ = u_γ/n gives D = (c/h)·ln[1+(1+γ)z]/(1+γ). The magnitudes at z = 1 reproduce their table exactly (−0.059, −0.274, −0.505 mag at γ = 0.1, 0.5, 1.0); on the same 960 rows those histories score χ² = 896.67, 1045.06 and 1298.32, against PF-1's 871.55. The direction and the size of their diagnostic both stand.

**A driven variant.** A linear potential V′(n) = −λ makes the net force a free parameter. Scanning the net deceleration â = (Ng − λ)/(Kh²) — the drain minus the drive, so â > 0 means the companions win and â < 0 means the potential does:

| â | λ / Ng | Δμ(z = 1) | χ² |
|---:|---:|---:|---:|
| +0.20 | (drain exceeds drive) | −0.113 | 925.02 |
| 0 | 1.00 | 0 | 871.55 |
| −0.10 | 2.53 | +0.063 | 852.58 |
| −0.20 | 4.05 | +0.132 | 841.25 |
| −0.30 | 5.58 | +0.207 | 840.43 |
| −0.40 | 7.10 | +0.292 | 854.46 |

So a potential that funds the companions and drives the field a little faster than they drain it improves PF-1's brightness fit by about 31 in χ², to within 4 of flat FLRW's 836.51 — with two free parameters against FLRW's one. **This is a fit to exposed data, not a prediction**, and the potential's stored energy is still unexplained. It is, however, the first place in this programme where the redshift history's known discrepancy moves in the right direction for a stated physical reason.

## Calculation 2: the measured shift, against atoms

The coordinate shift is the easy part: the histories give 1 + z = n_o/n_e, and the arrival-interval stretch equals it (V6). What the observer measures is (n_o/n_e)·(ν_atom,o/ν_atom,e), and the atomic factor is the completion's business.

**The identity.** Over the minimal family ε ∝ n^a, μ ∝ n^c, m_e ∝ n^b, e ∝ n^d, with the nucleus and matter's own limiting speed fixed:

- the measured shift exponent is p = b + 4d − 2a + (a+c)/2;
- rulers scale as n^(a−b−2d) and the Coulomb fine-structure constant as n^(2d−a);
- the hyperfine-to-optical ratio scales as n^(c+8d+2b−3a).

Requiring fixed rulers and a fixed fine-structure constant forces d = a/2 and b = 0. The gross structure is then automatically fixed, and the hyperfine-to-optical ratio drifts at exactly **2p·ṅ/n**. With ṅ/n = 7.6×10⁻¹¹ per year and the Yb⁺/Cs bound of 10⁻¹⁶ per year, **p ≤ 6.6×10⁻⁷**. A homogeneous propagation field therefore cannot deliver the measured redshift while leaving rulers, gross structure and α alone — not marginally, but by six orders of magnitude. Every completion must break one of those three, and the calculation says which:

| Completion | p | rulers | worst clock-ratio drift (per year) | O1 | O2 | O3 | O4 |
|---|---:|---|---:|:-:|:-:|:-:|:-:|
| M1 matched action, fixed charges and masses | −1 | fixed | 1.5×10⁻¹⁰ (fine/gross) | no | yes | no | yes |
| M2 fixed electrostatics, Z = 1/n | +1 | fixed | 1.5×10⁻¹⁰ (hyperfine/optical) | yes | yes | no | yes |
| M3 matched action, electron mass ∝ n² | +1 | shrink as 1/n | 1.5×10⁻¹⁰ | yes | yes | no | no |
| M4 universal clock, dτ = dt/n | 0 | fixed | 0 | no | yes | yes | yes |
| M5 CC-1 co-scaling | +1 | counts grow | 0 | yes | yes | yes | no |

No completion passes all four gates, so the "physical optical predictions" milestone is not met, and the bridge inherits the problem rather than repairing it — as the owner expected. The brightness numbers above are quoted under M2's convention, the only tested completion with the full shift and a fixed gross structure.

## What the source would hand to RC-2a (specified, not run)

- **One production epoch**, not a rate: everything is born at t\*, with momenta drawn from exp(−πk²/k\*²) and a number density N.
- **A mass that keeps changing**: m_C(t) = g(n(t) − n\*), so d(mv)/dt = −m∇Φ, which is a drag −(ṁ/m)v of 0.10 per Gyr today at the reference point, and the gravitating mass is the current one. A population whose mass grows sinks: this is the opposite of a fresh cold birth at every epoch, and RC-2a's engine does not do it yet.
- **The donor's depletion is the field's own energy loss**, uniform in space, so RC-2a's "daughters minus depletion" accounting carries over unchanged.
- **The transfer** Q_C = N g ṅ c² = 2.8×10⁻²⁵ W/m³ at the reference point, tied to the redshift rate.

## Diagnostics that limit the scope

- **Renormalization.** The one-loop scale m_C⁴(n)/64π² that the renormalized V(n) absorbs is 2×10⁷³ times the field's kinetic energy at the reference point. The calculation declares V; it does not explain why that cancellation holds.
- **Spatial physics.** Companions attract each other through n gradients as well as gravity, with a ratio c⁴/(4πG K v_n²(1−n\*)²) = 6×10⁻⁴ at the reference point (v_n = c), falling as 1/R_E, and the index a galaxy's companion overdensity would dig is δn ≈ 2×10⁻¹⁰ at a 200 km/s potential depth. Both are small, which is what makes the homogeneous treatment defensible — but no gradient, lensing, fifth-force or equivalence-principle calculation was done.
- **Gravity of the homogeneous energies.** A field carrying 10³–10⁴ times the critical density, with the stiff equation of state of a rolling scalar, has to sit in a static background. That is the unsolved problem the shared-interaction study already named, and nothing here solves it.

## What this shows, and what it does not

**Shown.**
- One interaction can carry photon energy into the field and the field's energy into companions, with feedback and an exactly closed ledger, and the production spectrum, the abundance in time, the speeds and the travel distance all follow from it rather than being chosen.
- The products are cold, easily: mass growth with conserved momentum is a cooling mechanism that needs no expansion.
- The photons' contribution is a part in a million unless the history begins hot, in which case the required beginning is explicit: 10⁷–10⁸ in redshift, 10⁷–10⁹ K.
- The brightness cost of draining the field, and the stored energy that keeps it acceptable: 145 times the companions' energy, 1.4×10⁴ times the critical density.
- A quantitative form of the atomic-reference problem: p ≤ 6.6×10⁻⁷ within the minimal family.

**Not shown.**
- That this companion is the companion of the earlier reports, or that the abundance is right: the abundance is a comparison scale from 2B-F1, and RC-2a stage 1 has already shown that scale was chosen by a boundary.
- Anything about structure: the calculation is homogeneous, and a cold uniform medium at this density is Jeans-unstable, as 2B-F1 reported.
- Any improvement to the measured redshift, which still fails by six orders of magnitude in this family.
- That the driven variant is anything but a two-parameter fit to exposed supernova data.
- The RC-2a interface, which is specified above and deferred by the owner's own order.

## Reproduce

```sh
python research_work/results/shared-field-bridge/bridge.py    # about a minute on one core
python research_work/results/shared-field-bridge/checks.py    # suite job: V1, V2, V3, the controls, the identity and the mapping
```

The driver compares its output with the archived results and overwrites them only with `--canonical`. Wall-clock timings are excluded from the comparison.
