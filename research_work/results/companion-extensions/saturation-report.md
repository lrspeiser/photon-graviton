# Local donor saturation: a nonuniform stopping profile

13 September 2026. Hypothetical endpoint closure, not a derived transport mechanism.

## Outcome

A single density threshold can replace the uniform contraction with an ordered, mass-conserving donor profile. All six coupled orbit calculations converge with positive energy release. However, every case has larger stellar-speed RMS error than its previous uniform-contraction counterpart. The simple saturation closure is therefore not adopted as an improvement to the reference model.

The reason this matters: reaching the right average amount of settling does not guarantee that gravity ends up at the right radii. This calculation lets the distribution change rather than assuming that a new stopping rule preserves the old successful shape.

## Rule and provenance

**Known spherical volume and mass conservation:** for initial shell edges r_i and r_(i+1), V_i=4*pi*(r_(i+1)^3-r_i^3)/3 and donor density is m_i/V_i.

**Proposed packing closure for this project, using familiar capacity mathematics:**

\[
V_{f,i}=\min(V_i,m_i/\rho_{\rm sat})\quad(m_i>0),\qquad
V_{f,i}=V_i\quad(m_i=0).
\]

Shells remain ordered. Their final inner volumes are the cumulative sum of preceding final volumes, retaining the original innermost boundary. A representative packet retains its original fractional position in shell volume. This is not the shortcut r_f/r_i=(rho_i/rho_sat)^(1/3), which generally fails to account for the changing positions of neighboring shells.

**Algebraic consequence of this closure:** occupied shell density becomes max(rho_i,rho_sat). Initially denser shells are not expanded; the threshold arrests further compression, rather than enforcing an absolute ceiling on every initial state. Empty-shell volume is preserved by assumption. Other components can occupy the same volume, so this is a donor-density condition, not a ceiling on total gravity or total mass density.

The exact-third inventory and previous phenomenological donor selection remain fixed. Calibrating the donor-mass-weighted mean r_f/r_i to the previous fitted value 0.803193303 in baseline I gives

\[
\rho_{\rm sat}=2.2177344\times10^7\ M_\odot\,{\rm kpc}^{-3}
=0.022177344\ M_\odot\,{\rm pc}^{-3}.
\]

This number is an inverse calibration, not a measured companion density or a first-principles constant. No new speed optimization is used, but the previous contraction target itself came from a fit to these observations. The threshold is frozen for baseline II and all receiver choices. Baselines I/II describe the same Galaxy with different ordinary-matter assumptions; this is a sensitivity transfer, not a new-galaxy prediction.

Baseline II gives mean contraction 0.811565. Individual donor radius ratios span approximately 0.576–1 in I and 0.571–1 in II. Thus the proposed stop is substantially nonuniform even though the mean is near the previous value.

## All results

The coupled solver recalculates self-gravity, adjusts unforced circular orbits, and allocates donor angular momentum to the declared receiving band. The energy is the same Newtonian circular-kinetic plus external and self-binding ledger as the preceding calculation. Positive release is an energy requirement for outgoing channels, not a luminosity prediction or proof that such channels exist.

| Ordinary baseline | Initial receiver band (kpc) | Inner RMS (km/s) | Outer RMS (km/s) | All RMS (km/s) | Previous uniform all RMS | Release (10^50 J) |
|---|---:|---:|---:|---:|---:|---:|
| I | 15–30 | 6.00 | 8.40 | 7.23 | 6.27 | 3.725 |
| I | 30–60 | 6.00 | 8.36 | 7.21 | 6.28 | 4.569 |
| I | 60–120 | 6.00 | 8.36 | 7.21 | 6.28 | 4.944 |
| II | 15–30 | 6.42 | 11.98 | 9.47 | 9.26 | 3.506 |
| II | 30–60 | 6.42 | 9.78 | 8.19 | 8.10 | 4.372 |
| II | 60–120 | 6.42 | 9.78 | 8.19 | 8.10 | 4.759 |

No-settling controls give all-bin RMS 6.34/10.91 for I/II. Saturation therefore also worsens the original baseline-I prediction while still improving baseline II. These reused circular-speed estimates are not individual-star speeds or untouched data. The raw output's `imposed_scores` field means the nonuniform donor profile with the extended population held fixed; the uniform comparison in this table comes from the preceding coupled-torque results.

## Checks and limitations

Doubling source bins from 8,192 to 16,384 and force nodes from 16,384 to 32,768, without recalibrating the threshold, changes predicted speeds by less than 0.084 km/s and released energy by less than 0.045%. Positive shell volumes, ordered radii, the occupied-shell density identity, and coupled angular-momentum residuals are checked. These numerical tolerances do not measure observational or ordinary-matter uncertainty.

The previous ideal thermal-condensate interpretation remains unsupported. This closure describes where donor packets would stop under an assumed packing rule; circular orbits provide the assumed support afterward. It supplies no local torque law, timescale, dynamical accessibility, stability proof, or pressure energy. Unlike the explicit restoring-energy branch, it does not assign a new stored-energy term merely for switching transport off. A real interaction could still have such a cost and would need to include it.

The finite model discards extremely small donor weights under the inherited packet cutoff and retains empty shells. Consequently this is not a continuum proof about disconnected donor regions. A physical version must specify whether companions can cross those regions and whether available states depend on donor density, all companion density, or the ordinary-matter environment.

## Decision and next discriminating work

Keep the exact-third reference and the uniform endpoint results; do not promote this constant donor-density rule. Its failure identifies the needed distinction: saturation must generate the radial shape as well as the mean contraction.

A useful next candidate is environment-dependent capacity, for example a fixed ratio of companion occupancy to available states supplied by ordinary matter. That would be a new hypothesis requiring one explicit shared formula and transfer to a different galaxy, not a separate threshold at each radius. A phase-space capacity is another distinct option because it involves allowed positions and velocities rather than spatial density alone. Both must retain the energy and angular-momentum accounting already required. Neither is presently established.

Reproduce with `python research_work/results/companion-extensions/saturation.py`. See [protocol](saturation-protocol.md), [source](saturation.py), and [complete results](saturation-results.json).
