# Conversion-first exploratory research pass

The user permits conversion without special void/time stretching as its cause. The original time idea was intended to cause light stretching, including for nearby galaxies; it did not remove local redshift. This pass does not establish a microscopic mechanism or complete the research goal.

## Galaxy redshift result

Fitted alpha = 0.0002488993286 per Mpc = 7.631288108e-05 per million light-years. The equivalent c*alpha is 74.6181 km/s/Mpc (a conversion coefficient, not an assumed expansion rate).
Training sky-tile bootstrap 95% interval for c*alpha: 72.417 to 76.887 km/s/Mpc.
All 164 fixed published distances are unchanged, spanning 10.20 to 93.20 Mpc. The primary error model uses a fixed illustrative 300 km/s residual scale. Distance uncertainties are retained as provenance, not fitted.

| Partition | Count | Exponential RMSE km/s | Linear control RMSE km/s | Archived-rate RMSE km/s | Exponential nominal 95% coverage |
|---|---:|---:|---:|---:|---:|
| train (previously exposed) | 104 | 457.577 | 457.240 | 459.419 | 0.865 |
| validation (previously exposed) | 35 | 437.065 | 437.213 | 450.240 | 0.857 |
| test (previously exposed) | 25 | 415.414 | 414.659 | 429.124 | 0.800 |

Paired test-tile bootstrap interval for exponential minus linear RMSE: -2.218 to 3.825 km/s. Negative favors the exponential benchmark. This interval is descriptive on reused data, not a new blind test.
The exponential law approaches the linear control at these short distances. The test difference interval includes zero, so this diagnostic does not resolve their curvature. Nominal 95% coverage is only 80% on the 25 reused test objects; the fixed 300 km/s error model is not a demonstrated adequate description. A distance-redshift trend cannot identify companions, their identity, or a first-principles cause. Catalog frame corrections, shared distance calibrations, motions and selection are not fully modeled. No rotation, lensing, spectrum or timing observation was fitted here.

## Predictions from the same rate

| Conversion z | Energy transferred | Fixed-condition duration ratio | Static flux / unshifted-source flux |
|---|---:|---:|---:|
| 0.01 | 0.990% | 1 | 0.990099 |
| 0.1 | 9.091% | 1 | 0.909091 |
| 1.0 | 50.000% | 1 | 0.500000 |
| 2.0 | 66.667% | 1 | 0.333333 |

These are scenarios, including high-z extrapolations beyond the fitted nearby sample. The flux formula assumes unchanged photon count, isotropic emission, Euclidean geometry, unchanged clocks and fixed speed. A separate duration stretch by S would reduce flux by another factor S. It cannot be inserted without deriving a timing mechanism.

## Capture, retention and energy

All 18 ODE/matrix-exponential ledger comparisons passed, including decay into an explicit receiving reservoir. They establish consistent bookkeeping for chosen rates, not a microscopic interaction.
Capture probability is 1-exp(-optical_depth). An optical depth of 3 gives about 95% capture, but neither the optical depth nor a capture cross section is derived. If captured power is constant for age T and storage lifetime is t_d, the stored fraction of all captured energy is (1-exp(-T/t_d))/(T/t_d). At T=t_d it is 63.2%; at T=10 t_d it is about 10%.
Conservative gravitational focusing can bend paths toward a well. It is not an absorption law and does not by itself establish a permanent deposited reservoir. A claim that companions seek larger wells needs a ray/transport law and a defined absorbing or bound state. Ordinary freely propagating radiation and an added attractive interaction must not be silently interchanged.

## Microscopic cause and timing requirements

A time-independent resonant three-mode Hamiltonian can exchange one higher-frequency photon for a lower-frequency photon plus a companion without special clock behavior. Twelve explicit matrix-exponential checks conserve energy. This is a limited quantum-mode toy: its conversion is reversible/oscillatory, gives a line mixture, and does not establish a nonzero physical photon-graviton coupling or the observed rate. See inelastic-mode-derivation.md for the derivation and missing steps.
For one initial massless photon, one lower-energy massless photon and one massless companion in otherwise empty space, energy/momentum balance gives |p_initial-p_final|^2 = (Delta E/c)^2 + 2 E_initial E_final (1-cos(theta))/c^2. A single companion with energy Delta E has squared momentum (Delta E/c)^2. Thus this simple channel requires collinear photon emission. This does not establish a nonzero transition rate; a background, recoil target or different dispersion needs separate accounting.
A time-stationary delay kernel has I_out(t)=integral K_R(u) I_in(t-u) du. Shifting an input pulse shifts its output by the same amount. Two otherwise identical input pulses therefore retain their centroid separation, even if each is broadened. The numerical delay checks confirm that statement for an illustrative exponential delay distribution. Constant-speed drift is its zero-width limit.
Consequently, adding the same random travel delays cannot generally stretch every feature of an arbitrary light curve by 1+z. A completed conversion-only explanation needs another derived time-transfer effect, source-population explanation, or explicitly time-dependent/nonlinear interaction. This is not a requirement to restore special void clocks, nor a proof against every conversion mechanism.
Transient duration and spectral-aging observations remain targets: see [Blondin et al.](https://arxiv.org/abs/0804.3595) and [DES time-dilation analysis](https://arxiv.org/abs/2406.05050). Their expansion interpretation is not adopted here; their measurement reductions must be audited before a candidate-specific likelihood. No time-dilation data fit was performed in this pass.

## Work-through of the eight checklist stages

1. Redshift forward model: specified in protocol.json, with fixed distances and constant-rate benchmark.
2. Spectra/counts/energy/timing: conditional predictions and mathematical checks computed; physical interaction missing.
3. Training fit: executed on recovered training rows with frozen protocol and descriptive bootstrap.
4. New validation: not available; historical validation/test rows evaluated and explicitly labeled reused.
5. Source map: per-photon conversion computed; joint luminosity histories and corrected environments missing.
6. Capture/storage: rate/lifetime scenarios checked; spatial cross sections and supported state missing.
7. Motion/lensing: prior conditional diagnostic retained; no unified response derived in this pass.
8. Full register: all 32 requirements mapped in coverage.json; none newly claimed complete at full scope.

## Next actionable work

Construct a candidate inelastic interaction with its receiving/recoil sector, and derive its joint energy, direction, polarization and arrival-time kernel. Use the fitted loss coefficient as a target rather than a derived constant. In parallel in the research sequence, obtain unexposed spectra and transient time-series with measurement provenance. A supported deposit law and joint gravity/light-propagation action are needed before fitting extra gravity. No assumption of expansion, Big Bang or an independent dark-matter source is required by this work plan.

Inputs, predictions, metrics and checks are in results.json, predictions.csv, checks.json and coverage.json alongside this report. The runnable source and protocol are preserved in the repository.
