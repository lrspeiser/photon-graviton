# Coverage, missing inputs and next outcomes

This is a comprehensive register of the project's declared requirements, not a claim to have calculated every observation. Status words distinguish measured-summary comparisons, synthetic mathematics, inherited constraints and absent predictions. None of the model branches is globally validated.

## Electromagnetic coverage

| Band/observable | Available evidence used now | What is still needed |
|---|---|---|
| Radio line centroids and widths | Methanol: two independent centroid groups, tied third line, fourth excluded for profile mismatch; six maser-distance summaries | Raw spectra/visibility products, shared-component fits, covariances, laboratory-frequency error propagation and paired optical/UV sightlines |
| Radio pulse timing | Published timing framework retained; no new pulsar fit | Arrival-time files and timing models, distances/proper motions, dispersion and intrinsic spin/acceleration priors; predicted time-field signature distinguishable from fitted spin |
| Microwave monopole | All 43 FIRAS residual channels recomputed | Full covariance, calibration priors and foreground likelihood; physical source/history and photon-number law |
| Microwave temperature/polarization sky | Existing 83 TT and 66 EE binned diagnostics retained, not rerun as a new fit | Perturbation/source model and angular transfer functions, actual likelihood; no imported standard-model best-fit sky as a prediction |
| Infrared continuum and lines | Synthetic propagation; prior DustPedia/source calculations retained | Calibrated same-source spectra, dust extinction/emission, passbands and source model |
| Optical galaxy redshifts | 164 group predictions with fixed distances and original exposure status | Individual spectral provenance, velocity/covariance model, truly new systems; 24 ELVES candidates remain pending, not validated |
| Optical supernova spectral evolution | All 35 published aging-rate summaries checked and recomputed | Original spectra, independent source-template calibration, shared uncertainties, selection and source-evolution model |
| Optical transient photometry | DES observer-time cache audited previously; inference pilot remains failed | Reliable marginalized likelihood, injection recovery and coverage, selection and bandpass treatment before fitting observed curves |
| Ultraviolet absorption | Formula predictions only in this pass | H2/metal/Lyman spectra and laboratory lines, common component identity, wavelength calibration and intergalactic absorption model |
| X-ray spectra | Formula predictions only in this pass | Calibrated line/continuum observations with response matrices; ionization, broadening, absorption and emitted spectrum model |
| Gamma-ray timing/spectra | GRB090510 observed summary; GW170817 comparison | Photon event files, detector response, independently justified intrinsic lag/energy model, high-energy attenuation; host redshift is not a gamma-ray rest-frequency measurement |
| Local redshift and null tests | Clock summaries, path-length predictions and existing atomic-response failures | Actual geometry/clock/matter completion; resonators, optical and microwave links, solar spectra and finite-exposure pipeline |
| All-band image and polarization fidelity | Independent-jump linewidth allowance and nondispersive ansatz | Derived scattering angle, frequency diffusion, polarization/birefringence and finite packet transfer kernel; images/interferometry likelihood |

The seven representative frequencies in seven-band-predictions.csv are not astronomical sources. One constant fractional rate is mathematically defined at every positive frequency; its physical applicability must still be demonstrated, including geometric-optics breakdown and material interactions. Missing UV/X-ray or paired broad-band data cannot be turned into observational passes by drawing a smooth curve.

## All original requirements retained

| ID | Requirement | Current status and concrete next outcome |
|---|---|---|
| R01 | Operational definitions | Source/receiver clocks, physical path, energy and event stretch distinguished; build one global field and measurement convention |
| R02 | Microscopic conversion | Rate postulate only; derive positive, gauge-consistent interaction and inverse channels |
| R03 | Conservation | Scalar accounts verified; derive stress-energy exchange, recoil and driver work in completed theory |
| R04 | Quantum/classical consistency | Not established for T1; count modes and test stability, causality and interactions |
| R05 | Companion identity | Unresolved; same-metric massless waves contradict loss-free travel in T1; specify field and dispersion |
| R06 | Source production | Historical catalog work retained; derive band-dependent production and source histories |
| R07 | Redshift | 164+6 exposed predictions and radio contrast executed; independent distance/spectral validation pending |
| R08 | Event timing | 35 published summaries and 60 prescribed-field cases; independent full transient inference pending |
| R09 | Brightness/counts | Flux Jacobian derived; photon-removal brightness penalty explicit; no joint photometric fit |
| R10 | Angular distance/surface brightness | Flat diagnostic only; derive ray-bundle geometry and rulers under candidate |
| R11 | Spectral/image fidelity | Radio chromaticity and limited jump-noise check; UV/X-ray, polarization and image kernel missing |
| R12 | Clocks/cavities | Conditional summaries and inherited atomic failures; full matter coupling and apparatus models missing |
| R13 | Local gravity | Ordinary baseline retained; T1 transition forces and complete solar/pulsar predictions missing |
| R14 | Capture/retention | Desired postulates and historical toy models; no T1 capture cross-section/lifetime |
| R15 | Halo profiles | No newly predicted deposit distribution; derive supported extended states rather than impose target map |
| R16 | Galaxy rotation/scaling | Historical 149-galaxy/3150-point benchmark retained; T1 supplies no source-derived force curve yet |
| R17 | Other galaxy dynamics | No completed predictions for streams, dwarfs or vertical gravity; derive same gravitational response |
| R18 | Lensing | No joined T1 metric/stress prediction; calculate both lensing and motion from same fields |
| R19 | Clusters | Prior conditional empirical comparisons retained; no new cluster-only parameter allowed |
| R20 | Merging systems | No deposit dynamics; calculate stripping, lag and offsets |
| R21 | Environment | Illustrative kappa profiles only; reconstruct independent environmental input and test held-out paths |
| R22 | Stellar/thermal physics | Extra sink heating exposed, not solved; compute source/receiver heat and cooling |
| R23 | Gravitational waves | Relative-speed summary check only; derive ordinary tensor and companion propagation/production separately |
| R24 | Compact objects | No completed capture/ringdown/accretion predictions |
| R25 | CMB spectrum | 43-channel conditional comparison executed; source and photon-number/thermalizer completion absent |
| R26 | CMB angular structure | Earlier bandpowers retained; no anisotropy/polarization generation model |
| R27 | Structure/BAO | No nonexpanding growth/history solution or independent ruler |
| R28 | Abundances/chronology | No complete abundance or stellar-age history; T1 past boundary explicitly recorded |
| R29 | Background/topology/drift | Exact kinematic zero-drift construction, finite domain; not a self-consistent background |
| R30 | Radiation/entropy | Bath normalization and sink energy calculated; no global radiation equilibrium or fuel history |
| R31 | Identifiability | Distance, source-evolution and covariance caveats explicit; no combined all-data probability |
| R32 | Distinguishing predictions | Candidate differences in timing, chromaticity, flux, drift and companion loss are explicit; freeze new tests after mechanism derivation |

## Data acquisition order

1. Complete a common dynamical candidate first enough to predict an observable; otherwise more data cannot identify an unspecified coefficient for every effect.
2. Repair the existing DES inference gate and validate it on synthetic injections before using real fluxes. Preserve failed thresholds.
3. Acquire independent same-component radio-to-optical/UV spectra with covariance; keep discrepant-profile lines as diagnostics rather than silently discarding them.
4. Acquire calibrated optical/IR source photometry and full FIRAS likelihood inputs to test one photon-number/sink law jointly.
5. Add UV, X-ray and gamma-ray response products plus local clock/pulsar arrival data. Choose sample and source-lag/atomic assumptions before comparing outcomes.
6. Resume the linked gravity/formation calculation once companion coupling and capture are specified. The earlier deferred supply question is not resolved by fitting redshift.

No file has been identified as sufficient by itself to close these physics gaps. Several are missing physical equations, not missing downloads. The active research goal remains incomplete.
