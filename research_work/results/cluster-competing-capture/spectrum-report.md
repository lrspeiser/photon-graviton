# Spectral check of equal-bolometric histories

The positive emission histories in history-report.md match present photon and companion energy, not necessarily photon spectra. This follow-up assigns every burst the same illustrative 6000 K blackbody spectrum and propagates photon frequency and energy by exp(-alpha c age), conserving photon count. This is a synthetic homogeneous snapshot, not an event-arrival-time calculation, stellar population synthesis or observational fit.

For spectral energy per logarithmic frequency S_emit, the known change of variables gives S_obs(nu)=exp(-x) S_emit(nu exp(x)) with x=alpha c age. Finite bursts retain their earlier width 0.02 in x. The familiar Planck distribution is an assumed source spectrum, not a novel emission law.

Five wavelength intervals partition all frequencies: below 0.4 micrometers, 0.4-0.7, 0.7-5, 5-1000, and above 1000 micrometers. Direct Planck integrals and integration over each burst reproduce both histories' identical total photon energies to relative 1e-9. For the old-burst age x=10 alternative, visible-band energy changes from 0.320416 to 0.430348 common energy units, despite equal bolometric totals. Half the sum of absolute band differences, normalized by the common total, is 0.11130. The x=30 and x=100 alternatives similarly give 0.11127.

Therefore these particular equal-bolometric histories can be distinguished by spectral information under the stipulated identical source spectra. The result does not establish that finite bands uniquely reconstruct arbitrary histories with evolving stellar spectra, dust and other sources. All data are synthetic; no observed background or final holdout is fitted. It identifies multiwavelength information as an additional constraint on old emission, alongside fuel inventories. All six objectives remain open.

Run `python research_work/results/cluster-competing-capture/spectrum.py` after history.py. The input hash and complete band results are preserved in spectrum-results.json.
