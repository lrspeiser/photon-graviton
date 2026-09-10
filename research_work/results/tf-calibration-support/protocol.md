# Parent-catalog calibration support

Apply the original measurement-availability, inclination and width-error cuts to all 10,737 staged rows. Compare eligible parent measurements with the 73 indicator calibrators, using corrected apparent magnitude m=icmag-Ai, reconstructed x=log10(Wmx/sin Inc)-2.5 and inclination. No predicted distances, redshifts or outcome-based cuts enter this diagnostic.

Before computing support, declare: report min/max and 5/50/95 percentiles; count parent rows outside each calibrator range; use the convex hull in (x,m), then (x,m,Inc), only as an empirical interpolation diagnostic. Being inside the hull does not establish complete coverage, unbiased selection, physical suitability or a calibrated distance posterior. Being outside in apparent magnitude need not invalidate a universal TF luminosity law, but makes measurement/selection-error transport less directly checked.

Report fixed 2-mag apparent-magnitude bins and counts of calibrator representation. Retain all exclusions with reasons. No inverse-probability weights will be inferred from 73/(parent count): this is a crossmatched convenience calibration sample and the relevant sampling probabilities are unknown. Do not mistake absence of calibrators for a physical void.

Verify all calibrators reproduce original x and m to numerical precision, every calibrator lies in its own hull, and all 10,737 rows have a unique documented disposition. Archive the parent support flags, aggregate results and source hashes. This addresses transportability of the error calibration, not the physical redshift mechanism or a new holdout.
