# Coarse-region robustness protocol

Pre-run protocol, 9 September 2026. Same 164 exposed groups, adopted distances and three previously fixed candidates. No new model, environment proxy, nuisance motion or optimized partition. Formula provenance remains established functional forms with parameters fitted here; the smooth rate is an empirical diagnostic, not a derived physical law.

Partition deterministically into eight regions: four right-ascension quadrants crossed with north/south declination hemispheres, with boundaries at RA 0/90/180/270/360 and declination 0. Assign an entire existing sky tile to its region; verify no historical tile crosses a region. Exclude each coarse region once and fit only others. No tuning or selection across candidates. Compare pooled and per-region RMS, MAE, median absolute residual and bias. Retain all objects, regardless of region size or residual.

Compute 2000 paired coarse-region residual bootstrap differences versus constant exponential, seed 2026090904. This is descriptive robustness, not independent validation; eight blocks give limited uncertainty resolution and shared calibration remains. Do not compare raw RMS with the fine-tile result as though the underlying evaluation tasks were identical.

Also report training parameter ranges across the eight excluded regions and boundary hits to expose instability. Report actual observed distance coverage of each held-out region and whether its objects exceed the training distance range. No claim of calibrated intervals or a physical precision floor.
