# Constant opacity revision: declared comparison

Before executing this revision, use the existing Pantheon+ standardized magnitudes and full STAT+SYS covariance. All rows were previously exposed by shared_interaction_test; no blind claim. Preserve final reserved datasets elsewhere.

Fix alpha=0.0002488993286382367/Mpc from the existing nearby-galaxy fit and b=1 from the timing requirement. Add the optional, unadopted whole-photon loss law P=exp(-epsilon alpha D), epsilon>=0. Known flux accounting then gives F=L P/[4 pi D^2 S^2]. This changes photon-number conservation and must be reported as a revision, not the previous model.

Calibrate absolute magnitude only on IS_CALIBRATOR=1, treating CEPH_DIST as stipulated geometric distance moduli. Include the proposed propagation attenuation at those distances, using S=exp(alpha D) rather than assigning all their measured local redshift to conversion. Fit the single epsilon using noncalibrator 0.1<=zHD<0.3 rows. Freeze it and evaluate zHD>=0.3. Retain the old conditional frame convention: D=ln(1+zHD)/alpha and S=1+zHEL for high-redshift flux, with the same latter factor for P. No expansion distance is used, but released velocity/standardization corrections retain model dependence.

Use generalized least squares and full calibration/evaluation cross covariance. Propagate the fitted epsilon's linear sampling covariance into the more distant residuals, including cross terms. Also report baseline scores with epsilon=0 and mean residuals in the existing redshift bins. Do not fit per-bin corrections or select another split after results. Quoted uncertainties remain conditional on the published reduction and fixed stipulated distances/rate; this is not an independent raw-photometry or bolometric measurement analysis.
