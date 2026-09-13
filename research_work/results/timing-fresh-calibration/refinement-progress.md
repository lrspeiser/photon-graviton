# Numerical refinement: first completed case

The previous turn added an input/result integrity audit while the verified live refinement process ran. The current turn obtains the first completed case and audits its numerical result. This is partial evidence, not a completed batch.

Case split_gaussian-snr5-b0-seed1913 was selected by the inherited worst-error rule. Increasing shape quadrature from 2048 to 8192 samples with a different scramble changes the fitted exponent by 0.00197427 (limit 0.05). The posterior-weighted mean centered likelihood-curve difference is 0.00412085 (limit 0.1). Both frozen numerical checks pass for this case.

The auditor verifies input/code hashes, base and refined arrays, the selected case, parameter bounds and support, and reproduces the saved best objective with zero discrepancy. The historical one-case audit snapshot is preserved in refinement-first-case-audit.json. It is not process-liveness evidence and need not represent later run progress. The active process was separately polled through its live execution handle.

Seven selected cases remain at this checkpoint. A single selected passing case establishes neither complete numerical convergence nor coverage. The width grid is unchanged, actual source evolution and brightness selection are untested here, and the joint physical redshift-duration-brightness mechanism remains unresolved. No original calibration cases are replaced by this refinement; no observed flux or final holdout was used. All six objectives remain open.
