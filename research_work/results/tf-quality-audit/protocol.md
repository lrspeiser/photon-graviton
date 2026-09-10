# Input quality audit of the exposed TF calibration

Use the same 73 indicator calibrators and their already exposed, group-excluded predictions. Before inspecting joined outcomes, declare image-quality strata QSflag=0..5, and summaries for QSflag>=4 and QSflag=5. These are diagnostics, not a newly validated exclusion rule. Do not clip residuals, change alpha, refit the luminosity relation, or replace any distance.

Retrieve only PGC, QSflag, QWflag, e_Ai, logWmxi and e_logWmxi from J/ApJ/902/145/table1. Compare published inclination-corrected widths and errors to our reconstruction. The initial formula uses rounded input widths/inclinations, so report discrepancies rather than requiring exact equality. Preserve all joined rows.

For p=1, compare residual magnitude against a diagonal measurement-only scale sqrt(sigma_mu^2+0.05^2+e_Ai^2+a^2*e_logWmxi^2), where the published conservative non-u photometry uncertainty is 0.05 mag and a is the already fitted slope. This omits slope/intercept uncertainty, intrinsic scatter, calibration covariance and selection. It is not a statistical rejection test or full likelihood. Report the median scale and count beyond three such scales, with this limitation explicit.

Check the most discrepant object's physical suitability against independent published kinematic observations. This targeted check is post hoc and cannot certify the remaining sample or justify deleting this object alone. If environmental interactions compromise linewidth distances, record the potential confounding of the intended void test. Do not use a paper's TF-derived distance or redshift-based cluster distance as an independent calibration.
