# Eight-region cross-validation robustness

The same three predeclared candidates were fitted while excluding each of eight sky regions (four RA quadrants crossed with two declination hemispheres). Region populations are 21,19,16,17,20,31,27,13. Every original sky tile remains inside a single excluded region and all 164 groups are retained. No new parameter, nuisance correction or model selection was introduced.

Out-of-fold RMS errors in c times redshift units are 450.36 km/s for the linear control, 450.72 for constant exponential and 461.63 for the smooth-rate candidate. The smooth-minus-constant descriptive paired-region bootstrap interval is [4.17,17.13] km/s; the linear-minus-constant interval is [-2.57,1.79]. With only eight exposed blocks and shared calibration, these are robustness summaries rather than independent significance or coverage guarantees.

The added distance flexibility again worsens prediction in this comparison. Constant-rate fitted c*alpha ranges from 72.77 to 74.23 across excluded regions; smooth endpoint coefficients range from 71.07-77.67 and 64.32-76.74. Full fold ranges, biases, absolute errors and training-distance extrapolation counts are saved. Formula structures are established mathematics or empirical diagnostics; no microscopic rate was derived by fitting them.

Reproduce with `python research_work/results/redshift-priority/coarse_sky.py`. The protocol predates execution. No preferred more-flexible law is selected, and fresh sample quarantine remains incomplete.
