# Restoring-field robustness result

Adding an everywhere-positive quadratic restoring potential to the homogeneous finite-radiation model prevents persistent one-way rolling. The field returns energy to the photon bath and the reference redshift rate changes sign. This exposes an assumption behind the earlier positive-rate solution: its field potential was exactly flat.

| Restoring frequency m | Maximum optical factor n | First return time |
|---|---:|---:|
| 0.01 | 8.26209 | 367.39167 |
| 0.03 | 3.12996 | 143.57027 |
| 0.1 | 1.42195 | 54.02145 |

All numbers are dimensionless. K=g=c0=1 and initial radiation density is 0.003. The m=0 control keeps rolling during the same maximum integration window. This is not a fit to galaxy data or a physical timescale constraint.

## Checks

Both numerical tolerances agree within 5.7e-13 in state over the sampled intervals. Relative energy errors are below 1.1e-14. The calculated first maxima agree with the independently derived energy turning points; the first return takes twice the first ascent time. Each positive-m case restores the initial radiation energy at return and exhibits both positive and negative rates. The script uses an additional maximum step restriction, so unusually small tolerance differences are consistent with both runs resolving the smooth problem at essentially the same step size. Energy and analytic checks provide independent validation of the integration.

Reproduce with `python research_work/results/distributed-restoring/run.py`. Protocol was written before execution. The established potential form, established energy method and internally derived expressions are distinguished in derivation.md; originality is unverified.

## Next implication

The next spatial model must specify whether restoring forces operate only in dense regions and how that dependence follows from matter coupling. It must allow radiation and field energy to move, rather than assuming a uniform bath. This test does not rule out a weakly restored void sector, continuous source-driven evolution or other interactions, but it prevents promoting the zero-potential rolling case to a generic explanation. Atomic clocks, field normalization, fresh observational validation and the full research goal remain unresolved.
