# Exact one-third lens transfer and paired galaxy assessment

The exact-one-third law has now been run on the six existing lens systems, with all galaxy capture constants and optical parameters frozen. Stellar mass and orbital anisotropy are fitted only to inner stellar measurements, as before. This exposed-sample test uses conditional luminosity and morphology proxies.

| Population proxy | Outer motion residual-square sum | Lens-angle fractional RMS | Orbit boundary count |
|---|---:|---:|---:|
| attenuated_Chabrier | 44.32216 | 13.92308% | 0 |
| attenuated_Salpeter | 42.70398 | 14.09369% | 0 |

Exact one-third does not repair the existing lens mismatch. Its 13.923/14.094% errors remain close to the fitted-exponent values 13.938/14.110% and above the earlier intercepted-source value 13.112%. Better outer-motion errors than the older source model do not establish joint agreement.

## Paired galaxy results

| Split | Error measure | Galaxies | Companion closer than MOND | Companion closer than baryons |
|---|---|---:|---:|---:|
| validation | rms_kms | 29 | 11 | 26 |
| validation | log_rms | 29 | 12 | 29 |
| test | rms_kms | 31 | 9 | 30 |
| test | log_rms | 31 | 13 | 31 |

These counts describe how often the candidate helps, complementing the previously reported aggregate RMS. They are not discovery probabilities, and a win count does not weight the size of a miss. Do not select only winning galaxies or retune on comparison targets. Per-galaxy errors are retained in cross-test-audit-results.json.

Verification: exact q=1/3; matching capture-file hash; unchanged optical and photometric hashes; identical observed stellar velocities and geometry; probabilities inside (0,1); at least one successful optimizer start for every system/population.

See third-retention-lensing-protocol.md for the pre-run specification. No cluster-transfer or mechanism goal is completed.
