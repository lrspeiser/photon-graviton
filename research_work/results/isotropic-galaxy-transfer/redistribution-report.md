# Conservative redistribution of the exact-one-third model

All variants preserve the original one-third retention, luminosity proxy and capture constants. Only shared redistribution parameters were fitted on 89 training galaxies. Predictions use the same 29 validation/31 test galaxies and six lens systems; all sets are previously exposed. This completes the bounded family and one boundary follow-up in redistribution-protocol.md.

| Variant | Validation RMS km/s | Test RMS km/s | Validation log RMS | Test log RMS | Chabrier lens RMS | Salpeter lens RMS |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 32.4950 | 23.5910 | 0.115401 | 0.091075 | 13.92308% | 14.09369% |
| shared | 32.4977 | 23.5951 | 0.115398 | 0.091077 | 13.92387% | 14.09453% |
| partial | 31.7489 | 23.1414 | 0.111998 | 0.086648 | 14.13192% | 14.31233% |
| retention_conditioned | 32.5198 | 23.6349 | 0.115360 | 0.091101 | 13.92825% | 14.09844% |
| compact_partial | 31.6661 | 22.5906 | 0.112271 | 0.084544 | 13.91492% | 14.08269% |

## Outcome and iteration

The shared scale remains almost unchanged (s=0.99984), and the retention-conditioned scale also remains close to one. Neither gives a useful improvement. Partial redistribution initially selects f=0.01572 and s=0.25 at its allowed boundary: about 1.6% moves inward, not outward. It improves galaxy errors but worsens lens angles. An outward local optimum exists, but was not the best multi-start fit; it is not the selected result.

A separately declared extension to s>=1/16 finds an interior optimum: f=0.00369849 and s=0.10217473 (approximately). Thus about 0.37% of the existing deposit inventory is concentrated to about one-tenth its original radius. All six optimizer starts converge to essentially the same optimum, without a bound hit. This is an exploratory post-boundary extension and adds two shared parameters; it is not a newly derived law or proof of a central particle population.

The compact follow-up lowers validation/test speed RMS from 32.495/23.591 to 31.666/22.591 km/s. It remains worse than the matched MOND benchmark 26.876/16.398 km/s. Its lens errors remain about 14%; changes at the hundredth-of-a-percentage-point scale do not resolve the lens discrepancy or establish statistical improvement. The inner-star nuisance fits can compensate for redistribution, so lens response need not track a fixed-stellar-mass intuition.

Keep the compact variant as a modest rotation improvement candidate, not a replacement established across all tests. Do not continue adding radial parameters to these exposed objects. The bounded redistribution exercise has found no substantial joint rotation/lensing solution. Coma free-profile fits cannot distinguish a global scale from their already fitted scale, and do not validate this extension.

## Formulas and provenance

eta=X^(1/3)/(1+X^(1/3)) is unchanged. For each baseline deposited density rho0, rho_new(r)=(1-f)*rho0(r)+f*rho0(r/s)/s^3. Consequently M_new(<r)=(1-f)*M0(<r)+f*M0(<r/s). These are known conservative scaling/mixture identities; interpreting them as companion migration is our proposed hypothesis. They preserve total deposit inventory, not automatically total field-plus-kinetic energy. Work, release channels and support still require a physical model. No ages, external supply amplitudes or per-galaxy capture constants were changed.

## Verification

Input hashes, q=1/3, unchanged photometric/optical inputs, observed stellar data, geometry and retention mappings were checked. Every lens fit has a successful optimizer start. Galaxy radial/angular refinement changes aggregate RMS by less than 0.001 km/s; the unchanged baseline reproduces the prior result within 0.05 km/s. Direct finite-volume quadrature independently checks the conservative mixture identity to relative error 3.33e-16. These numerical checks do not establish a physical mechanism or significance of tiny improvements.

Commands: redistribution.py; redistribution.py --compact-followup; lensing.py --redistribution=shared (and partial, retention_conditioned, compact_partial); redistribution-report.py. Every candidate output is preserved separately. All six research goals remain open.
