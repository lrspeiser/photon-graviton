from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'summary.json').read_text())
failed=[r for r in d['comparisons'] if not r['passes']]
lines=['# Source-direction refinement','',
'The source law, radial launch speed, fixed ordinary bar and integration horizon are unchanged. This compares deterministic angular quadratures, not different physical theories or fresh observational fits. Each launch sphere is normalized separately.','',
'## Population fractions at approximately 244 Myr','',
'| Launch radius (kpc) | Angular grid | Cohort ever entered (%) | Continuous population ever entered (%) | Continuous population inside now (%) |',
'|---:|---|---:|---:|---:|']
for r in d['rows']:
    if r['T']==.25:
        lines.append(f"| {r['R']:g} | {r['n']}x{2*r['n']} | {100*r['cohort_ever_entered']:.3f} | {100*r['continuous_ever_entered']:.3f} | {100*r['continuous_inside']:.4f} |")
lines+=['',f"All {len(d['checks'])} new orbit status checks pass: **{d['all_orbit_checks_pass']}**. Every new orbit is checked for Jacobi drift. Tighter-tolerance pairs cover every eighth direction (and any invariant failure), not every new orbit. The two grids add64 and144 trajectories across the two spheres. L40 is used throughout; prior L40/L64 agreement on the6x12 grid is not a new field-order test for every trajectory.",
'',f"{len(failed)} of {len(d['comparisons'])} adjacent-grid statistical comparisons fail the predeclared development gates. The complete differences, including relative residence changes, are in summary.json. These are not confidence intervals.",
'','| Comparison | R (kpc) | T | Failed statistic | Absolute change | Limit |',
'|---|---:|---:|---|---:|---:|']
for r in failed:
    lines.append(f"| {r['grid_a']} to {r['grid_b']} | {r['R']:g} | {r['T']:g} | {r['metric']} | {r['absolute_change']:.6f} | {r['limit']:g} |")
lines+=['','## What this means','',
'A central encounter flag jumps at a grazing orbit, whereas time inside tends to zero. Consequently the encounter statistic can be harder to integrate over directions. Both matter for describing trajectories, but instantaneous gravity follows the actual spatial distribution, not the fraction with a past visit. A pass at the loose0.005 absolute residence threshold does not establish a small relative uncertainty.','',
'At the final epoch,8-to12 residence changes are19.6 percent for the1kpc source and33.2 percent for the3kpc source relative to the finer values. In contrast, their total angular source normalizations change only0.0000904 percent and0.277 percent. Thus source normalization can be well integrated while a trajectory-dependent statistic remains poorly resolved. Do not repair this by tuning the injection amplitude.','',
'The [population-to-gravity derivation](population-to-gravity.md) uses that spatial distribution directly. It applies known transport and Newtonian formulas under the declared cold-particle postulate; none of those standard formulas is claimed unique. Source amplitude and physical capture funding remain unspecified. The companion/time redshift mechanism is not tested here.','',
'Next assess the three-dimensional force coefficients under direction refinement. Even agreement of enclosed central mass does not establish agreement of local forces in a nonspherical distribution. A self-gravitating formation calculation, full ordinary-matter field and joint observational tests remain required. All nine goals remain active.','',
'Reproduce with run.py8 and run.py12 (pass the number as a separate argument), followed by export.py and document.py. Raw data retain primary and tighter-tolerance trajectories. No discarded failing direction, tuned source parameter or opened holdout is involved.']
(HERE/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8',newline='\n')
