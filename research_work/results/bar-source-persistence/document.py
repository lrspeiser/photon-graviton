"""Render the retained numerical summary into the checkpoint report."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'summary.json').read_text())
s=json.loads((HERE/'spherical-residence.json').read_text())
lines=['## Results','',
       'Percentages below refer to a continuously supplied population at approximately 244 Myr. The source on each launch sphere is normalized independently.','',
       '| Launch radius | Model/grid | Ever visited center (%) | Currently inside (%) |',
       '|---|---|---:|---:|']
for R in (1.,3.):
    control=next(a for a in s if a['R']==R)['refined']['statistics'][-1]
    lines.append(f"| {R:g} kpc | Spherical control | {100*control['previously_entered_fraction']:.3f} | {100*control['central_residence_fraction']:.4f} |")
    for grid,L in [('coarse',64),('fine',64),('fine',40)]:
        a=next(a for a in d['rows'] if a['R_kpc']==R and a['T']==.25 and a['grid']==grid and a['L']==L)
        lines.append(f"| {R:g} kpc | Bar {grid}, L{L} | {100*a['continuous_ever_entered']:.3f} | {100*a['continuous_inside']:.4f} |")
failed=[g for g in d['statistical_gates'] if not g['passes']]
lines+=['',f"All {len(d['orbit_checks'])} final orbit-tolerance gates pass. {len(failed)} of {len(d['statistical_gates'])} source/field statistical gates fail; they are retained below. Passing the absolute residence threshold of 0.5 percentage points does **not** establish small relative error: the 1 kpc coarse/fine residence estimates differ by about a factor of 2.5 at the final epoch.",
        '', '| Failed comparison | Radius (kpc) | T | Statistic | Absolute difference | Limit |',
        '|---|---:|---:|---|---:|---:|']
for g in failed:
    lines.append(f"| {g['kind']} | {g['R_kpc']:g} | {g['T']:g} | {g['metric']} | {g['difference']:.6f} | {g['limit']:g} |")
lines+=['',f"The corrected audit records {len(d['naive_crossing_disagreements'])} orbit/tolerance instances where naive boundary-event counts differ from turning-point-bracketed crossings. These are numerical detection differences, not separate physical events or failed final tolerance gates. See the exact records in summary.json.",
        '', 'First entry after an earlier radial turn occurs in the retained data (the weighted fractions are in summary.csv). Thus the previous first-pass avoidance result cannot establish lasting avoidance. At both launch radii the sampled bar residence fractions are below the spherical control, but the actual formed population remains uncomputed.']
p=HERE/'report.md';text=p.read_text(encoding='utf8')
if '<!-- RESULTS -->' not in text:
    start=text.index('## Results');end=text.index('## Meaning and next goal')
    text=text[:start]+'<!-- RESULTS -->\n\n'+text[end:]
p.write_text(text.replace('<!-- RESULTS -->','\n'.join(lines)),encoding='utf8',newline='\n')
