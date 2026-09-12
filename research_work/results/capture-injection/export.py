from pathlib import Path
import json,csv

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'results.json').read_text())
assert all(x['passes'] for x in data['checks'])
rows=[]
for pair in data['results']:
    x=pair['refined']
    for b in x['binding']:
        rows.append(dict(kappa0=x['kappa0'],c_over_v0=b['c_over_v0'],bound_energy_fraction=b['bound_energy_fraction'],
                         bound_rest_fraction=b['bound_rest_fraction'],rest_energy_fraction=x['rest_energy_fraction'],
                         bulk_kinetic_energy_fraction=x['kinetic_energy_fraction'],central_unresolved_source_fraction=x['unresolved_central_source_fraction']))
with (HERE/'comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
lines=['## Final checkpoint','',
'| kappa0 | Initially bound source energy, c/v0=300 | c/v0=1000 | c/v0=3000 | Source energy converted to bulk kinetic energy |',
'|---|---:|---:|---:|---:|']
for pair in data['results']:
    x=pair['refined'];f=[100*b['bound_energy_fraction'] for b in x['binding']]
    lines.append(f"| {x['kappa0']:g} | {f[0]:.6g}% | {f[1]:.6g}% | {f[2]:.6g}% | {100*x['kinetic_energy_fraction']:.6g}% |")
checks=data['checks']
lines+=['',
f"All declared gates pass. The largest change in bound-energy fraction across radial resolutions is {100*max(max(c['bound_fraction_changes']) for c in checks):.4f} percentage points. The largest absolute beta disagreement between continuity and direct angular integration is {max(c['continuity_error'] for c in checks):.5g}; angular refinement changes beta by at most {max(c['angular_error'] for c in checks):.5g}. Local rest-plus-kinetic energy and momentum identities pass to floating-point tolerance.",
'',
'At c/v0=1000, the middle-opacity example puts only about0.152 percent of its instantaneous captured-energy source into initially bound cold cohorts under this rule. In the weakest-opacity example all sampled source regions pass that initial binding check. Both still inject radial flow rather than circular support.',
'',
'Zero entries retain the finite-grid qualification above. No source luminosity or physical age is inferred. [comparison.csv](comparison.csv) includes rest-mass fractions and the unresolved central-shell source bounds.','']
p=HERE/'report.md';s=p.read_text(encoding='utf8').split('## Final checkpoint')[0]
p.write_text(s+'\n'+'\n'.join(lines),encoding='utf8',newline='\n')
print('\n'.join(lines))
