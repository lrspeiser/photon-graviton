from pathlib import Path
import json,csv

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'results.json').read_text())
assert all(c['passes'] for x in data['checks'] for c in x['comparisons'])
x=data['cases'][-1]
rows=[dict(T=0,kind='source at injection',**{f'fraction_inside_{r}':v for r,v in zip(('.1','1','3','10'),x['source_enclosed_fractions'])},**{f'projected_fraction_inside_{r}':v for r,v in zip(('.1','1','3','10'),x['source_projected_fractions'])},mean_kinetic=x['source_mean_kinetic'],mean_mechanical=x['source_mean_mechanical'])]
for o in x['outputs']:
    rows.append(dict(T=o['T'],kind='moving population',**{f'fraction_inside_{r}':v for r,v in zip(('.1','1','3','10'),o['enclosed_fractions'])},**{f'projected_fraction_inside_{r}':v for r,v in zip(('.1','1','3','10'),o['projected_fractions'])},mean_kinetic=o['mean_kinetic'],mean_mechanical=o['mean_mechanical']))
with (HERE/'comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
lines=['## Final checkpoint','',
'| Population | Inside r0.1 | Inside r1 | Inside r3 | Inside r10 | Mean kinetic energy per added rest mass |',
'|---|---:|---:|---:|---:|---:|']
for r in rows:
    label='Source at injection' if r['T']==0 else f"Moving, T={r['T']:g}"
    f=[100*r[f'fraction_inside_{key}'] for key in ('.1','1','3','10')]
    lines.append(f"| {label} | {f[0]:.5f}% | {f[1]:.5f}% | {f[2]:.5f}% | {f[3]:.5f}% | {r['mean_kinetic']:.7f} |")
ratio=rows[-1]['fraction_inside_.1']/rows[0]['fraction_inside_.1']
lines+=['','The source row gives the spatial baseline for permanent in-place capture and the kinetic energy at injection. It is not a stationary population with that kinetic energy.', '',f'At T30 the innermost enclosed fraction is {ratio:.2f} times the in-place value. The mean mechanical energy remains approximately -0.502400 in all moving-population samples; kinetic growth is balanced by a more negative potential energy. The finite-time response is not a stationary solution.','',
'[comparison.csv](comparison.csv) retains all values. The accumulated mass coefficient at T30 is approximately10.2752 per unit source amplitude; actual added mass is eta times this coefficient. Setting eta=1 and treating this as negligible added mass would be inconsistent with the approximation.','']
lines+=['','Projected enclosed fractions (the input to a conditional spherical weak-lensing calculation):','',
'| Population | Inside b0.1 | Inside b1 | Inside b3 | Inside b10 |',
'|---|---:|---:|---:|---:|']
for row in rows:
    label='Source at injection' if row['T']==0 else f"Moving, T={row['T']:g}"
    f=[100*row[f'projected_fraction_inside_{key}'] for key in ('.1','1','3','10')]
    lines.append(f'| {label} | {f[0]:.5f}% | {f[1]:.5f}% | {f[2]:.5f}% | {f[3]:.5f}% |')
ratio_proj=rows[-1]['projected_fraction_inside_.1']/rows[0]['projected_fraction_inside_.1']
lines+=['',f'At fixed added mass, the innermost projected fraction grows {ratio_proj:.2f}-fold by T30, compared with {ratio:.2f}-fold for true spherical enclosed mass. These are changes in the added component, not multipliers for the total galaxy force or lensing.', '']
p=HERE/'report.md';s=p.read_text(encoding='utf8').split('## Final checkpoint')[0]
p.write_text(s+'\n'+'\n'.join(lines),encoding='utf8',newline='\n')
print('\n'.join(lines))
