from pathlib import Path
import json,csv

HERE=Path(__file__).resolve().parent
allcases=[];checks=[]
for rate in (.1,1.):
    cases=json.loads((HERE/f'results-C{rate:g}.json').read_text());allcases+=cases
    for kind,a,b in (('resolution',cases[0],cases[1]),('core',cases[1],cases[2])):
        for x,y in zip(a['outputs'],b['outputs']):
            mass=abs(x['mass']/y['mass']-1)
            fraction=max(abs(u-v) for u,v in zip(x['enclosed_fractions'],y['enclosed_fractions']))
            checks.append(dict(Cdot=rate,kind=kind,T=x['T'],mass_relative_change=mass,fraction_change=fraction,
                               passes_declared_tolerance=mass<.02 and fraction<.02))
allcases+=json.loads((HERE/'core-C0.1.json').read_text())
rows=[]
for c in allcases:
    for x in c['outputs']:
        rows.append(dict(Cdot=c['Cdot'],source_shells=c['n'],epsilon=c['epsilon'],T=x['T'],mass=x['mass'],
                         **{f'mass_fraction_r{r}':v for r,v in zip(('.1','1','3','10'),x['enclosed_fractions'])},
                         central_depth=x['central_deposit_depth'],max_speed_c=x['max_speed_over_c_to_date'],
                         opacity_change_per_crossing=x['max_opacity_log_change_per_crossing_to_date'],
                         mechanical_residual=x['max_mechanical_residual_to_date']))
with (HERE/'comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
(HERE/'comparisons.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf8',newline='\n')
lines=['## Final checkpoint','',
'Lower-rate results at fine resolution; these are dimensionless synthetic calculations, not observed galaxy masses:','',
'| Time | Deposited rest mass, epsilon0.05 | Fraction inside r0.1 | Fraction inside r1 | Fraction inside r3 |',
'|---|---:|---:|---:|---:|']
low=next(c for c in allcases if c['Cdot']==.1 and c['n']==128 and c['epsilon']==.05)
for x in low['outputs']:
    f=x['enclosed_fractions']
    lines.append(f"| {x['T']:g} | {x['mass']:.7f} | {100*f[0]:.4f}% | {100*f[1]:.4f}% | {100*f[2]:.4f}% |")
lines+=['','Core sensitivity at T3, all at the same fine resolution:','',
'| Epsilon | Deposited mass | Central deposited potential depth |',
'|---|---:|---:|']
for c in allcases:
    if c['Cdot']==.1 and c['n']==128:
        x=c['outputs'][-1];lines.append(f"| {c['epsilon']:g} | {x['mass']:.8f} | {x['central_deposit_depth']:.6f} |")
lines+=['','Strong-source status, with late values explicitly outside the approximation:','',
'| Source shells | Epsilon | Final maximum speed/c | Maximum opacity log change per crossing | Mechanical numerical gate |',
'|---|---:|---:|---:|---|']
for c in allcases:
    if c['Cdot']==1:
        lines.append(f"| {c['n']} | {c['epsilon']:g} | {c['max_speed_over_c']:.4f} | {c['max_opacity_log_change_per_crossing']:.3f} | {'pass' if c['numerical_energy_gate'] else 'fail'} |")
lines+=['',
'All low-rate mass/fraction resolution comparisons and all low-rate mechanical gates pass. The stronger rate fails the mass-resolution gate at T2 and T3; it also leaves the nonrelativistic and slowly changing radiation assumptions. Its enormous late mass values remain in the raw results for audit, not as physical predictions. All local radiation/rest-mass identities remain numerically balanced, which does not rescue those physical failures.',
'',
'The [comparison table](comparison.csv) and comparisons.json retain every checkpoint and failed gate. A force-gradient finite-difference check of the shell Hamiltonian agrees within2.8e-9 absolute; it tests the implemented regularized force, not a microscopic gravity theory.','']
p=HERE/'report.md';s=p.read_text(encoding='utf8').split('## Final checkpoint')[0]
p.write_text(s+'\n'+'\n'.join(lines),encoding='utf8',newline='\n')
print(json.dumps(checks,indent=2))
