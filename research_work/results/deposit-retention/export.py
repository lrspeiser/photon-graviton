from pathlib import Path
import json,csv

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'results.json').read_text())
assert all(x['passes'] for x in data)
rows=[]
for item in data:
    x=item.get('refined_again',item['refined'])
    rows.append(dict(kappa0=x['kappa0'],mass=x['mass'],K=x['K_circular'],
                     U_self=x['U_self'],U_external=x['U_external'],
                     mechanical_energy=x['K_circular']+x['U_self']+x['U_external'],
                     core_bound_speed_ratio=x['isotropic_mean_speed_squared_over_escape_squared'][0],
                     outward_density_increase=bool(x['outward_density_increases']),
                     refinement_change=item['max_relative_change']))
with (HERE/'comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
lines=['## Final checkpoint','',
'| kappa0 | Deposited mass | Required orbital K | U_self | U_external | Mean squared speed / escape squared at r0.1, formal isotropic closure |',
'|---|---:|---:|---:|---:|---:|']
for r in rows:
    lines.append(f"| {r['kappa0']:g} | {r['mass']:.7f} | {r['K']:.7f} | {r['U_self']:.7f} | {r['U_external']:.7f} | {r['core_bound_speed_ratio']:.5g} |")
lines+=['',
'The high-opacity core ratio is approximately6512, whereas any isotropic population consisting entirely of bound particles requires a ratio no larger than1. Its density increases by approximately6803-fold between r0.1 and r3. These discrepancies are much larger than the numerical changes. The lower-opacity profiles pass the sampled necessary inequalities; their complete distribution functions have not been validated.',
'',
f"All three final512/1024-shell comparisons pass the unchanged1-percent gate; the largest difference is {100*max(r['refinement_change'] for r in rows):.4f} percent. The independent kinetic-energy integrals agree to floating-point precision. This tests the support quadrature and virial identity, not physical formation or collective stability.",
'',
'Specific circular angular momenta and speeds at r0.1,0.3,1,3,10 are retained in results.json. [comparison.csv](comparison.csv) provides energy and necessary-condition diagnostics; NPZ files retain the generated profiles. Mechanical energies are in G M_b^2/a units and cannot be read as a photon-supply budget without the missing physical scale and assembly history.','']
p=HERE/'report.md';s=p.read_text(encoding='utf8').split('## Final checkpoint')[0]
p.write_text(s+'\n'+'\n'.join(lines),encoding='utf8',newline='\n')
print('\n'.join(lines))
