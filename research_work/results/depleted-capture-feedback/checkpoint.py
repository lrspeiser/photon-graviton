"""Build the report checkpoint only after both refinement batches finish."""
import csv
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
for R in (30,100):
    assert len(json.loads((HERE/f'refinement-{R}.json').read_text()))==5
rows=list(csv.DictReader((HERE/'comparison.csv').open(encoding='utf8')))
lines=['## Final numerical checkpoint','',
       '| p | C | kappa0 | Mass, R30 | Mass, R100 | Mass ratio | All diagnostic gates, R30 / R100 |',
       '|---|---:|---:|---:|---:|---:|---|']
for p,C in ((4,1.),(6,1.),(6,100.)):
    for k in (.1,10.,1000.):
        a,b=[next(r for r in rows if float(r['R'])==R and int(r['p'])==p and float(r['C'])==C and float(r['kappa0'])==k) for R in (30,100)]
        ma,mb=float(a['mass']),float(b['mass'])
        states=['pass' if r['numerical_gate_passed']=='True' else 'unresolved' for r in (a,b)]
        lines.append(f'| {p} | {C:g} | {k:g} | {ma:.7g} | {mb:.7g} | {mb/ma:.5f} | {states[0]} / {states[1]} |')
unresolved=[r for r in rows if r['numerical_gate_passed']!='True']
ledger=max(float(r['energy_ledger_error']) for r in rows)
chord=max(float(r['chord_mass_error']) for r in rows)
lines+=['',f'{len(rows)-len(unresolved)} of18 cases pass the unchanged numerical gate after refinement. Maximum final incoming/transmitted/deposited ledger relative error is {ledger:.3g}; maximum independent chord-versus-shell absorption discrepancy is {chord:.3g}. These numerical identities do not include missing physical energy sectors.','']
for r in unresolved:
    lines.append(f"Unresolved: R={r['R']}, p={r['p']}, C={r['C']}, kappa0={r['kappa0']}, shells={r['shells']}; maximum relative diagnostic change={float(r['max_relative_change']):.5g}, squared force contribution at r1={float(r['vc2_r1']):.5g}. Do not treat the full case as numerically validated.")
lines+=['','Full values and status labels are in [comparison.csv](comparison.csv); initial and subsequent failed gates remain available in the JSON files. The mass ratios diagnose boundary dependence only over the tested domains; they are not an infinite-domain extrapolation.','']
p=HERE/'report.md';s=p.read_text(encoding='utf8').split('## Final numerical checkpoint')[0]
p.write_text(s+'\n'+'\n'.join(lines),encoding='utf8',newline='\n')
print('\n'.join(lines))
