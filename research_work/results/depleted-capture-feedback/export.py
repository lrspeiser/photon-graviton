"""Export every case, explicitly retaining unresolved refinement status."""
import csv
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
base=json.loads((HERE/'results.json').read_text())
extra=[]
for R in (30,100):
    path=HERE/f'refinement-{R}.json'
    if path.exists():
        extra.extend(json.loads(path.read_text()))
rows=[]
for pair,check in zip(base['results'],base['refinement_checks']):
    result=pair['refined'];passed=check['passes'];change=check['max_relative_change']
    matches=[x for x in extra if all(x[k]==result[k] for k in ('R','p','C','kappa0'))]
    if matches:
        last=matches[0]['trials'][-1]
        result=last['result'];passed=last['passes'];change=max(last['relative_changes'])
    rows.append({k:result[k] for k in ('R','p','C','kappa0','B','shells','mass','central_depth','capture_fraction','energy_ledger_error','chord_mass_error')}
                |dict(vc2_r1=result['vc2'][0],vc2_r3=result['vc2'][1],vc2_r10=result['vc2'][2],
                      numerical_gate_passed=passed,max_relative_change=change))
with (HERE/'comparison.csv').open('w',encoding='utf8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
print(json.dumps(rows,indent=2))
