"""Weighted linear-response force coefficients, with explicit numerical status."""
from pathlib import Path
import json,csv
import numpy as np
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'orbits.json').read_text());assert len(d['records'])==d['expected']==36
rows=[];checks=[]
for i,r in enumerate(d['records']):
    for k in r['result']['kernels']:
        checks.append(dict(probe=i+1,T=k['T'],passes=k['final_check']['passes'],
                           initial=k['initial_check'],final=k['final_check']))
for R in (1.,3.):
    group=[r for r in d['records'] if r['R']==R];norm=sum(r['raw_weight'] for r in group)
    for j,T in enumerate((.05,.1,.25)):
        potential=np.zeros(5);force=np.zeros((5,3))
        for r in group:
            k=r['result']['kernels'][j];a=k['extra_refined'] or k['refined'];w=r['raw_weight']/norm
            potential+=w*np.array(a['potential']);force+=w*np.array(a['acceleration'])
        for i,point in enumerate(d['points']):
            rows.append(dict(source_R_kpc=R,T=T,x_kpc=point[0],y_kpc=point[1],z_kpc=point[2],
                             potential_per_GM=potential[i],ax_per_GM=force[i,0],ay_per_GM=force[i,1],az_per_GM=force[i,2]))
out=dict(rows=rows,checks=checks,all_orbits_pass=all(r['passes'] for r in d['records']),
         source_angular_convergence_established=False,units={'potential_per_GM':'1/kpc','acceleration_per_GM':'1/kpc^2'})
(HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
with (HERE/'summary.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps(dict(all_orbits_pass=out['all_orbits_pass'],failed=[r for r in checks if not r['passes']],
                     final_epoch=[r for r in rows if r['T']==.25]),indent=2))
