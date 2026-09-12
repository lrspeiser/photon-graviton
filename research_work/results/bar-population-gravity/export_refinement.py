"""Compare population force coefficients at fixed physical postulates."""
from pathlib import Path
import json,csv
import numpy as np
HERE=Path(__file__).resolve().parent
tables=[];raw_checks=[]
for n,name,N in [(6,'orbits.json',36),(8,'orbits8.json',64)]:
    d=json.loads((HERE/name).read_text());assert len(d['records'])==d['expected']==N
    raw_checks.extend(dict(n=n,probe=i+1,passes=r['passes']) for i,r in enumerate(d['records']))
    for R in (1.,3.):
        group=[r for r in d['records'] if r['R']==R];norm=sum(r['raw_weight'] for r in group)
        for j,T in enumerate((.05,.1,.25)):
            p=np.zeros(5);a=np.zeros((5,3))
            for r in group:
                k=r['result']['kernels'][j];v=k['extra_refined'] or k['refined'];w=r['raw_weight']/norm
                p+=w*np.array(v['potential']);a+=w*np.array(v['acceleration'])
            for i,point in enumerate(d['points']):
                tables.append(dict(n=n,R=R,T=T,point=point,potential=float(p[i]),acceleration=a[i].tolist()))
comparisons=[]
for b in [r for r in tables if r['n']==8]:
    a=next(r for r in tables if r['n']==6 and r['R']==b['R'] and r['T']==b['T'] and r['point']==b['point'])
    dp=abs(a['potential']-b['potential'])/max(abs(b['potential']),.01)
    da=float(np.linalg.norm(np.array(a['acceleration'])-b['acceleration'])/max(np.linalg.norm(b['acceleration']),.01))
    comparisons.append(dict(R=b['R'],T=b['T'],point=b['point'],potential_scaled_change=dp,
                            acceleration_scaled_change=da,potential_pass=dp<.02,acceleration_pass=da<.05,
                            coarse=a,refined=b))
out=dict(tables=tables,comparisons=comparisons,orbit_checks=raw_checks,
         all_orbit_checks_pass=all(r['passes'] for r in raw_checks),
         all_source_gates_pass=all(r['potential_pass'] and r['acceleration_pass'] for r in comparisons))
(HERE/'source-refinement.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
csvrows=[dict(n=r['n'],R=r['R'],T=r['T'],x=r['point'][0],y=r['point'][1],z=r['point'][2],
              potential_per_GM=r['potential'],ax_per_GM=r['acceleration'][0],ay_per_GM=r['acceleration'][1],az_per_GM=r['acceleration'][2]) for r in tables]
with (HERE/'source-refinement.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(csvrows[0]),lineterminator='\n');w.writeheader();w.writerows(csvrows)
print(json.dumps(dict(all_orbits_pass=out['all_orbit_checks_pass'],potential_failures=sum(not r['potential_pass'] for r in comparisons),
                     force_failures=sum(not r['acceleration_pass'] for r in comparisons),
                     final_epoch=[r for r in comparisons if r['T']==.25]),indent=2))
