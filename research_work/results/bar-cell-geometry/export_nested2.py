from pathlib import Path
import json,csv
import numpy as np
HERE=Path(__file__).resolve().parent
records=[];points=None
for layers in (256,512):
    for R in (1,3):
        d=json.loads((HERE/f'volume2-t{layers}-R{R}.json').read_text(encoding='utf8'))
        assert len(d['records'])==3;records.extend(d['records']);points=d['points']
parent=[]
for R in (1,3):
    d=json.loads((HERE/f'volume-t512-R{R}.json').read_text(encoding='utf8'));assert d['points']==points
    parent.extend(d['records'])
comparisons=[]
for b in [r for r in records if r['layers']==512]:
    for kind,candidates in [('age',[r for r in records if r['layers']==256]),('angular',parent)]:
        a=next(r for r in candidates if r['R']==b['R'] and r['T']==b['T'])
        for i,p in enumerate(points):
            dp=abs(a['potential'][i]-b['potential'][i])/max(abs(b['potential'][i]),.01)
            da=float(np.linalg.norm(np.array(a['acceleration'][i])-b['acceleration'][i])/max(np.linalg.norm(b['acceleration'][i]),.01))
            comparisons.append(dict(kind=kind,R=b['R'],T=b['T'],point=p,potential_change=dp,force_change=da,
                                     potential_pass=dp<.02,force_pass=da<.05))
geometry=json.loads((HERE/'nested2-geometry.json').read_text(encoding='utf8'))
old_geometry=json.loads((HERE/'nested-geometry.json').read_text(encoding='utf8'))
rows=[]
for r in records:
    for i,p in enumerate(points):
        rows.append(dict(R=r['R'],T=r['T'],layers=r['layers'],x=p[0],y=p[1],z=p[2],potential_per_GM=r['potential'][i],
                         ax_per_GM=r['acceleration'][i][0],ay_per_GM=r['acceleration'][i][1],az_per_GM=r['acceleration'][i][2]))
with (HERE/'nested2-summary.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
out=dict(records=records,comparisons=comparisons,geometry=geometry,parent_geometry=old_geometry,
         maximum_mass_error=max(abs(r['mass']-1) for r in records))
(HERE/'nested2-summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
for kind in ('age','angular'):
    r=[r for r in comparisons if r['kind']==kind]
    print(kind,'potential fails',sum(not a['potential_pass'] for a in r),'force fails',sum(not a['force_pass'] for a in r),'of',len(r))
print('Final geometry',[r['records'][-1] for r in geometry])
