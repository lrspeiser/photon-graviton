"""Compare adaptive time resolution and fixed-parent-mass spatial refinement."""
from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent
runs=[]
for n in (256,512):
 d=json.loads((H/f'volume-adaptive-t{n}-R3.json').read_text());assert len(d['records'])==1
 r=d['records'][0];assert r['R']==3 and r['T']==.1 and r['layers']==n
 runs.append(r)
points=d['points'];prior=json.loads((H/'volume2-t512-R3.json').read_text());assert prior['points']==points
old=next(r for r in prior['records'] if r['T']==.1);fine=runs[-1];comparisons=[]
for kind,coarse in [('age',runs[0]),('angular_fixed_parent_mass',old)]:
 for i,p in enumerate(points):
  dp=abs(coarse['potential'][i]-fine['potential'][i])/max(abs(fine['potential'][i]),.01)
  da=float(np.linalg.norm(np.array(coarse['acceleration'][i])-fine['acceleration'][i])/max(np.linalg.norm(fine['acceleration'][i]),.01))
  comparisons.append(dict(kind=kind,point=p,potential_change=dp,force_change=da,potential_pass=dp<.02,force_pass=da<.05))
out=dict(runs=runs,previous=old,comparisons=comparisons)
(H/'adaptive-summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
for r in comparisons:print(r)
