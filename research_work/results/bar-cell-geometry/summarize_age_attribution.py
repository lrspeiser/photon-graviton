from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent
a=json.loads((H/'age-attribution-n1.json').read_text())['records'][0]
b=json.loads((H/'age-attribution-n2.json').read_text())['records'][0]
d=np.array([y['acceleration'] for y in b['bins']])-np.array([x['acceleration'] for x in a['bins']]);net=d.sum(axis=0)
projection=d@net/(net@net);normshare=np.linalg.norm(d,axis=1)/np.linalg.norm(d,axis=1).sum()
out=dict(net_force_difference=net.tolist(),bins=[dict(age_min=x['age_min'],age_max=x['age_max'],mass=x['mass'],delta=v.tolist(),signed_projection_fraction=float(p),absolute_difference_norm_fraction=float(q)) for x,v,p,q in zip(a['bins'],d,projection,normshare)],projection_sum=float(sum(projection)))
assert abs(sum(projection)-1)<1e-12
(H/'age-attribution-summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
