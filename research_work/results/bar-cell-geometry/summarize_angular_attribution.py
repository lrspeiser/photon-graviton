"""Rank common angular regions without deleting or renormalizing their mass."""
from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent
a=json.loads((H/'angular-attribution-n1.json').read_text())['records'];b=json.loads((H/'angular-attribution-n2.json').read_text())['records']
d=np.array([y['acceleration'] for y in b])-np.array([x['acceleration'] for x in a]);net=d.sum(axis=0)
order=np.argsort(-np.linalg.norm(d,axis=1));q=np.linalg.norm(d,axis=1);q/=q.sum();proj=d@net/(net@net)
assert abs(q.sum()-1)<1e-12 and abs(proj.sum()-1)<1e-12
rows=[dict(parent_face=int(i),norm_fraction=float(q[i]),signed_projection=float(proj[i]),mass_first=a[i]['mass'],mass_second=b[i]['mass'],force_difference=d[i].tolist()) for i in order]
counts={str(f):int(np.searchsorted(np.cumsum(q[order]),f)+1) for f in (.5,.9,.99)}
out=dict(net_difference=net.tolist(),counts_for_norm_fraction=counts,ranked_regions=rows,
         top90_source_mass=float(sum(r['mass_second'] for r in rows[:counts['0.9']])),total_bin_mass=float(sum(r['mass_second'] for r in rows)))
(H/'angular-attribution-summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
print(counts,out['top90_source_mass'],out['total_bin_mass'])
