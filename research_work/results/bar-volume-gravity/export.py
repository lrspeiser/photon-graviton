from pathlib import Path
import json,csv
import numpy as np
HERE=Path(__file__).resolve().parent
records=[];points=None
for n in (6,8):
    for layers in (64,128,256,512):
        d=json.loads((HERE/f'volume-n{n}-t{layers}.json').read_text(encoding='utf8'))
        assert len(d['records'])==6;points=d['points'];records.extend(d['records'])
def compare(a,b,kind):
    rows=[]
    for i,point in enumerate(points):
        dp=abs(a['potential'][i]-b['potential'][i])/max(abs(b['potential'][i]),.01)
        da=float(np.linalg.norm(np.array(a['acceleration'][i])-b['acceleration'][i])/max(np.linalg.norm(b['acceleration'][i]),.01))
        rows.append(dict(kind=kind,source_a=a['source_n'],source_b=b['source_n'],layers_a=a['layers'],layers_b=b['layers'],
                         R=b['R'],T=b['T'],point=point,potential_change=dp,force_change=da,potential_pass=dp<.02,force_pass=da<.05))
    return rows
comparisons=[]
for b in records:
    if b['layers']>64:
        a=next(r for r in records if r['source_n']==b['source_n'] and r['layers']==b['layers']//2 and r['R']==b['R'] and r['T']==b['T'])
        comparisons+=compare(a,b,'age')
    if b['source_n']==8 and b['layers']==512:
        a=next(r for r in records if r['source_n']==6 and r['layers']==512 and r['R']==b['R'] and r['T']==b['T'])
        comparisons+=compare(a,b,'angular')
rows=[]
for r in records:
    for i,p in enumerate(points):
        rows.append(dict(source_n=r['source_n'],layers=r['layers'],R=r['R'],T=r['T'],x=p[0],y=p[1],z=p[2],
                         potential=r['potential'][i],ax=r['acceleration'][i][0],ay=r['acceleration'][i][1],az=r['acceleration'][i][2]))
with (HERE/'summary.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
out=dict(records=records,comparisons=comparisons,mass_error=max(abs(r['mass']-1) for r in records),
         cell_checks=json.loads((HERE/'checks.json').read_text()),cache_checks=json.loads((HERE/'cache-checks.json').read_text()))
(HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
for kind in ('age','angular'):
    subset=[r for r in comparisons if r['kind']==kind and (kind=='angular' or r['layers_b']==512)]
    print(kind,len(subset),'potential fails',sum(not r['potential_pass'] for r in subset),'force fails',sum(not r['force_pass'] for r in subset),
          'max force change',max(r['force_change'] for r in subset))
print('Mass error',out['mass_error'])
