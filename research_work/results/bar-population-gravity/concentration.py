"""Diagnose contribution concentration without removing any source direction."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
rows=[]
for n,name in [(6,'orbits.json'),(8,'orbits8.json')]:
    d=json.loads((HERE/name).read_text())
    for R in (1.,3.):
        group=[r for r in d['records'] if r['R']==R];norm=sum(r['raw_weight'] for r in group)
        for j,T in enumerate((.05,.1,.25)):
            for k,point in enumerate(d['points']):
                terms=[]
                for r in group:
                    a=r['result']['kernels'][j];a=a['extra_refined'] or a['refined']
                    terms.append(np.array(a['acceleration'][k])*r['raw_weight']/norm)
                terms=np.array(terms);magnitudes=np.linalg.norm(terms,axis=1)
                total=np.linalg.norm(terms.sum(axis=0));absolute_sum=magnitudes.sum()
                ranked=np.argsort(-magnitudes)[:3]
                rows.append(dict(n=n,R=R,T=T,point=point,cancellation_ratio=float(total/absolute_sum),
                    largest_absolute_share=float(magnitudes.max()/absolute_sum),
                    largest_relative_to_net=float(magnitudes.max()/total),
                    top_three=[dict(mu=group[i]['mu'],phi=group[i]['phi'],weighted_force=terms[i].tolist(),
                                    absolute_share=float(magnitudes[i]/absolute_sum)) for i in ranked]))
(HERE/'concentration.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps([r for r in rows if r['T']==.25 and r['n']==8],indent=2))
