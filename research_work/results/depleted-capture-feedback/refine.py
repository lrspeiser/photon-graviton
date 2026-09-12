"""Retain failed initial gates and refine affected cases without relaxing them."""
import json
import sys
import numpy as np
from run import HERE, Model

radius=float(sys.argv[1])
data=json.loads((HERE/'results.json').read_text())
output=[]
for pair,check in zip(data['results'],data['refinement_checks']):
    old=pair['refined']
    if old['R']!=radius or check['passes']:
        continue
    trials=[]
    for n in (512,1024):
        new=Model(radius,n,8).run(old['p'],old['C'],old['kappa0'],True)
        differences=(abs(np.array(old['diagnostics'])/new['diagnostics']-1)).tolist()
        trial=dict(result=new,relative_changes=differences,passes=max(differences)<.01)
        trials.append(trial)
        print(radius,old['p'],old['C'],old['kappa0'],n,differences,flush=True)
        if trial['passes']:
            break
        old=new
    output.append(dict(p=old['p'],C=old['C'],kappa0=old['kappa0'],R=radius,trials=trials))
    (HERE/f'refinement-{radius:g}.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8',newline='\n')
