"""Follow up the angular bias in the weak-source candidate; no retuning."""
import json
import numpy as np
from run import HERE, Model

base=json.loads((HERE/'results.json').read_text())
results=[]
model=Model(100,256,16)
for k in (.1,10.,1000.):
    old=next(x['refined'] for x in base['results'] if x['refined']['R']==100 and x['refined']['p']==6 and x['refined']['C']==1 and x['refined']['kappa0']==k)
    new=model.run(6,1.,k,True)
    changes=abs(np.array(old['diagnostics'])/new['diagnostics']-1).tolist()
    results.append(dict(kappa0=k,result=new,relative_changes=changes,passes=max(changes)<.01))
    print(k,changes,flush=True)
(HERE/'angular-evolution-check.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
assert all(r['passes'] for r in results)
