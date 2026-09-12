from pathlib import Path
import json,hashlib
import numpy as np
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
rows=[]
for n in (6,8):
    d=json.loads((HERE/f'prepared{n}.json').read_text(encoding='utf8'));assert len(d['records'])==n*n
    for i,r in enumerate(d['records']):
        p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
        c=np.load(p);t=c['ages'];indices=np.unique(np.r_[np.arange(0,len(t),2),len(t)-1])
        a=CubicHermiteSpline(t,c['position'],c['velocity'])
        b=CubicHermiteSpline(t[indices],c['position'][indices],c['velocity'][indices])
        mid=(t[:-1]+t[1:])/2
        dx=float(np.max(np.linalg.norm(a(mid)-b(mid),axis=1)))
        dv=float(np.max(np.linalg.norm(a(mid,1)-b(mid,1),axis=1)))
        rows.append(dict(n=n,probe=i+1,position_change=dx,velocity_change=dv,passes=dx<1e-6 and dv<.01))
(HERE/'cache-checks.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
print('Pass',all(r['passes'] for r in rows),'max position',max(r['position_change'] for r in rows),'max velocity',max(r['velocity_change'] for r in rows))
