from pathlib import Path
import sys,json
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole
from axis import evaluate
bar=FastMultipole.load(ROOT/'research_work/data-cache/bar-field/bar-L40.npz')
rows=[]
for z in (-3.,-1.,-.001,.001,1.,3.):
    p,a=evaluate(bar,[0,0,z])
    for angle in (1e-4,1e-6):
        pp,aa=bar.evaluate([abs(z)*angle,0,z])
        dp=float(abs(pp[0]-p[0])/max(abs(p[0]),1.))
        da=float(np.linalg.norm(aa[0]-a[0])/max(np.linalg.norm(a[0]),1.))
        rows.append(dict(z=z,angle=angle,potential_change=dp,force_change=da))
        if angle==1e-6:assert dp<1e-9 and da<1e-5
(HERE/'axis-checks.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
print('Axis checks pass; fine-angle maximum force change',max(r['force_change'] for r in rows if r['angle']==1e-6))
