from pathlib import Path
from itertools import permutations
import json
import numpy as np
from scipy.special import roots_legendre
from tetra import gravity
HERE=Path(__file__).resolve().parent
tets=[]
for order in permutations(range(3)):
    v=[np.full(3,-1.)]
    for i in order:
        a=v[-1].copy();a[i]=1.;v.append(a)
    tets.append(v)
ref=json.loads((HERE.parent/'stream-force-quadrature/results.json').read_text())['reference']
p,a=gravity(tets,np.ones(6)/6,[[.137,.219,0]])
ep=abs(p[0]-ref['potential'])/abs(ref['potential'])
ea=float(np.linalg.norm(a[0]-np.r_[ref['force'],0])/np.linalg.norm(ref['force']))
assert max(ep,ea)<1e-9
regular=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1.]])
pc,ac=gravity([regular],[1],[[0,0,0]])
assert np.linalg.norm(ac)<1e-12
point=np.array([3.,2.,4.]);pe,ae=gravity([regular],[1],[point])
checks=[]
for n in (16,32):
    z,w=roots_legendre(n);z=(z+1)/2;w=w/2
    u,v,t=np.meshgrid(z,z,z,indexing='ij');wu,wv,wt=np.meshgrid(w,w,w,indexing='ij')
    bary=np.stack([(1-u)*(1-v)*(1-t),u,(1-u)*v,(1-u)*(1-v)*t],axis=-1)
    x=bary@regular;weights=6*wu*wv*wt*(1-u)**2*(1-v)
    dx=x-point;r=np.linalg.norm(dx,axis=-1)
    pp=-np.sum(weights/r);aa=np.sum(weights[...,None]*dx/r[...,None]**3,axis=(0,1,2))
    dp=float(abs(pp-pe[0])/abs(pe[0]));da=float(np.linalg.norm(aa-ae[0])/np.linalg.norm(ae[0]))
    assert max(dp,da)<1e-9
    checks.append(dict(nodes=n,potential_error=dp,force_error=da,mass=float(weights.sum())))
thin=[]
for height in (.01,.0001,.000001):
    tet=regular.copy();tet[:,2]*=height
    pt,at=gravity([tet],[1],[point])
    xx=bary@tet;dx=xx-point;rr=np.linalg.norm(dx,axis=-1)
    pp=-np.sum(weights/rr);aa=np.sum(weights[...,None]*dx/rr[...,None]**3,axis=(0,1,2))
    dp=float(abs(pp-pt[0])/abs(pp));da=float(np.linalg.norm(aa-at[0])/np.linalg.norm(aa))
    thin.append(dict(height=height,potential_relative_error=dp,force_relative_error=da))
    assert max(dp,da)<1e-7
out=dict(thin_checks=thin,cube_potential_error=float(ep),cube_force_error=ea,center_force=ac.tolist(),external=checks,passes=True)
(HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
print(out)
