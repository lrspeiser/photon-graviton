"""Evaluate the new mesh using the same actual source8 reference trajectories."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SYM=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])
TIMES=np.array([0,.005,.01,.025,.05,.1,.25])
mesh=json.loads((HERE/'mesh.json').read_text(encoding='utf8'))
vertices=np.array(mesh['vertices']);faces=np.array(mesh['faces'])
inverses=np.linalg.inv(np.transpose(vertices[faces],(0,2,1)))
lookup=np.array(mesh['lookup']);signs=np.array(mesh['signs'])
results=[]
for R in (1.,3.):
    d=json.loads((HERE/f'prepared-R{R:g}.json').read_text(encoding='utf8'))
    assert len(d['records'])==d['expected']==len(mesh['directions'])
    assert all(r['passes'] for r in d['records'])
    paths=[]
    for r in d['records']:
        p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
        c=np.load(p);paths.append(CubicHermiteSpline(c['ages'],c['position'],c['velocity'])(TIMES))
    C=np.stack(paths,axis=1)[:,lookup]*signs[None,:,:]
    old=json.loads((HERE.parent/'bar-volume-gravity/prepared8.json').read_text(encoding='utf8'))
    directions=[];positions=[];weights=[]
    for r in [r for r in old['records'] if r['R']==R]:
        p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
        c=np.load(p);s=CubicHermiteSpline(c['ages'],c['position'],c['velocity'])
        directions.extend(np.array(r['direction'])*SYM);positions.append(s(TIMES)[:,None,:]*SYM)
        weights.extend([r['raw_weight']/4]*4)
    weights=np.array(weights)/sum(weights);F=np.concatenate(positions,axis=1)
    predicted=[]
    for direction in directions:
        alpha=np.einsum('ijk,k->ij',inverses,direction)
        candidates=np.flatnonzero(np.all(alpha>=-1e-10,axis=1)&(alpha.sum(axis=1)>0));assert len(candidates)>0
        i=int(candidates[0]);b=alpha[i]/sum(alpha[i])
        predicted.append(np.sum(C[:,faces[i]]*b[None,:,None],axis=1))
    error=np.linalg.norm(np.stack(predicted,axis=1)-F,axis=2)
    records=[]
    for j,T in enumerate(TIMES):
        records.append(dict(T=float(T),rms_kpc=float(np.sqrt(sum(weights*error[j]**2))),
                            rms_over_R=float(np.sqrt(sum(weights*error[j]**2))/R),maximum_kpc=float(max(error[j])),
                            weight_above_002R=float(sum(weights[error[j]>.02*R]))))
    results.append(dict(R=R,records=records,sample_errors=error.tolist()))
(HERE/'nested-geometry.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
for r in results:print(r['R'],r['records'])
