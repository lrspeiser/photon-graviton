from pathlib import Path
import json,hashlib
import numpy as np
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SYM=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])
TIMES=np.array([0,.005,.01,.025,.05,.1,.25])
def load(n,R):
    d=json.loads((HERE.parent/f'bar-volume-gravity/prepared{n}.json').read_text(encoding='utf8'))
    assert len(d['records'])==n*n and all(r['passes'] for r in d['records'])
    directions=[];positions=[];weights=[]
    for r in [r for r in d['records'] if r['R']==R]:
        p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
        c=np.load(p);s=CubicHermiteSpline(c['ages'],c['position'],c['velocity'])
        directions.extend(np.array(r['direction'])*SYM)
        positions.append(s(TIMES)[:,None,:]*SYM[None,:,:])
        weights.extend([r['raw_weight']/4]*4)
    return np.array(directions),np.concatenate(positions,axis=1),np.array(weights)/sum(weights)

results=[]
for R in (1.,3.):
    coarse,C,_=load(6,R);fine,F,weights=load(8,R)
    faces=np.sort(ConvexHull(coarse).simplices,axis=1)
    inverses=np.linalg.inv(np.transpose(coarse[faces],(0,2,1)))
    assigned=[];bary=[];ray_errors=[]
    for direction in fine:
        alpha=np.einsum('ijk,k->ij',inverses,direction)
        candidates=np.flatnonzero(np.all(alpha>=-1e-10,axis=1)&(alpha.sum(axis=1)>0))
        assert len(candidates)>0
        index=int(candidates[0]);a=alpha[index]/sum(alpha[index])
        assert min(a)>-1e-10 and abs(sum(a)-1)<1e-12
        ray=a@coarse[faces[index]]
        ray_errors.append(float(np.linalg.norm(ray/np.linalg.norm(ray)-direction)))
        assigned.append(index);bary.append(a)
    bary=np.array(bary);assigned=np.array(assigned)
    predicted=np.sum(C[:,faces[assigned],:]*bary[None,:,:,None],axis=2)
    error=np.linalg.norm(predicted-F,axis=2)
    records=[]
    for j,T in enumerate(TIMES):
        records.append(dict(T=float(T),rms_kpc=float(np.sqrt(sum(weights*error[j]**2))),
                            maximum_kpc=float(max(error[j])),rms_over_R=float(np.sqrt(sum(weights*error[j]**2))/R),
                            weight_above_002R=float(sum(weights[error[j]>.02*R]))))
    cells=[]
    for index in np.unique(assigned):
        mask=assigned==index;w=weights[mask]
        cells.append(dict(face=int(index),vertices=coarse[faces[index]].tolist(),tested_directions=int(sum(mask)),
                          tested_source_weight=float(sum(w)),final_rms_kpc=float(np.sqrt(sum(w*error[-1,mask]**2)/sum(w))),
                          peak_error_kpc=float(np.max(error[:,mask]))))
    results.append(dict(R=R,records=records,cells=sorted(cells,key=lambda r:-r['final_rms_kpc']),
                        maximum_ray_reconstruction_error=max(ray_errors),assigned_faces=assigned.tolist(),
                        sample_errors=error.tolist(),sample_weights=weights.tolist()))
(HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
for r in results:print(r['R'],r['records'])
