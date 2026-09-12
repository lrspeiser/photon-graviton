"""Test known six-node quadratic positional interpolation with actual midpoint trajectories."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SYM=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])
TIMES=np.array([0,.005,.01,.025,.05,.1,.25])
mesh=json.loads((HERE/'mesh2.json').read_text(encoding='utf8'))
vertices=np.array(mesh['vertices']);faces=np.array(mesh['faces'])
inverses=np.linalg.inv(np.transpose(vertices[faces],(0,2,1)))
lookup=np.array(mesh['lookup']);signs=np.array(mesh['signs'])
parent=np.array(mesh['parent_faces'])
parent_inverse=np.linalg.inv(np.transpose(vertices[parent],(0,2,1)))
quadratic_nodes=[]
for i,(a,b,c) in enumerate(parent):
    ab=faces[4*i,1];ac=faces[4*i,2];bc=faces[4*i+1,2]
    for node,j,k in ((ab,a,b),(bc,b,c),(ac,a,c)):
        expected=vertices[j]+vertices[k];expected/=np.linalg.norm(expected)
        assert np.linalg.norm(expected-vertices[node])<1e-12
    quadratic_nodes.append([a,b,c,ab,bc,ac])
quadratic_nodes=np.array(quadratic_nodes)
def shape(b):
    a,c,d=b
    return np.array([a*(2*a-1),c*(2*c-1),d*(2*d-1),4*a*c,4*c*d,4*a*d])
nodal=np.array([[1,0,0],[0,1,0],[0,0,1],[.5,.5,0],[0,.5,.5],[.5,0,.5]])
assert np.max(abs(np.array([shape(b) for b in nodal])-np.eye(6)))<1e-12
for b in np.vstack([nodal,[.2,.3,.5],[1/3,1/3,1/3]]):assert abs(shape(b).sum()-1)<1e-12

def predict(direction,C):
    alpha=np.einsum('ijk,k->ij',parent_inverse,direction)
    candidates=np.flatnonzero(np.all(alpha>=-1e-10,axis=1)&(alpha.sum(axis=1)>0));assert len(candidates)
    i=int(candidates[0]);b=alpha[i]/sum(alpha[i])
    return np.sum(C[:,quadratic_nodes[i]]*shape(b)[None,:,None],axis=1),i

results=[]
selected=(float(sys.argv[1]),) if len(sys.argv)>1 else (1.,3.)
assert all(R in (1.,3.) for R in selected)
for R in selected:
    d=json.loads((HERE/f'prepared2-R{R:g}.json').read_text(encoding='utf8'))
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
        prediction,_=predict(direction,C);predicted.append(prediction)
    error=np.linalg.norm(np.stack(predicted,axis=1)-F,axis=2)
    records=[]
    for j,T in enumerate(TIMES):
        records.append(dict(T=float(T),rms_kpc=float(np.sqrt(sum(weights*error[j]**2))),
                            rms_over_R=float(np.sqrt(sum(weights*error[j]**2))/R),maximum_kpc=float(max(error[j])),
                            weight_above_002R=float(sum(weights[error[j]>.02*R]))))
    probes=json.loads((HERE/f'unsampled-probes-R{R:g}.json').read_text(encoding='utf8'))
    probe_results=[]
    for probe in probes['records']:
        direction=np.array(probe['direction'])
        predicted,i=predict(direction,C)
        actual=np.array(probe['actual_positions'])
        probe_results.append(dict(original_face=probe['face'],parent_face=i,
                                  error_over_R=(np.linalg.norm(predicted-actual,axis=1)/R).tolist()))
    results.append(dict(R=R,records=records,sample_errors=error.tolist(),targeted_probe_results=probe_results))
(HERE/(f'quadratic-geometry-R{selected[0]:g}.json' if len(selected)==1 else 'quadratic-geometry.json')).write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
for r in results:print(r['R'],r['records'])
