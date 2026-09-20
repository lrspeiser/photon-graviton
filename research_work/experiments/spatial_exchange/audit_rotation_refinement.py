"""Independent endpoint and Lagrange-product reconstruction for SE-RR1."""
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from audit_local_rotor import reconstruct
from audit_transfer_refinement import outside

ROOT=Path(__file__).resolve().parent


def integral(field,radius,count,rot,order):
    theta=np.arange(count)*2*np.pi/count
    pos=np.array([radius*np.cos(theta),radius*np.sin(theta),np.zeros(count)]).T@rot.T
    tangent=np.array([-radius*np.sin(theta),radius*np.cos(theta),np.zeros(count)]).T@rot.T
    h=8/field.shape[0];grid=(pos+4)/h;base=np.floor(grid).astype(int);t=grid-base
    nodes=(0,1) if order==1 else (-1,0,1,2)
    weights={node:np.prod([(t-other)/(node-other) for other in nodes if other!=node],axis=0) for node in nodes}
    values=np.zeros_like(pos)
    for i,j,k in itertools.product(nodes,repeat=3):
        index=(base+np.array([i,j,k]))%field.shape[0]
        weight=weights[i][:,0]*weights[j][:,1]*weights[k][:,2]
        values+=weight[:,None]*field[index[:,0],index[:,1],index[:,2]]
    return float(np.sum(values*tangent)*2*np.pi/count)


def main():
    out=ROOT/'rotation-refinement-v1';r=json.loads((out/'results.json').read_text());manifest=json.loads((out/'manifest.json').read_text())
    checks=[hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==v for group in ('hashes','inputs') for f,v in manifest[group].items()]
    data=np.load(out/'rotated.npz');reference=np.load(ROOT/'local-transfer-refinement-v1'/'n64.npz');trace=data['trace'];rot=np.array(r['rotation']);h=8/r['n'];maximum=0.
    for key,t in (('initial',0),('final',-1)):
        state=data[key];values=reconstruct(state,h,2.);shift=state.copy();shift[3]-=state[0];values[:3]=reconstruct(shift,h,2.)[:3]
        field=state.copy();field[2:]=0;ja=reconstruct(field,h,2.)[3:6]
        error=float(max(np.max(abs(values-trace[t,1:10])),np.max(abs(ja-trace[t,11:14])),np.max(abs(outside(state,h,rot[:,2])-trace[t,14:20]))))
        maximum=max(maximum,error);checks.append(error<1e-12)
    measurements=[]
    for item in r['measurements']:
        x=dict(order=item['order'],radius=item['radius'])
        for label,field,matrix in (('rotated',data['final'][0],rot),('unrotated',reference['final'][0],np.eye(3))):
            x[label]=[integral(field,item['radius'],count,matrix,item['order']) for count in (256,512)]
            checks.append(np.max(abs(np.array(x[label])-item[label]))<1e-12)
        measurements.append(x)
    old=json.loads((ROOT/'loop-interpolation-v1'/'results.json').read_text())
    comparison_pass=[]
    for row in r['comparisons']:
        entry=next(x for x in measurements if x['order']==row['order'] and x['radius']==1.5)
        error=abs(entry['rotated'][1]-entry['unrotated'][1])/max(abs(entry['unrotated'][1]),1e-10)
        before=next(x['rotation_difference'] for x in old['comparisons'] if x['order']==row['order']);passed=error<.02 and error<before
        checks.extend((abs(error-row['relative_difference'])<1e-12,before==row['previous_difference'],passed==row['passed']));comparison_pass.append(passed)
    e=float(np.max(abs(trace[:,1]-trace[0,1]))/trace[0,1]);j=float(np.max(np.linalg.norm(trace[:,4:7]-trace[0,4:7],axis=1))/np.linalg.norm(trace[0,4:7]));edge=float(np.max(trace[:,10])/trace[0,1])
    jref=reference['trace'][-1,13];turn=abs(float(ja@rot[:,2])-jref)/max(abs(jref),1e-10)
    numerical=bool(np.isfinite(data['final']).all() and np.isfinite(trace).all() and e<1e-5 and j<.01 and edge<1e-5)
    checks.extend((abs(e-r['energy_drift'])<1e-12,abs(j-r['angular_drift'])<1e-12,abs(edge-r['edge_fraction'])<1e-12,
                   abs(turn-r['angular_rotation_difference'])<1e-12,numerical==r['numerical_passed'],r['passed']==(numerical and turn<.02 and all(comparison_pass))))
    result=dict(passed=bool(all(checks)),checks=len(checks),maximum_endpoint_error=maximum,campaign_passed=r['passed'],scope='Independent endpoint and loop reconstruction; intermediate extrema from recorded traces.')
    (ROOT/'rotation-refinement-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()
