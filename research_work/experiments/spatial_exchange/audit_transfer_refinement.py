"""Independent endpoint/loop/outer-region audit; supports partial archives."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_local_rotor import reconstruct

ROOT=Path(__file__).resolve().parent


def loop(a,h,r,count,rot):
    angle=np.arange(count)*2*np.pi/count
    pos=np.stack((r*np.cos(angle),r*np.sin(angle),np.zeros(count)),axis=-1)@rot.T
    tangent=np.stack((-r*np.sin(angle),r*np.cos(angle),np.zeros(count)),axis=-1)@rot.T
    grid=(pos+4)/h; index=np.floor(grid).astype(int); fraction=grid-index
    values=np.zeros((count,3))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                shift=np.array([i,j,k]); node=index+shift
                weight=np.prod(np.where(shift,fraction,1-fraction),axis=-1)
                values+=weight[:,None]*a[node[:,0],node[:,1],node[:,2]]
    return float(np.sum(values*tangent)*2*np.pi/count)


def outside(state,h,axis):
    n=state.shape[1]; points=np.arange(n)
    x=np.stack(np.meshgrid(*(points*h-4 for _ in range(3)),indexing='ij'),axis=-1)
    def density(q,p):
        gradient=[(np.take(q,(points+1)%n,axis=j)-np.take(q,(points-1)%n,axis=j))/(2*h) for j in range(3)]
        momentum=-np.stack([np.einsum('...a,...a->...',p,g) for g in gradient],axis=-1)
        return (np.cross(x,momentum)+np.cross(q,p))@axis
    field=density(state[0],state[1]); total=field+density(state[2],state[3]); r2=np.sum(x*x,axis=-1)
    return np.array([np.sum(value[r2>r*r])*h**3 for r in (1.2,1.5,2.) for value in (field,total)])


def main():
    out=ROOT/'local-transfer-refinement-v1'; manifest=json.loads((out/'manifest.json').read_text())
    checks=[hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==v for f,v in manifest['hashes'].items()]
    rows={}; maximum=0.
    for name in ('n40','n48','n64','n48-rotation'):
        if not (out/(name+'.json')).exists():continue
        row=json.loads((out/(name+'.json')).read_text()); rows[name]=row
        data=np.load(out/(name+'.npz')); trace=data['trace']; h=8/row['n']; rot=np.array(row['rotation'])
        checks.extend((np.max(abs(rot.T@rot-np.eye(3)))<1e-12,abs(np.linalg.det(rot)-1)<1e-12))
        for key,t in (('initial',0),('final',-1)):
            state=data[key]; result=reconstruct(state,h,2.); shifted=state.copy();shifted[3]-=state[0]
            result[:3]=reconstruct(shifted,h,2.)[:3]
            field=state.copy();field[2:]=0; ja=reconstruct(field,h,2.)[3:6]
            error=float(max(np.max(abs(result-trace[t,1:10])),np.max(abs(ja-trace[t,11:14])),np.max(abs(outside(state,h,rot[:,2])-trace[t,14:20]))))
            maximum=max(maximum,error);checks.append(error<1e-12)
        checks.append(abs(float(ja@rot[:,2])-row['A_J'])<1e-12)
        checks.append(np.max(abs(trace[-1,14:20]/row['initial_J']-row['outer_angular_fractions']))<1e-12)
        e=np.max(abs(trace[:,1]-trace[0,1]))/trace[0,1];j=np.max(np.linalg.norm(trace[:,4:7]-trace[0,4:7],axis=1))/np.linalg.norm(trace[0,4:7])
        edge=np.max(trace[:,10])/trace[0,1]
        checks.extend((abs(e-row['energy_drift'])<1e-12,abs(j-row['angular_drift'])<1e-12,abs(edge-row['edge_fraction'])<1e-12))
        for entry in row['loops']:
            for label,count in (('coarse',128),('fine',256)):
                checks.append(abs(loop(data['final'][0],h,entry['radius'],count,rot)-entry[label])<1e-12)
        passed=bool(np.isfinite(trace).all() and np.isfinite(data['final']).all() and e<1e-5 and j<.01 and edge<1e-5 and max(abs(r['coarse']-r['fine']) for r in row['loops'])<1e-5)
        checks.append(passed==row['passed'])
    complete=(out/'summary.json').exists();campaign=None
    if complete:
        summary=json.loads((out/'summary.json').read_text());checks.append(summary['runs']==list(rows.values()))
        old=json.loads((ROOT/'local-transfer-v1'/'space.json').read_text());comparison=[]
        for label,c,m,f in [('circulation',old['loops'][1]['fine'],rows['n48']['loops'][1]['fine'],rows['n64']['loops'][1]['fine']),('angular',old['A_Jz'],rows['n48']['A_J'],rows['n64']['A_J'])]:
            previous=abs(c-m)/max(abs(m),1e-10);err=abs(m-f)/max(abs(f),1e-10)
            comparison.append(dict(name=label+'-space',relative_difference=err,previous_difference=previous,passed=err<.05 and err<previous))
        for label,b,r in [('circulation',rows['n48']['loops'][1]['fine'],rows['n48-rotation']['loops'][1]['fine']),('angular',rows['n48']['A_J'],rows['n48-rotation']['A_J'])]:
            err=abs(b-r)/max(abs(b),1e-10);comparison.append(dict(name=label+'-rotation',relative_difference=err,passed=err<.02))
        checks.append(comparison==summary['comparisons']);campaign=summary['passed'];checks.append(campaign==all(r['passed'] for r in list(rows.values())+comparison))
    result=dict(passed=bool(all(checks)),checks=len(checks),completed_runs=len(rows),declared_runs=4,complete=complete,campaign_passed=campaign,maximum_endpoint_error=maximum)
    (ROOT/'local-transfer-refinement-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()
