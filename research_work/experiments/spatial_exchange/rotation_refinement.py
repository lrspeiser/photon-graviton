"""SE-RR1: finer rotated evolution with both declared measurement methods."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from local_transfer import rhs,ledger
from refine_local_transfer import rotation,outer
from loop_interpolation import measure

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'rotation-refinement-v1';out.mkdir(exist_ok=False)
    dependencies=('rotation_refinement.py','local_transfer.py','local_rotor.py','run_local_rotor.py','refine_local_transfer.py','loop_interpolation.py')
    reference=ROOT/'local-transfer-refinement-v1'/'n64.npz';prior=ROOT/'loop-interpolation-v1'/'results.json'
    manifest=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in dependencies},
                  inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (reference,prior)})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    n=64;h=8/n;dt=.01;rot=rotation();axis=rot[:,2]
    x=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)),indexing='ij'),axis=-1)
    f=.03*np.maximum(1-np.sum(x*x,axis=-1)/1.2**2,0)**3
    state=np.zeros((4,n,n,n,3));state[2]=f[...,None]*rot[:,0];state[3]=f[...,None]*rot[:,1];initial=state.copy()
    trace=[[0.,*ledger(state,x,h,2.,1.),*outer(state,x,h,axis)]]
    for step in range(100):
        k1=rhs(state,h,2.,1.);k2=rhs(state+.5*dt*k1,h,2.,1.)
        k3=rhs(state+.5*dt*k2,h,2.,1.);k4=rhs(state+dt*k3,h,2.,1.)
        state+=dt*(k1+2*k2+2*k3+k4)/6
        trace.append([(step+1)*dt,*ledger(state,x,h,2.,1.),*outer(state,x,h,axis)])
    trace=np.array(trace);np.savez_compressed(out/'rotated.npz',initial=initial,final=state,trace=trace)
    original=np.load(reference);previous=json.loads(prior.read_text());measurements=[];comparisons=[]
    for order in (1,3):
        for radius in (1.,1.5,2.):
            values=[measure(state[0],h,radius,count,rot,order) for count in (256,512)]
            unrotated=[measure(original['final'][0],h,radius,count,np.eye(3),order) for count in (256,512)]
            measurements.append(dict(order=order,radius=radius,rotated=values,unrotated=unrotated))
            if radius==1.5:
                error=abs(values[1]-unrotated[1])/max(abs(unrotated[1]),1e-10)
                before=next(r['rotation_difference'] for r in previous['comparisons'] if r['order']==order)
                comparisons.append(dict(order=order,relative_difference=error,previous_difference=before,passed=error<.02 and error<before))
    e=float(np.max(abs(trace[:,1]-trace[0,1]))/trace[0,1]);j=float(np.max(np.linalg.norm(trace[:,4:7]-trace[0,4:7],axis=1))/np.linalg.norm(trace[0,4:7]))
    edge=float(np.max(trace[:,10])/trace[0,1]);ja=float(trace[-1,11:14]@axis);jref=float(original['trace'][-1,13])
    angular=abs(ja-jref)/max(abs(jref),1e-10)
    numerical=bool(np.isfinite(state).all() and np.isfinite(trace).all() and e<1e-5 and j<.01 and edge<1e-5)
    result=dict(n=n,dt=dt,rotation=rot.tolist(),energy_drift=e,angular_drift=j,edge_fraction=edge,
                A_J=ja,angular_rotation_difference=angular,measurements=measurements,comparisons=comparisons,
                numerical_passed=numerical,passed=numerical and angular<.02 and all(r['passed'] for r in comparisons))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
