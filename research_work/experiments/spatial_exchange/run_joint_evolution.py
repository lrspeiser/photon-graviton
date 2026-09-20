"""SE-JE1: moving matter and a backreacting photon start with empty fields."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from joint_transfer import JointTransfer
from local_rotor import curl,derivative
from point_probes import evaluate as probe
from refine_local_transfer import rotation

ROOT=Path(__file__).resolve().parent


def initial(model,rot,speed):
    theta=np.arange(6)*np.pi/3
    positions=np.vstack((.7*np.column_stack((np.cos(theta),np.sin(theta),np.zeros(6))),[-1.,1.,0.]))@rot.T
    momenta=np.vstack((speed*np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(6))),[1e-4,0.,0.]))@rot.T
    return np.concatenate((np.zeros(2*model.nf),positions.ravel(),momenta.ravel()))


def metrics(model,y,rot,coords):
    value=model.evaluate(y,False);f,m,q,p=model.unpack(y);phi=f[0]
    A,Pi,Q,P=[np.moveaxis(x,0,-1) for x in (f[1:4],m[1:4],f[4:7],m[4:7])]
    C=np.exp(2*model.g*phi);a=C*C;beta=C[...,None]*model.kappa*model.cq*model.eta*A/np.sqrt(1+model.eta**2*np.sum(A*A,axis=-1))[...,None]
    K=P-model.direct*A-model.epsilon*np.cross(curl(A,model.h),Q)
    current=-np.stack([np.sum(m*derivative(f,j+1,model.h),axis=0) for j in range(3)],axis=-1)
    spin=np.cross(A,Pi)+np.cross(Q,P)
    momentum=np.sum(current,axis=(0,1,2))*model.dv+np.sum(p,axis=0)
    angular=np.sum(np.cross(coords,current)+spin,axis=(0,1,2))*model.dv+np.sum(np.cross(q,p),axis=0)
    def sector_j(field,mom):
        density=-np.stack([np.sum(mom*derivative(field,j,model.h),axis=-1) for j in range(3)],axis=-1)
        return np.sum(np.cross(coords,density)+np.cross(field,mom),axis=(0,1,2))*model.dv
    eq=.5*np.sum(a[...,None]*K*K+model.Omega**2*Q*Q,axis=-1)
    for j in range(3):
        plus=(np.roll(Q,-1,j)-Q)/model.h;minus=(Q-np.roll(Q,1,j))/model.h
        eq+=model.cq**2*.25*np.sum(plus*plus+minus*minus,axis=-1)-beta[...,j]*np.sum(K*derivative(Q,j,model.h),axis=-1)
    v=probe(q[-1],p[-1],0.,model.length,np.sqrt(C),C,np.moveaxis(beta,-1,0))['velocity']
    edge=np.max(abs(f[:,np.max(abs(coords),axis=-1)>=4]))
    return dict(energy=value['energy'],field=value['field'],matter=value['matter'],momentum=momentum.tolist(),angular=angular.tolist(),
                A_angular=sector_j(A,Pi).tolist(),Q_angular=sector_j(Q,P).tolist(),KQ_energy=float(np.sum(eq)*model.dv),
                cone_error=value['cone_error'],edge_amplitude=float(edge),photon_velocity=v.tolist(),bend=float(np.arctan2(v@rot[:,1],v@rot[:,0])))


def main():
    out=ROOT/'joint-evolution-v1';out.mkdir(exist_ok=False)
    files=('joint_transfer.py','run_joint_evolution.py','local_rotor.py','point_probes.py','refine_local_transfer.py')
    manifest=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    configs=[('free',24,.02,2.,1.,.2,False,True),('baseline',24,.02,0.,0.,.2,False,False),('combined',24,.02,2.,1.,.2,False,False),
             ('slow-baseline',24,.02,0.,0.,.0007,False,False),('slow-combined',24,.02,2.,1.,.0007,False,False),
             ('time',24,.01,2.,1.,.2,False,False),('space',32,.01,2.,1.,.2,False,False),('rotation',24,.02,2.,1.,.2,True,False)]
    rows=[]
    for name,n,dt,eps,direct,speed,turn,free in configs:
        print('start',name,flush=True);start=time.monotonic();rot=rotation() if turn else np.eye(3)
        model=JointTransfer(n=n,length=10.,masses=[1.]*6+[0.],epsilon=eps,direct=direct,g=0. if free else .08,eta=0. if free else .08)
        coords=np.stack(np.meshgrid(*(np.arange(n)*model.h-5 for _ in range(3)),indexing='ij'),axis=-1)
        y=initial(model,rot,speed);original=y.copy();trace=[dict(time=0.,**metrics(model,y,rot,coords))]
        f,m,q,p=model.unpack(y);particles=[np.stack((q.copy(),p.copy()))]
        for step in range(round(2/dt)):
            k1=model.evaluate(y);k2=model.evaluate(y+.5*dt*k1);k3=model.evaluate(y+.5*dt*k2);k4=model.evaluate(y+dt*k3)
            y+=dt*(k1+2*k2+2*k3+k4)/6
            trace.append(dict(time=(step+1)*dt,**metrics(model,y,rot,coords)));f,m,q,p=model.unpack(y);particles.append(np.stack((q.copy(),p.copy())))
        energy_error=max(abs(r['energy']-trace[0]['energy']) for r in trace)/trace[0]['energy']
        j0=np.array(trace[0]['angular']);angular_error=max(np.linalg.norm(np.array(r['angular'])-j0) for r in trace)/max(np.linalg.norm(j0),1e-8)
        pdiff=max(np.linalg.norm(np.array(r['momentum'])-trace[0]['momentum']) for r in trace)
        cone=max(r['cone_error'] for r in trace);edge=max(r['edge_amplitude'] for r in trace)
        checks=dict(finite=bool(np.isfinite(y).all()),energy=energy_error<1e-5,angular=bool(angular_error<.01),cone=cone<1e-10,edge=edge<1e-5)
        free_error=None
        if free:
            _,_,q0,p0=model.unpack(original);v0=p0/np.sqrt(model.masses**2+np.sum(p0*p0,axis=1))[:,None]
            free_error=float(max(np.max(abs(q-q0-2*v0)),np.max(abs(p-p0))))
            checks['free']=bool(np.max(abs(y[:2*model.nf]))<1e-12 and free_error<1e-10)
        row=dict(name=name,n=n,dt=dt,epsilon=eps,direct=direct,speed=speed,rotation=rot.tolist(),free=free,checks=checks,passed=all(checks.values()),
                 energy_drift=energy_error,angular_drift=float(angular_error),momentum_drift=float(pdiff),free_error=free_error,bend=trace[-1]['bend']-trace[0]['bend'],
                 trace=trace,seconds=time.monotonic()-start)
        np.savez_compressed(out/(name+'.npz'),initial=original,final=y,particles=np.array(particles))
        (out/(name+'.json')).write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
        print(json.dumps({k:row[k] for k in ('name','passed','energy_drift','angular_drift','bend','seconds')}),flush=True)
    by={r['name']:r for r in rows};base=by['combined'];comparisons=[]
    for name,tolerance in (('time',.001),('space',.05)):
        error=abs(base['bend']-by[name]['bend'])/max(abs(by[name]['bend']),1e-8)
        comparisons.append(dict(name=name,relative_error=error,passed=error<tolerance))
    v=np.array(by['rotation']['trace'][-1]['photon_velocity'])@rotation();u=np.array(base['trace'][-1]['photon_velocity'])
    error=float(np.arctan2(np.linalg.norm(np.cross(v,u)),v@u)/max(abs(base['bend']),1e-8))
    comparisons.append(dict(name='rotation',relative_error=error,passed=error<.02))
    summary=dict(runs=[{k:v for k,v in row.items() if k!='trace'} for row in rows],comparisons=comparisons,passed=all(r['passed'] for r in rows+comparisons))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(dict(passed=summary['passed'],comparisons=comparisons)),flush=True)


if __name__=='__main__':main()
