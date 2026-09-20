"""SR-2 declared fixed-physics convergence campaign; first evidence immutable."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import time
from types import SimpleNamespace
import numpy as np
from evolution import rk4
from spatial_model import SpatialEvolution
from resolution_model import ResolutionEvolution

ROOT=Path(__file__).resolve().parent


def save(path,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(obj,indent=2,default=convert)+'\n',encoding='utf-8')


def controls():
    rng=np.random.default_rng(20260921);equivalence=[];gradients=[]
    for i in range(36):
        cfg=dict(s=i%2,eta=(0.,.08)[(i//2)%2],n=(12,24)[(i//4)%2],length=8.,radius=1.2,angle=.37*(i%3))
        compact=ResolutionEvolution(order=2,**cfg);original=SpatialEvolution(**cfg)
        y=compact.initial();f,pi,q,p=compact.unpack(y)
        envelope=np.exp(-np.sum(compact.x**2,axis=-1)/2)
        f[:]=rng.normal(0,.03,f.shape)*envelope;pi[:]=rng.normal(0,.02,pi.shape)*envelope
        q[:]+=rng.normal(0,.04,q.shape)
        a=compact.rhs(y);b=original.rhs(y)
        rhs_error=float(np.max(abs(a-b))/max(1.,np.max(abs(b))))
        ma=compact.metrics(y);mb=original.metrics(y)
        energy_error=abs(ma['total']-mb['total'])/max(1.,abs(mb['total']))
        cone_error=abs(ma['cone_residual']-mb['cone_residual'])
        equivalence.append(dict(index=i,config=cfg,rhs_error=rhs_error,energy_error=energy_error,cone_error=cone_error,passed=max(rhs_error,energy_error,cone_error)<1e-11))
    for order in (2,4):
        for i in range(24):
            model=ResolutionEvolution(order=order,s=i%2,n=12,length=8,radius=1.2)
            y=model.initial();f,pi,q,p=model.unpack(y)
            envelope=np.exp(-np.sum(model.x**2,axis=-1)/2)
            f[:]=rng.normal(0,.07,f.shape)*envelope;pi[:]=rng.normal(0,.04,pi.shape)*envelope
            q[:]+=rng.normal(0,.03,q.shape)
            fd,pd,qd,ppd=model.unpack(model.rhs(y))
            gradient=np.concatenate(((-model.dv*(pd+model.gamma*fd)).ravel(),(model.dv*fd).ravel(),(-ppd).ravel(),qd.ravel(),[0.]))
            direction=rng.normal(size=y.size);direction[-1]=0.;direction/=np.linalg.norm(direction);eps=2e-5
            numerical=(model.energy(y+eps*direction)-model.energy(y-eps*direction))/(2*eps)
            exact=np.dot(gradient,direction);error=abs(numerical-exact)/max(1.,abs(exact))
            gradients.append(dict(order=order,index=i,error=error,passed=error<2e-6))
    n=16;h=2*np.pi/n;theta=2*np.pi*np.arange(n)/n
    mesh=np.stack(np.meshgrid(theta,theta,theta,indexing='ij'))
    symbol=np.sum((np.sin(mesh)*(4-np.cos(mesh))/(3*h))**2,axis=0)+h**6/64*np.sum((16*np.sin(mesh/2)**4/h**4)**2,axis=0)
    mask=np.ones((n,n,n),bool);mask[0,0,0]=False
    fourier=dict(nonzero_modes=int(mask.sum()),minimum_nonzero_stiffness=float(symbol[mask].min()),zero_mode=float(symbol[0,0,0]),passed=bool(np.all(symbol[mask]>0) and symbol[0,0,0]==0))
    spectral=[]
    for mode in ((1,2,3),(7,4,1),(8,0,0),(8,8,8)):
        model=ResolutionEvolution(order=4,n=n,length=2*np.pi,eta=0)
        y=model.initial();f,pi,q,p=model.unpack(y)
        phase=sum(mode[j]*(model.x[...,j]+np.pi) for j in range(3));f[1]=.1*np.cos(phase)
        field=model.metrics(y)['field']
        expected=.5*model.dv*np.sum(f[1]**2)*(symbol[mode]+model.omega**2)
        error=abs(field/expected-1)
        spectral.append(dict(mode=mode,relative_error=error,passed=error<1e-12))
    derivative=[]
    for k in (1,3):
        errors=[]
        for n in (32,64,128):
            h=2*np.pi/n;x=np.arange(n)*h
            actual=ResolutionEvolution.d4(SimpleNamespace(dx=h),np.sin(k*x),0)
            errors.append(float(np.sqrt(np.mean((actual-k*np.cos(k*x))**2))))
        ratios=[errors[i]/errors[i+1] for i in range(2)]
        derivative.append(dict(mode=k,errors=errors,ratios=ratios,passed=min(ratios)>12))
    return dict(equivalence=equivalence,gradients=gradients,fourier=fourier,spectral=spectral,derivative=derivative,
                passed=all(x['passed'] for x in equivalence+gradients+spectral+derivative) and fourier['passed'])


def configurations():
    return [dict(name='o2-n28',order=2,n=28,length=14.,dt=.02),
            dict(name='o4-n28',order=4,n=28,length=14.,dt=.02),
            dict(name='o4-n40',order=4,n=40,length=15.,dt=.01),
            dict(name='o4-n56',order=4,n=56,length=14.,dt=.01),
            dict(name='o4-n70',order=4,n=70,length=14.,dt=.01),
            dict(name='o2-n56',order=2,n=56,length=14.,dt=.01),
            dict(name='o2-n70',order=2,n=70,length=14.,dt=.01),
            dict(name='o4-n56-time',order=4,n=56,length=14.,dt=.005),
            dict(name='o4-n56-rotated',order=4,n=56,length=14.,dt=.01,angle=np.pi/3)]


def main():
    out=ROOT/'resolution-v1';out.mkdir(exist_ok=False)
    names=('resolution-protocol.md','resolution_model.py','run_resolution.py','spatial_model.py','evolution.py','model.py')
    save(out/'manifest.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},numpy=np.__version__,configurations=configurations()))
    checks=controls();save(out/'controls.json',checks);print('controls',checks['passed'],flush=True)
    if not checks['passed']:raise RuntimeError('SR-2 controls failed; evidence preserved')
    rows=[]
    for cfg in configurations():
        print('start',cfg['name'],flush=True);wall=time.perf_counter();cpu=time.process_time()
        model=ResolutionEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')})
        y=model.initial();initial=y.copy();trace=[];steps=round(2.5/cfg['dt']);stride=max(1,round(.05/cfg['dt']))
        extent=0.
        for step in range(steps+1):
            q=model.unpack(y)[2];extent=max(extent,float(np.max(abs(q))+model.radius))
            if step%stride==0 or step==steps:
                metric=model.metrics(y);metric['time']=step*cfg['dt'];trace.append(metric)
            if step<steps:y=rk4(model,y,cfg['dt'])
        clearance=cfg['length']/2-1-extent-2.5*model.max_characteristic_seen
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
        cone=max(t['cone_residual'] for t in trace);edge=max(t['edge_amplitude'] for t in trace)
        bend=trace[-1]['probe_angle']-trace[0]['probe_angle']
        v0=np.array(trace[0]['probe_velocity']);vf=np.array(trace[-1]['probe_velocity'])
        unsigned=float(np.arctan2(np.linalg.norm(np.cross(v0,vf)),np.dot(v0,vf)))
        entry=dict(config=cfg,initial=trace[0],final=trace[-1],bend=bend,unsigned_bend=unsigned,ledger_drift=drift,cone_error=cone,edge_amplitude=edge,clearance=clearance,
                   maximum_source_extent=extent,maximum_characteristic=model.max_characteristic_seen,
                   momentum_drift=max(np.linalg.norm(np.array(t['momentum'])-trace[0]['momentum']) for t in trace),
                   angular_drift=max(np.linalg.norm(np.array(t['angular'])-trace[0]['angular']) for t in trace),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
        entry['passed']=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0
        np.savez_compressed(out/(cfg['name']+'.npz'),initial=initial,final=y)
        save(out/(cfg['name']+'.json'),dict(summary=entry,trace=trace));rows.append(entry)
        print('done',cfg['name'],'passed',entry['passed'],'bend',bend,'seconds',entry['wall_seconds'],flush=True)
    by={r['config']['name']:r for r in rows}
    b40=abs(by['o4-n40']['bend']);b56=abs(by['o4-n56']['bend']);b70=abs(by['o4-n70']['bend'])
    fine=abs(b70-b56)/b70;earlier=abs(b56-b40);latest=abs(b70-b56)
    method=abs(by['o2-n70']['bend']/by['o4-n70']['bend']-1)
    temporal=abs(by['o4-n56-time']['bend']/by['o4-n56']['bend']-1)
    rotmodel=ResolutionEvolution(n=56,length=14,angle=np.pi/3)
    vb=np.array(by['o4-n56-rotated']['final']['probe_velocity'])@rotmodel.rot
    vu=np.array(by['o4-n56']['final']['probe_velocity'])
    rotation_angle=float(np.arctan2(np.linalg.norm(np.cross(vb,vu)),np.dot(vb,vu)))
    rotation=rotation_angle/abs(by['o4-n56']['bend'])
    old=json.loads((ROOT/'spatial-v2/s1-fast.json').read_text())['summary']
    replay=max(abs(by['o2-n28']['bend']-old['bend']),abs(by['o2-n28']['final']['total']-old['final']['total']))
    accuracy=dict(finest_relative_bend_difference=fine,preceding_absolute_difference=earlier,latest_absolute_difference=latest,
                  finest_method_relative_difference=method,time_relative_bend_difference=temporal,rotation_direction_error_radians=rotation_angle,rotation_error_relative_to_bend=rotation,
                  replay_error=replay,passed=fine<.01 and latest<earlier and method<.01 and temporal<.001 and rotation<.01 and replay<1e-10)
    save(out/'summary.json',dict(runs=rows,accuracy=accuracy,numerical_passed=all(r['passed'] for r in rows),accuracy_passed=accuracy['passed']))
    print('complete',json.dumps(accuracy),flush=True)


if __name__=='__main__':main()
