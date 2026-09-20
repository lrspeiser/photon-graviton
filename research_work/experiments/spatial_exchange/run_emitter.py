"""SE-2 preregistered controls and ten mixed-emitter spatial runs."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from emitter_model import ExchangeEvolution,rk4,COMMON
from model import ExchangeEvolution as OriginalEvolution

ROOT=Path(__file__).resolve().parent


def save(path,obj):
    def cast(value):
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)
    path.write_text(json.dumps(obj,indent=2,default=cast)+'\n',encoding='utf-8')


def controls():
    rng=np.random.default_rng(20260924);rows=[]
    for index in range(48):
        cfg=dict(emitter_angle=(0.,np.pi/4,np.pi/2)[index%3],chi=(-50.,0.,50.,200.)[index%4],mix=(-.5,.5)[(index//4)%2],emission=(0.,.4)[(index//8)%2],n=12,length=10.,radius=1.2)
        model=ExchangeEvolution(**cfg);y=model.initial();f,pi,q,p,Q,P=model.unpack(y)
        envelope=np.exp(-np.sum(model.x**2,axis=-1)/2)
        f[:]=rng.normal(0,.03,f.shape)*envelope;pi[:]=rng.normal(0,.02,pi.shape)*envelope
        q[:]+=rng.normal(0,.03,q.shape);Q[:]+=rng.normal(0,.03,6);P[:]=rng.normal(0,.03,6)
        fd,pd,qd,ppd,Qd,Pd=model.unpack(model.rhs(y))
        gradient=np.concatenate(((-model.dv*(pd+model.gamma*fd)).ravel(),(model.dv*fd).ravel(),(-ppd).ravel(),qd.ravel(),-Pd,Qd,[0.]))
        direction=rng.normal(size=y.size);direction[-1]=0;direction/=np.linalg.norm(direction);eps=2e-5
        numerical=(model.energy(y+eps*direction)-model.energy(y-eps*direction))/(2*eps);exact=float(np.dot(gradient,direction))
        error=abs(numerical-exact)/max(1.,abs(exact));metric=model.metrics(y)
        # Frozen local coefficients: six copies of the same longitudinal block.
        local=rng.normal(0,.03,4);U=.08*local[0];a=np.exp(4*U);c=np.exp(2*U)
        beta=c*.04*local[1:]/np.sqrt(1+.08**2*np.dot(local[1:],local[1:]))
        n=rng.normal(size=3);n/=np.linalg.norm(n);bn=np.dot(beta,n)
        principal=np.array([[bn,-1.],[-a,bn]]);sym=np.diag([a,1.])
        eigerr=float(np.max(abs(np.sort(np.linalg.eigvals(principal))-(bn+np.array([-c,c])))))
        symerr=float(np.max(abs(sym@principal-(sym@principal).T)))
        k=cfg['mix']*np.tanh(U/.001);m2=.75**2*np.exp(2*cfg['chi']*U)
        massmatrix=m2*np.array([[k*k,-k],[-k,1.]])
        mineig=float(np.linalg.eigvalsh(massmatrix).min())
        equivalence=0.
        if cfg['emitter_angle']==0:
            old=OriginalEvolution(**{k:v for k,v in cfg.items() if k!='emitter_angle'})
            equivalence=max(float(np.max(abs(old.rhs(y)-model.rhs(y)))),abs(old.energy(y)-model.energy(y)))
        rows.append(dict(equivalence_error=equivalence,index=index,config=cfg,gradient_error=error,minimum_mass=metric['minimum_mass'],cone_error=metric['cone_error'],principal_eigenvalue_error=eigerr,symmetrizer_error=symerr,mixing_matrix_minimum_eigenvalue=mineig,
                         passed=equivalence<1e-12 and error<2e-6 and metric['minimum_mass']>0 and metric['cone_error']<1e-10 and eigerr<1e-10 and symerr<1e-10 and mineig>-1e-12))
    return dict(rows=rows,passed=all(r['passed'] for r in rows))


def configurations():
    common=dict(n=32,length=16.,radius=.9,dt=.02,mix=.5)
    variants=[dict(name=f'emitter-{label}-chi{chi}',emitter_angle=theta,chi=float(chi))
              for label,theta in [('mixed',np.pi/4),('Y',np.pi/2)] for chi in (0,50,200)]
    variants += [dict(name='no-excitation',emitter_angle=np.pi/2,chi=200.,excitation=0.),
                 dict(name='time-refinement',emitter_angle=np.pi/2,chi=200.,dt=.01),
                 dict(name='space-refinement',emitter_angle=np.pi/2,chi=200.,n=40,dt=.01),
                 dict(name='rotation',emitter_angle=np.pi/2,chi=200.,angle=np.pi/3)]
    return [dict(common,**v) for v in variants]


def main():
    out=ROOT/'emitter-v1';out.mkdir(exist_ok=False)
    manifest=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('emitter-protocol.md','emitter_model.py','run_emitter.py','model.py')},helper_source=dict(path='research_work/experiments/common_cone/model.py',sha256=hashlib.sha256((COMMON/'model.py').read_bytes()).hexdigest()),configurations=configurations(),numpy=np.__version__)
    save(out/'manifest.json',manifest);check=controls();save(out/'controls.json',check)
    print('controls',check['passed'],flush=True)
    if not check['passed']:raise RuntimeError('SE-1 controls failed; first evidence preserved')
    rows=[]
    for cfg in configurations():
        print('start',cfg['name'],flush=True);wall=time.perf_counter();cpu=time.process_time()
        model=ExchangeEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')});y=model.initial();initial=y.copy()
        trace=[];extent=0.;steps=round(4/cfg['dt']);stride=round(.1/cfg['dt'])
        for step in range(steps+1):
            extent=max(extent,float(np.max(abs(model.unpack(y)[2]))+model.radius))
            if step%stride==0 or step==steps:
                metric=model.metrics(y);metric['time']=step*cfg['dt'];trace.append(metric)
            if step<steps:y=rk4(model,y,cfg['dt'])
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
        cone=max(t['cone_error'] for t in trace);edge=max(t['edge_amplitude'] for t in trace)
        clearance=model.length/2-1-extent-4*model.maximum_characteristic
        null=True
        if cfg['name'] in ('no-emission','no-excitation'):null=max(max(t['radiation_amplitude'],t['companion_amplitude']) for t in trace)<1e-12
        if cfg['name']=='emission-only':null=max(t['companion_amplitude'] for t in trace)<1e-12
        row=dict(config=cfg,initial=trace[0],final=trace[-1],relative_energy_drift=drift,cone_error=cone,edge_amplitude=edge,clearance=clearance,null_control_passed=null,
                 maximum_source_extent=extent,maximum_characteristic=model.maximum_characteristic,
                 momentum_drift=max(np.linalg.norm(np.array(t['momentum'])-trace[0]['momentum']) for t in trace),angular_drift=max(np.linalg.norm(np.array(t['angular'])-trace[0]['angular']) for t in trace),
                 wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,passed=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0 and null)
        np.savez_compressed(out/(cfg['name']+'.npz'),initial=initial,final=y)
        save(out/(cfg['name']+'.json'),dict(summary=row,trace=trace));rows.append(row)
        print('done',cfg['name'],'passed',row['passed'],'X',row['final']['radiation_propagation'],'Y',row['final']['companion_propagation'],flush=True)
    by={r['config']['name']:r for r in rows};comparisons=[];base=by['emitter-Y-chi200']
    base_model=ExchangeEvolution(chi=200.,emitter_angle=np.pi/2);base_y=np.load(out/'emitter-Y-chi200.npz')['final'];base_q=base_model.unpack(base_y)[2]
    for name in ('time-refinement','space-refinement','rotation'):
        row=by[name];cfg=row['config'];model=ExchangeEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')})
        y=np.load(out/(name+'.npz'))['final'];f,pi,q,p,Q,P=model.unpack(y)
        qback=q@model.rot;difference=float(np.max(np.linalg.norm(qback-base_q,axis=1)))
        energy=max(abs(row['final'][key]-base['final'][key])/max(abs(base['final'][key]),1e-12) for key in ('radiation_propagation','companion_propagation','mixing_potential'))
        if name=='sign-mirror':
            bf,bpi,bq,bp,bQ,bP=base_model.unpack(base_y);matter=max(np.max(abs(q-bq)),np.max(abs(p-bp)),np.max(abs(Q-bQ)),np.max(abs(P-bP)))
            field_error=abs(row['final']['field']/base['final']['field']-1)
            passed=matter<1e-7 and field_error<1e-8
        else:
            matter=None;field_error=None;passed=difference<(1e-3 if name=='time-refinement' else .01) and energy<(.01 if name=='time-refinement' else .05)
        comparisons.append(dict(name=name,position_difference=difference,channel_energy_relative_difference=energy,matter_state_difference=matter,total_field_relative_difference=field_error,passed=passed))
    result=dict(runs=rows,comparisons=comparisons,passed=all(r['passed'] for r in rows+comparisons))
    save(out/'summary.json',result);print('complete',result['passed'],flush=True)


if __name__=='__main__':main()
