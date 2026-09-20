"""Preregistered SR-1 controls and eight short coupled evolutions."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from evolution import Evolution, rk4
from spatial_model import SpatialEvolution, local

ROOT = Path(__file__).resolve().parent


def plain(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    raise TypeError(type(x).__name__)


def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,default=plain)+'\n',encoding='utf-8')


def controls():
    rng = np.random.default_rng(20260920)
    records = []
    for index in range(240):
        s = (index//60)%2; lam = .5*((index//30)%2); m = index%2
        f = rng.normal(0,.4,4); pi = rng.normal(0,.3,4)
        grad = rng.normal(0,.3,(3,4)); p = rng.normal(size=3)
        result = local(f,pi,grad,p,m,s,lam)
        errors = []
        for variable,derivative in zip((f,pi,grad,p),result[1:5]):
            for k in range(variable.size):
                old = variable.flat[k]; eps = 2e-6
                variable.flat[k] = old+eps; plus = local(f,pi,grad,p,m,s,lam)[0]
                variable.flat[k] = old-eps; minus = local(f,pi,grad,p,m,s,lam)[0]
                variable.flat[k] = old
                numerical = (plus-minus)/(2*eps); exact = derivative.flat[k]
                errors.append(abs(numerical-exact)/max(1.,abs(exact)))
        a,d,c,beta = result[-1]
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        bn = np.dot(beta,n)
        # First-order state (Pi, longitudinal gradient): y_t=-K y_x.
        K = np.array([[bn,-d],[-a,bn]])
        S = np.diag([a,d])
        eigerr = np.max(abs(np.sort(np.linalg.eigvals(K))-(bn+np.array([-c,c]))))
        symmetry = np.max(abs(S@K-(S@K).T))
        # Full local Pi/gradient energy block for one field component.
        block = np.block([[np.array([[a]]),-beta[None]],[-beta[:,None],d*np.eye(3)]])
        minimum = float(np.linalg.eigvalsh(block).min())
        convexity = float(np.linalg.eigvalsh(result[5]).min())
        photon = local(f,pi,grad,p,0,s,lam)[4]
        cone = abs(np.linalg.norm(photon-beta)/c-1)
        passed = max(errors)<2e-6 and eigerr<1e-10 and symmetry<1e-10 and minimum>0 and convexity>-1e-10 and cone<1e-12
        records.append(dict(index=index,s=s,lam=lam,mass=m,gradient_error=max(errors),eigenvalue_error=eigerr,symmetry_error=symmetry,minimum_energy_eigenvalue=minimum,momentum_hessian_minimum=convexity,cone_error=cone,passed=passed))
    grid=[]; negative=[]
    for s in (0,1):
        for lam in (0.,.5):
            model=SpatialEvolution(s=s,lam=lam,n=12,length=8,radius=1.2)
            y=model.initial(); f,pi,q,p=model.unpack(y)
            envelope=np.exp(-np.sum(model.x**2,axis=-1)/2)
            f[:]=rng.normal(0,.15,f.shape)*envelope
            pi[:]=rng.normal(0,.12,pi.shape)*envelope
            q[:]+=rng.normal(0,.03,q.shape)
            def gradient(omit=False):
                dy=model.rhs(y,omit_self=omit); fd,pd,qd,ppd=model.unpack(dy)
                return np.concatenate(((-model.dv*(pd+model.gamma*fd)).ravel(),(model.dv*fd).ravel(),(-ppd).ravel(),qd.ravel(),[0.]))
            exact=gradient()
            for k in range(6):
                direction=rng.normal(size=y.size);direction[-1]=0;direction/=np.linalg.norm(direction)
                eps=2e-5
                numeric=(model.energy(y+eps*direction)-model.energy(y-eps*direction))/(2*eps)
                expected=np.dot(exact,direction)
                error=abs(numeric-expected)/max(1.,abs(expected))
                grid.append(dict(s=s,lam=lam,index=k,error=error,passed=error<2e-6))
            wrong=gradient(True); direction=wrong-exact;direction/=np.linalg.norm(direction)
            eps=2e-5; numeric=(model.energy(y+eps*direction)-model.energy(y-eps*direction))/(2*eps)
            good_error=abs(numeric-np.dot(exact,direction));bad_error=abs(numeric-np.dot(wrong,direction))
            negative.append(dict(s=s,lam=lam,correct_error=good_error,omitted_error=bad_error,passed=bad_error>max(100*good_error,1e-6)))
            if s==0 and lam==0:
                old=Evolution(n=12,length=8,radius=1.2)
                equivalence=dict(energy=abs(model.energy(y)-old.energy(y)),rhs=float(np.max(abs(model.rhs(y)-old.rhs(y)))))
                equivalence['passed']=max(equivalence.values())<1e-11
    return dict(local=records,grid=grid,negative=negative,cc2_equivalence=equivalence,passed=all(x['passed'] for x in records+grid+negative) and equivalence['passed'])


def main():
    out=ROOT/'spatial-v1';out.mkdir(exist_ok=False)
    sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('spatial-protocol.md','spatial_model.py','run_spatial.py','evolution.py','model.py')}
    save(out/'manifest.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources=sources,numpy=np.__version__))
    checks=controls();save(out/'controls.json',checks)
    print('controls',checks['passed'],flush=True)
    if not checks['passed']:raise RuntimeError('Controls failed; preserved first evidence')
    slow=200/299792.458;slowp=slow/np.sqrt(1-slow**2)
    configs=[]
    for s in (0,1):
        for tag,eta,p in (('scalar',0.,.2),('fast',.08,.2),('slow',.08,slowp)):
            configs.append(dict(name=f's{s}-{tag}',s=s,eta=eta,source_momentum=p,n=24,dt=.02))
    for n,dt,tag in ((24,.01,'time'),(32,.01,'space')):
        configs.append(dict(name='s1-fast-'+tag,s=1,eta=.08,source_momentum=.2,n=n,dt=dt))
    summaries=[]
    for cfg in configs:
        print('start',cfg['name'],flush=True)
        model=SpatialEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')})
        y=model.initial();initial=y.copy();trace=[]
        steps=round(2.5/cfg['dt'])
        for step in range(steps+1):
            metrics=model.metrics(y);metrics['time']=step*cfg['dt'];trace.append(metrics)
            if step<steps:y=rk4(model,y,cfg['dt'])
        maximum_speed=max(x['max_characteristic'] for x in trace)
        clearance=model.length/2-1-max(x['source_extent'] for x in trace)-2.5*maximum_speed
        drift=max(abs(x['ledger']-trace[0]['ledger']) for x in trace)/abs(trace[0]['ledger'])
        cone=max(x['cone_residual'] for x in trace)
        edge=max(x['edge_amplitude'] for x in trace)
        entry=dict(config=cfg,initial=trace[0],final=trace[-1],ledger_drift=drift,cone_error=cone,edge_amplitude=edge,clearance=clearance,
                   bend=trace[-1]['probe_angle']-trace[0]['probe_angle'],momentum_drift=max(np.linalg.norm(x['momentum']-trace[0]['momentum']) for x in trace),angular_drift=max(np.linalg.norm(x['angular']-trace[0]['angular']) for x in trace))
        entry['passed']=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0
        np.savez_compressed(out/(cfg['name']+'.npz'),initial=initial,final=y)
        save(out/(cfg['name']+'.json'),dict(summary=entry,trace=trace))
        summaries.append(entry);print('done',cfg['name'],'passed',entry['passed'],'bend',entry['bend'],flush=True)
    comparisons=[]
    for tag in ('time','space'):
        base=np.load(out/'s1-fast.npz')['final'];fine=np.load(out/f's1-fast-{tag}.npz')['final']
        bm=SpatialEvolution(s=1,n=24);fm=SpatialEvolution(s=1,n=24 if tag=='time' else 32)
        position=float(np.linalg.norm(bm.unpack(base)[2][-1]-fm.unpack(fine)[2][-1]))
        be=next(x['final']['field'] for x in summaries if x['config']['name']=='s1-fast')
        fe=next(x['final']['field'] for x in summaries if x['config']['name']==f's1-fast-{tag}')
        relative=abs(be-fe)/abs(fe)
        comparisons.append(dict(name=tag,probe_position_difference=position,field_energy_relative_difference=relative,passed=position<.01 and relative<.1))
    save(out/'summary.json',dict(runs=summaries,comparisons=comparisons,passed=all(x['passed'] for x in summaries+comparisons)))
    print('complete',all(x['passed'] for x in summaries+comparisons),flush=True)


if __name__=='__main__':main()
