"""CC-2 declared nonlinear 3D campaign. Refuses to overwrite first evidence."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import time
import json
import platform
import subprocess
from datetime import datetime,timezone
import numpy as np
from evolution import Evolution,rk4
from run_audit import HERE,ROOT,save,digest


def controls():
    model=Evolution(n=12,length=6)
    y=model.initial();rng=np.random.default_rng(20260921)
    y[:2*model.nf]=rng.normal(0,.03,2*model.nf)
    rhs=model.rhs(y);fd,dp,qv,pv=model.unpack(rhs)
    expected=np.concatenate(((-(dp+model.gamma*fd)*model.dv).ravel(),(fd*model.dv).ravel(),-pv.ravel(),qv.ravel(),[0.]))
    rows=[]
    def add(name,error,tolerance):rows.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<=tolerance)))
    indices=list(rng.choice(2*model.nf,32,replace=False))+list(range(2*model.nf,len(y)-1))
    for j in indices:
        eps=1e-8 if j>=len(y)-4 else 1e-6
        yp=y.copy();ym=y.copy();yp[j]+=eps;ym[j]-=eps
        num=(model.energy(yp)-model.energy(ym))/(2*eps)
        add(f'energy_gradient_{j}',abs(num-expected[j])/max(1,abs(expected[j])),2e-5)
    eps=1e-6
    rate=(model.energy(y+eps*rhs)-model.energy(y-eps*rhs))/(2*eps)+rhs[-1]
    add('damped_energy_ledger_derivative',abs(rate),2e-7)
    wrong=model.rhs(y,omit_self=True)
    wrongrate=(model.energy(y+eps*wrong)-model.energy(y-eps*wrong))/(2*eps)+wrong[-1]
    add('missing_self_source_negative_control',0 if abs(wrongrate)>1e-8 else 1,0)
    w=model.ingredients(y)[14]
    add('source_normalization',np.max(abs(w.sum(axis=(1,2,3))-1)),1e-12)
    return rows,dict(correct_ledger_rate=rate,omitted_self_source_rate=wrongrate)


def integrate(spec):
    dt=spec.get('dt',.02)
    model=Evolution(**{k:v for k,v in spec.items() if k not in ('id','dt')})
    y=model.initial();steps=round(2.5/dt)
    states=[];momenta=[];metrics=[];times=[]
    start=time.monotonic();extent=0
    for step in range(steps+1):
        f,pi,q,p=model.unpack(y)
        states.append(q.copy());momenta.append(p.copy());metrics.append(model.metrics(y));times.append(step*dt)
        extent=max(extent,float(np.max(abs(q))))
        if step<steps:y=rk4(model,y,dt)
    energy0=metrics[0]['ledger']
    drift=max(abs(m['ledger']-energy0) for m in metrics)/max(1,abs(energy0))
    cone=max(m['cone_residual'] for m in metrics)
    speed=max(m['max_characteristic'] for m in metrics)
    gates=dict(energy=drift<=1e-4,averaged_cone=cone<=1e-10)
    record=dict(spec=spec,seconds=time.monotonic()-start,times=times,metrics=metrics,energy_drift=drift,
                max_cone_residual=cone,conservative_boundary_clearance=model.length/2-extent-model.radius-2.5*speed,
                momentum_residual=max(np.linalg.norm(m['momentum']-metrics[0]['momentum']) for m in metrics),
                angular_residual=max(np.linalg.norm(m['angular']-metrics[0]['angular']) for m in metrics),
                gates=gates,passed=all(gates.values()))
    return record,dict(q=np.array(states),p=np.array(momenta),final_state=y)


def main():
    out=HERE/'evolution-v1';out.mkdir(exist_ok=False)
    files=[HERE/name for name in ('evolution-protocol.md','evolution.py','run_evolution.py','model.py','run_audit.py')]
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
         utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),numpy=np.__version__,hashes={p.relative_to(ROOT).as_posix():digest(p) for p in files}))
    tests,negative=controls();save(out/'controls.json',dict(checks=tests,negative_control=negative))
    if not all(t['passed'] for t in tests):raise AssertionError('Hamiltonian controls failed; first evidence preserved')
    print(f'{len(tests)} full nonlinear Hamiltonian controls pass; {negative}',flush=True)
    specs=[]
    for source in ('rest','rotating','reverse'):
        for eta in (0.,.08):
            for probe in ('none','photon','massive'):
                specs.append(dict(id=f'{source}-eta{eta:g}-{probe}',source=source,eta=eta,probe=probe))
    for source in ('rest','rotating','reverse'):
        base=dict(source=source,eta=.08,probe='photon')
        specs.extend([dict(base,id=f'{source}-time',dt=.01),dict(base,id=f'{source}-space',n=40,dt=.01),
                      dict(base,id=f'{source}-domain',n=48,length=18.)])
    for n,dt in ((32,.02),(40,.01)):
        specs.append(dict(id=f'rotated-n{n}',source='rotating',eta=.08,probe='photon',angle=np.pi/3,n=n,dt=dt))
    for radius in (.6,1.2):specs.append(dict(id=f'radius-{radius}',source='rotating',eta=.08,probe='photon',radius=radius))
    rows=[]
    for spec in specs:
        row,arrays=integrate(spec);rows.append(row)
        np.savez_compressed(out/(spec['id']+'.npz'),**arrays)
        save(out/'runs.json',rows)
        print(f'{len(rows)}/{len(specs)} {spec["id"]}: energy={row["energy_drift"]:.3g}, cone={row["max_cone_residual"]:.3g}, seconds={row["seconds"]:.1f}',flush=True)
    comparisons=[]
    for source in ('rest','rotating','reverse'):
        baseline=next(r for r in rows if r['spec']['id']==f'{source}-eta0.08-photon')
        q=np.load(out/(baseline['spec']['id']+'.npz'))['q']
        for kind in ('time','space','domain'):
            fine=next(r for r in rows if r['spec']['id']==f'{source}-{kind}')
            qf=np.load(out/(fine['spec']['id']+'.npz'))['q']
            err=np.linalg.norm(q-qf[::1 if kind=='domain' else 2],axis=2).max()
            energy=abs(baseline['metrics'][-1]['field']-fine['metrics'][-1]['field'])/max(fine['metrics'][-1]['field'],1e-10)
            comparisons.append(dict(source=source,kind=kind,position=err,field_energy_relative=energy,
                                    passed=bool(err<=(1e-4 if kind=='domain' else .01) and energy<=(1e-3 if kind=='domain' else .05))))
    for n in (32,40):
        baseline='rotating-eta0.08-photon' if n==32 else 'rotating-space'
        q=np.load(out/(baseline+'.npz'))['q']
        qr=np.load(out/(f'rotated-n{n}.npz'))['q']@Evolution(angle=np.pi/3).rot
        err=np.linalg.norm(q-qr,axis=2).max()
        comparisons.append(dict(kind='rotation',n=n,position=err,passed=bool(err<=.01)))
    save(out/'comparisons.json',comparisons)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(f'CC-2 complete: {sum(r["passed"] for r in rows)}/{len(rows)} evolution gates, {sum(r["passed"] for r in comparisons)}/{len(comparisons)} comparisons pass',flush=True)


if __name__=='__main__':main()
