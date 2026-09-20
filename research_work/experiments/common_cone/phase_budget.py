"""PF-0: reduced nonlinear feedback energy and threshold diagnostic."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from run_audit import HERE,ROOT,save,digest


def components(y,ratio):
    b,p=y[:2]
    return np.array([.5*np.exp(-2*b*b)*p*p,.5*b*b,ratio*np.exp(-b*b/2)])


def rhs(time,y,ratio,damping):
    b,p,q=y
    coefficient=np.exp(-2*b*b)
    velocity=coefficient*p
    return np.array([velocity,2*b*coefficient*p*p-b+ratio*b*np.exp(-b*b/2)-damping*velocity,damping*velocity**2])


def main():
    out=HERE/'phase-budget-v1';out.mkdir(exist_ok=False)
    sources=[HERE/name for name in ('phase-budget-protocol.md','phase_budget.py','run_audit.py')]
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                 hashes={p.relative_to(ROOT).as_posix():digest(p) for p in sources}))
    rng=np.random.default_rng(20260922);controls=[]
    for index in range(300):
        y=np.array([rng.uniform(-2,2),rng.uniform(-1,1),0.]);ratio=10**rng.uniform(-1,1)
        derivative=rhs(0,y,ratio,0)
        expected=np.array([-derivative[1],derivative[0]])
        fd=[]
        for j in range(2):
            shift=np.eye(3)[j]*1e-5
            fd.append((components(y+shift,ratio).sum()-components(y-shift,ratio).sum())/2e-5)
        error=np.max(abs(np.array(fd)-expected))/max(1,np.max(abs(expected)))
        damped=rhs(0,y,ratio,.05)
        ledger=(components(y+1e-6*damped,ratio).sum()-components(y-1e-6*damped,ratio).sum())/2e-6+damped[2]
        controls.append(dict(index=index,state=y,ratio=ratio,gradient_error=error,ledger_rate=ledger,
                             passed=bool(error<=2e-5 and abs(ledger)<=2e-7)))
    save(out/'controls.json',controls)
    if not all(c['passed'] for c in controls):raise AssertionError('PF-0 first control failure preserved')
    branches=[];runs=[]
    for ratio in (.5,.99,1,1+1e-6,1.01,1.1,2,10):
        amplitude=brentq(lambda b:1-ratio*np.exp(-b*b/2),0,10,xtol=1e-14) if ratio>1 else 0.
        expected=np.sqrt(2*np.log(ratio)) if ratio>1 else 0.
        potential=lambda b:.5*b*b+ratio*np.exp(-b*b/2)
        delta=1e-3
        curvature=(potential(amplitude+delta)-2*potential(amplitude)+potential(amplitude-delta))/delta**2
        reference=2*np.log(ratio) if ratio>1 else 1-ratio
        branch=dict(ratio=ratio,amplitude=amplitude,amplitude_reference=expected,curvature=curvature,curvature_reference=reference,
                    minimum_potential=potential(amplitude),U=-amplitude**2/2,
                    matter_fraction=np.exp(-amplitude**2/2),field_potential_fraction=.5*amplitude**2/ratio,
                    releasable_fraction=(ratio-potential(amplitude))/ratio,
                    passed=bool(abs(amplitude-expected)<=1e-9 and abs(curvature-reference)<=2e-5))
        branches.append(branch)
        for damping in (0,.05):
            growth=2*(ratio-1)/(np.sqrt(damping*damping+4*(ratio-1))+damping) if ratio>1 else 0.
            horizon=min(10000,max(200,40/max(growth,1e-3)))
            initial=np.array([1e-8,0.,0.]);energy0=components(initial,ratio).sum()
            times=np.unique(np.concatenate(([0.],np.geomspace(1e-4,horizon,499))))
            result=solve_ivp(lambda time,y:rhs(time,y,ratio,damping),(0,horizon),initial,method='Radau',
                             rtol=1e-9,atol=1e-12,t_eval=times,max_step=10)
            energy=np.array([components(y,ratio) for y in result.y.T])
            ledger=energy.sum(axis=1)+result.y[2]
            drift=np.max(abs(ledger-energy0))/max(1,abs(energy0))
            maxamp=np.max(abs(result.y[0]))
            sid=f't{ratio:g}-gamma{damping:g}'
            # Near-critical ratio gets an exact key to avoid rounding into t=1.
            sid=f't{ratio:.9g}-gamma{damping:g}'
            np.savez_compressed(out/(sid+'.npz'),times=result.t,state=result.y.T,energy=energy)
            runs.append(dict(id=sid,ratio=ratio,damping=damping,horizon=horizon,growth_rate=growth,
                             initial_energy=energy0,integration_success=result.success,message=result.message,
                             energy_ledger_drift=drift,max_amplitude=maxamp,final_amplitude=result.y[0,-1],
                             half_branch_reached=bool(ratio>1 and maxamp>=amplitude/2),final_components=energy[-1],
                             final_outlet_energy=result.y[2,-1],passed=bool(result.success and drift<=1e-6)))
            save(out/'runs.json',runs)
            print(f'{len(runs)}/16 {sid}: ledger={drift:.3g}, half branch={runs[-1]["half_branch_reached"]}',flush=True)
    save(out/'branches.json',branches)
    zero=solve_ivp(lambda time,y:rhs(time,y,2,.05),(0,200),np.zeros(3),method='Radau',rtol=1e-9,atol=1e-12)
    summary=dict(control_passes=sum(c['passed'] for c in controls),branch_passes=sum(b['passed'] for b in branches),
                 evolution_passes=sum(r['passed'] for r in runs),maximum_ledger_drift=max(r['energy_ledger_drift'] for r in runs),
                 exact_zero_seed_maximum=np.max(abs(zero.y)),
                 weak_U_band=[1e-8,1e-5],threshold_interval=[np.exp(1e-8),np.exp(1e-5)],
                 relative_threshold_width=np.expm1(1e-5-1e-8),
                 scope='Reduced seeded local mode only; not 3D source generation or stable galactic swirl.')
    save(out/'summary.json',summary)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(encoding='utf8'))


if __name__=='__main__':main()
