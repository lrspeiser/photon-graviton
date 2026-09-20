"""RC-1 coherent exchange with reciprocal environmental dynamics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
OMEGA=.2


def quantities(y,kappa0,A,delta0):
    Q,P,ar,ai,br,bi=y
    na=ar*ar+ai*ai;nb=br*br+bi*bi;cross=ar*br+ai*bi
    tangent=np.tanh(Q);derivative=1-tangent*tangent
    delta=delta0*tangent;kappa=kappa0*(1+A*tangent)
    source=P*P/2+OMEGA**2*Q*Q/2
    photon=(1+delta)*na;receiving=(1-delta)*nb;interaction=2*kappa*cross
    force=-OMEGA**2*Q-delta0*derivative*(na-nb)-2*kappa0*A*derivative*cross
    return dict(total=source+photon+receiving+interaction,norm=na+nb,source=source,photon=photon,receiving=receiving,
                interaction=interaction,force=force,delta=delta,kappa=kappa,converted_population=nb,cross=cross)


def rhs(t,y,kappa0,A,delta0,omit_reaction=False):
    Q,P,ar,ai,br,bi=y;terms=quantities(y,kappa0,A,delta0)
    delta=terms['delta'];kappa=terms['kappa']
    return np.array([P,-OMEGA**2*Q if omit_reaction else terms['force'],
                     (1+delta)*ai+kappa*bi,-(1+delta)*ar-kappa*br,
                     kappa*ai+(1-delta)*bi,-kappa*ar-(1-delta)*br])


def solve(cfg,reverse=None,omit=False):
    args=(cfg['kappa0'],cfg['A'],cfg['delta0'],omit)
    if reverse is None:
        initial=np.array([-1.,cfg['P0'],1.,0.,0.,0.]);interval=(0,100);times=np.linspace(0,100,1001)
    else:initial=reverse;interval=(100,0);times=[0]
    sol=solve_ivp(rhs,interval,initial,args=args,method='DOP853',rtol=1e-10,atol=1e-12,t_eval=times)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n',encoding='utf-8')


def main():
    out=ROOT/'reversible-conversion-v1';out.mkdir(exist_ok=False)
    manifest=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('reversible-conversion-protocol.md','reversible_conversion.py')})
    save(out/'manifest.json',manifest)
    gradients=[];rng=np.random.default_rng(20260922)
    for i in range(100):
        cfg=dict(kappa0=(.05,.2,.4)[i%3],A=(-.5,.5)[i%2],delta0=(0.,.3)[(i//2)%2])
        y=rng.normal(size=6);eps=1e-5
        plus=y.copy();minus=y.copy();plus[0]+=eps;minus[0]-=eps
        numerical=(quantities(plus,**cfg)['total']-quantities(minus,**cfg)['total'])/(2*eps)
        error=abs(numerical+quantities(y,**cfg)['force'])
        bound=1-np.sqrt(cfg['delta0']**2+(cfg['kappa0']*(1+abs(cfg['A'])))**2)
        gradients.append(dict(index=i,error=error,global_eigenvalue_lower_bound=bound,passed=error<2e-7 and bound>0))
    save(out/'gradient-controls.json',gradients)
    if not all(r['passed'] for r in gradients):raise RuntimeError('Derivative controls failed')
    arrays={};controls=[];matched={}
    for delta0,P0 in itertools.product((0.,.3),(0.,.2)):
        cfg=dict(kappa0=0.,A=0.,delta0=delta0,P0=P0);sol=solve(cfg);key=f'no-conversion-{len(controls)}'
        arrays[key]=sol.y;matched[(delta0,P0)]=sol.y
        population=sol.y[4]**2+sol.y[5]**2
        controls.append(dict(name=key,config=cfg,max_population=float(population.max()),passed=bool(population.max()<1e-12)))
    for kappa0 in (.05,.2,.4):
        cfg=dict(kappa0=kappa0,A=0.,delta0=0.,P0=0.);sol=solve(cfg);key=f'constant-{kappa0}'
        arrays[key]=sol.y;pop=sol.y[4]**2+sol.y[5]**2
        error=float(np.max(abs(pop-np.sin(kappa0*sol.t)**2)))
        controls.append(dict(name=key,config=cfg,population_error=error,passed=error<1e-8))
    rows=[];negatives=[]
    for kappa0,A,delta0,P0 in itertools.product((.05,.2,.4),(-.5,.5),(0.,.3),(0.,.2)):
        cfg=dict(kappa0=kappa0,A=A,delta0=delta0,P0=P0);key=f'run-{len(rows):02d}'
        sol=solve(cfg);arrays[key]=sol.y;terms=quantities(sol.y,kappa0,A,delta0)
        energy=float(np.max(abs(terms['total']-terms['total'][0]))/abs(terms['total'][0]))
        norm=float(np.max(abs(terms['norm']-1)))
        reverse=solve(cfg,reverse=sol.y[:,-1]);back=float(np.max(abs(reverse.y[:,-1]-sol.y[:,0])))
        qdifference=float(np.max(abs(sol.y[0]-matched[(delta0,P0)][0])))
        identity=delta0!=0 or (qdifference<1e-7 and np.max(abs(terms['cross']))<1e-8)
        rows.append(dict(name=key,config=cfg,relative_energy_drift=energy,norm_drift=norm,reversal_error=back,
                         peak_receiving_population=float(terms['converted_population'].max()),final_receiving_population=float(terms['converted_population'][-1]),
                         initial_energy=float(terms['total'][0]),final_energy=float(terms['total'][-1]),
                         source_energy_change=float(terms['source'][-1]-terms['source'][0]),final_interaction_energy=float(terms['interaction'][-1]),
                         max_environment_difference_from_no_conversion=qdifference,symmetric_no_force_identity=bool(identity),
                         passed=energy<1e-8 and norm<1e-8 and back<1e-7 and identity))
        if delta0!=0:
            wrong=solve(cfg,omit=True);arrays[key+'-omitted-reaction']=wrong.y
            wt=quantities(wrong.y,kappa0,A,delta0)['total'];error=float(np.max(abs(wt-wt[0]))/abs(wt[0]))
            negatives.append(dict(name=key+'-omitted-reaction',config=cfg,relative_energy_drift=error,detected=error>1e-5))
    np.savez_compressed(out/'trajectories.npz',time=np.linspace(0,100,1001),**arrays)
    G=6.67430e-11;mass=2e41;required_mass=(2e5)**2*(6e20-6e18)/G
    budgets=[]
    for s,epsilon,fraction in itertools.product((0,1),(1.,1e-3,1e-6),(0.,.1,.5,1.)):
        chi=1+s;required=required_mass/(epsilon*mass)
        required_g=None if fraction==0 else max(0.,(required-(1-fraction)*chi)/fraction)
        effective=(1-fraction)*chi+fraction*chi
        budgets.append(dict(s=s,available_fraction=epsilon,converted_energy_fraction=fraction,equal_response_active_coefficient=effective,
                            required_effective_coefficient=required,minimum_receiving_coefficient_if_other_fixed=required_g,
                            equal_response_budget_sufficient=effective>=required))
    result=dict(runs=rows,controls=controls,negative_controls=negatives,budgets=budgets,required_active_mass_kg=required_mass,
                passed=all(r['passed'] for r in rows+controls+gradients) and all(r['detected'] for r in negatives))
    save(out/'summary.json',result)
    print(json.dumps(dict(passed=result['passed'],runs=len(rows),max_energy_drift=max(r['relative_energy_drift'] for r in rows),max_norm_drift=max(r['norm_drift'] for r in rows),max_reversal_error=max(r['reversal_error'] for r in rows),peak_conversion=max(r['peak_receiving_population'] for r in rows))))
    if not result['passed']:raise RuntimeError('RC-1 acceptance failure; first archive preserved')


if __name__=='__main__':main()
