"""RC-1 unchanged equations and gates, stricter integration tolerances."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from reversible_conversion import ROOT,rhs,quantities,save


def solve(cfg,reverse=None,omit=False):
    initial=np.array([-1.,cfg['P0'],1.,0.,0.,0.]) if reverse is None else reverse
    sol=solve_ivp(rhs,(0,100) if reverse is None else (100,0),initial,
                  args=(cfg['kappa0'],cfg['A'],cfg['delta0'],omit),method='DOP853',rtol=1e-12,atol=1e-14,
                  t_eval=np.linspace(0,100,1001) if reverse is None else [0])
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def main():
    old=json.loads((ROOT/'reversible-conversion-v1/summary.json').read_text());oldarrays=np.load(ROOT/'reversible-conversion-v1/trajectories.npz')
    out=ROOT/'reversible-conversion-v2';out.mkdir(exist_ok=False)
    files=('reversible-conversion-protocol.md','reversible-conversion-amendment.md','reversible_conversion.py','refine_reversible_conversion.py')
    save(out/'manifest.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},original_summary_sha256=hashlib.sha256((ROOT/'reversible-conversion-v1/summary.json').read_bytes()).hexdigest()))
    arrays={};controls=[];matched={}
    for previous in old['controls']:
        cfg=previous['config'];name=previous['name'];sol=solve(cfg);arrays[name]=sol.y;population=sol.y[4]**2+sol.y[5]**2
        if name.startswith('no-conversion'):
            matched[(cfg['delta0'],cfg['P0'])]=sol.y
            controls.append(dict(name=name,config=cfg,max_population=float(population.max()),passed=bool(population.max()<1e-12)))
        else:
            error=float(np.max(abs(population-np.sin(cfg['kappa0']*sol.t)**2)))
            controls.append(dict(name=name,config=cfg,population_error=error,passed=error<1e-8))
    rows=[];negatives=[]
    for previous in old['runs']:
        cfg=previous['config'];name=previous['name'];sol=solve(cfg);arrays[name]=sol.y
        terms=quantities(sol.y,cfg['kappa0'],cfg['A'],cfg['delta0'])
        energy=float(np.max(abs(terms['total']-terms['total'][0]))/abs(terms['total'][0]));norm=float(np.max(abs(terms['norm']-1)))
        reverse=solve(cfg,reverse=sol.y[:,-1]);arrays[name+'-reversed']=reverse.y[:,-1]
        back=float(np.max(abs(reverse.y[:,-1]-sol.y[:,0])))
        qdifference=float(np.max(abs(sol.y[0]-matched[(cfg['delta0'],cfg['P0'])][0])))
        identity=cfg['delta0']!=0 or (qdifference<1e-7 and np.max(abs(terms['cross']))<1e-8)
        rows.append(dict(name=name,config=cfg,relative_energy_drift=energy,norm_drift=norm,reversal_error=back,
                         peak_receiving_population=float(terms['converted_population'].max()),final_receiving_population=float(terms['converted_population'][-1]),
                         initial_energy=float(terms['total'][0]),final_energy=float(terms['total'][-1]),
                         source_energy_change=float(terms['source'][-1]-terms['source'][0]),final_interaction_energy=float(terms['interaction'][-1]),
                         max_environment_difference_from_no_conversion=qdifference,symmetric_no_force_identity=bool(identity),
                         original_endpoint_difference=float(np.max(abs(sol.y[:,-1]-oldarrays[name][:,-1]))),
                         passed=energy<1e-8 and norm<1e-8 and back<1e-7 and identity))
        if cfg['delta0']!=0:
            wrong=solve(cfg,omit=True);key=name+'-omitted-reaction';arrays[key]=wrong.y
            wt=quantities(wrong.y,cfg['kappa0'],cfg['A'],cfg['delta0'])['total'];error=float(np.max(abs(wt-wt[0]))/abs(wt[0]))
            negatives.append(dict(name=key,config=cfg,relative_energy_drift=error,detected=error>1e-5))
    np.savez_compressed(out/'trajectories.npz',time=np.linspace(0,100,1001),**arrays)
    result=dict(runs=rows,controls=controls,negative_controls=negatives,budgets=old['budgets'],required_active_mass_kg=old['required_active_mass_kg'],passed=all(r['passed'] for r in rows+controls) and all(r['detected'] for r in negatives))
    save(out/'summary.json',result)
    print(json.dumps(dict(passed=result['passed'],max_energy_drift=max(r['relative_energy_drift'] for r in rows),max_norm_drift=max(r['norm_drift'] for r in rows),max_reversal_error=max(r['reversal_error'] for r in rows),max_original_endpoint_difference=max(r['original_endpoint_difference'] for r in rows))))
    if not result['passed']:raise RuntimeError('Refined RC-1 gate failed; archive preserved')


if __name__=='__main__':main()
