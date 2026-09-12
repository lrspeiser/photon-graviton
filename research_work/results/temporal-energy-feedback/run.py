from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
results=[]
for label,seed,radiation,capture in [('empty_seed',0.,1.,0.),('open_void',.01,1.,0.),('brighter_open_void',.01,10.,0.),('capture_cell',.01,1.,1.)]:
    # n, EM, GW, receiving field, deposits, exported receiving-field energy
    initial=np.array([1.,radiation,.1*radiation,seed,0.,0.])
    def rhs(t,y):
        n,em,gw,field,deposit,export=y;h=field/n
        return [field,-h*em,-h*gw,h*(em+gw)-(1+capture)*field,capture*field,field]
    sol=solve_ivp(rhs,[0,12],initial,rtol=2e-11,atol=1e-13,dense_output=True,max_step=.025)
    assert sol.success
    sample=sol.sol(np.linspace(0,12,1201));total=initial[1:].sum()
    err=float(np.max(abs(sample[1:].sum(axis=0)-total)));assert err<1e-8 and sample[1:].min()>-1e-12
    def arrive(te):
        return brentq(lambda to:quad(lambda t:1/sol.sol(t)[0],te,to,epsabs=1e-11)[0]-1,te,12,xtol=1e-12)
    te=.1;to=arrive(te);second=arrive(te+.001);S=float(sol.sol(to)[0]/sol.sol(te)[0]);eps=1e-5
    finite=(arrive(te+eps)-arrive(te-eps))/(2*eps)
    assert abs(finite-S)/S<1e-5
    final=sample[:,-1]
    assert abs(final[0]-1-final[5])<1e-8 # independently implied by n'=field and export'=field
    results.append(dict(case=label,initial_energy=float(total),final_n=float(final[0]),final_em=float(final[1]),final_gw=float(final[2]),final_field=float(final[3]),final_deposit=float(final[4]),exported_field_energy=float(final[5]),max_energy_error=err,
        test_path_length=1.,emission_time=te,em_arrival=to,gw_arrival=to,relative_messenger_delay=0.,frequency_stretch=S,event_interval_stretch=(second-to)/.001,infinitesimal_event_stretch_check=finite,
        effective_alpha_start=seed,effective_alpha_end=float(final[3]),effective_alpha_max=float(sample[3].max())))
assert results[0]['frequency_stretch']==1. and results[0]['final_deposit']==0.
(HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n');print(json.dumps(results,indent=2))
