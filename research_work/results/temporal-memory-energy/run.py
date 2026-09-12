from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
HERE=Path(__file__).resolve().parent
resets=[];growth=[]
for shape in ['linear','quadratic']:
    for K in [.05,.5,2.]:
        def U(n):return K*(n-1) if shape=='linear' else K*(n-1)**2/2
        def Up(n):return K if shape=='linear' else K*(n-1)
        for Ri in [.1,1.,10.]:
            ni=1.1
            def rhs(t,y):
                n,R,D=y;ndot=-(n-1);return [ndot,-ndot*R/n,-ndot*(Up(n)-R/n)]
            sol=solve_ivp(rhs,[0,12],[ni,Ri,0.],rtol=2e-11,atol=1e-13,max_step=.02,dense_output=True);assert sol.success
            n,R,D=sol.sol(np.linspace(0,12,1201));exact=U(ni)-U(n)-(Ri*ni/n-Ri)
            error=float(max(abs(D-exact)));assert error<1e-9
            rate=(n-1)*(Up(n)-R/n)
            limit=U(ni)-Ri*(ni-1)
            resets.append(dict(shape=shape,K=K,R_initial=Ri,final_deposit=float(D[-1]),complete_reset_budget=float(limit),minimum_deposit_rate=float(min(rate)),passive_rate_nonnegative_on_sample=bool(min(rate)>=-1e-12),endpoint_budget_check=error))
for K in [.05,.5,2.]:
    for capture in [0.,1.]:
        def rhs(t,y):
            n,R,T,D=y;return [T,-T*R/n,T*(R/n-K)-capture*T,capture*T]
        sol=solve_ivp(rhs,[0,30],[1.,1.,.01,0.],rtol=2e-11,atol=1e-13,max_step=.02,dense_output=True);assert sol.success
        n,R,T,D=sol.sol(np.linspace(0,30,1501));memory=K*(n-1)
        err=float(max(abs(R+T+D+memory-1.01)));assert err<1e-9 and min(min(R),min(T),min(D),min(memory))>-1e-12
        growth.append(dict(K=K,capture=capture,initial_receiving_growth_coefficient=1-K-capture,final_n=float(n[-1]),final_reference_radiation=float(R[-1]),final_receiving=float(T[-1]),final_memory_energy=float(memory[-1]),final_deposit=float(D[-1]),energy_error=err))
cycles=[]
for row in growth:
    if row['capture']!=0:continue
    restored_radiation=row['final_reference_radiation']*row['final_n']
    reset_budget=row['final_memory_energy']-(restored_radiation-row['final_reference_radiation'])
    assert abs(restored_radiation-1)<1e-9 and abs(reset_budget-(.01-row['final_receiving']))<1e-9
    cycles.append(dict(K=row['K'],radiation_after_complete_reset=restored_radiation,net_reset_budget=reset_budget,seed_minus_remaining_receiver=.01-row['final_receiving'],interpretation='Net endpoint budget only; nonnegative instantaneous capture rate is not guaranteed'))
result=dict(reset_diagnostics=resets,growth_diagnostics=growth,closed_cycle_checks=cycles,scope='Illustrative constitutive memory energy; reference ledger, not complete spacetime stress-energy')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(result,indent=2))
