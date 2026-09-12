from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
slow=json.loads((HERE/'results.json').read_text())['runs']
baseline=json.loads((HERE.parent/'radiation-powered-initiation/results.json').read_text())['runs']
K=2.;sigma=.5;epsilon=.01
rows=[]
for tau in [1.,10.,100.]:
    def rhs(x,y):
        FR,FT,FM,cap,time=y
        assert FM<sigma*K
        n=1/(1-FM/(sigma*K));R=n*FR;T=n*FT;M=K*(n-1)
        Q=(T-M)/tau;gamma=10*(np.exp(-x*x)+np.exp(-(x-10)**2))
        return [-epsilon*R,-Q-gamma*T+epsilon*R,Q,gamma*T,n]
    sol=solve_ivp(rhs,[0,10],[.22,0,0,0,0],rtol=1e-11,atol=1e-13,max_step=.01)
    assert sol.success
    flux_error=float(max(abs(sol.y[:4].sum(axis=0)-.22)));assert flux_error<1e-10
    FR,FT,FM,cap,travel=sol.y[:,-1];n=1/(1-FM/(sigma*K));z=1/n-1
    source=baseline if tau==1 else slow;label='initiation_capture' if tau==1 else 'tau'+str(int(tau))
    compare=[]
    for cells in [80,160,320]:
        sim=next(r for r in source if r['case']==label and r['cells']==cells)
        p=next(p for p in sim['probes'] if p['emission']==20)
        error=abs(p['measured_z']-z)
        if cells==320:assert error<1e-5
        compare.append(dict(cells=cells,evolved_z=p['measured_z'],absolute_difference=float(error)))
    rows.append(dict(response_time=tau,stationary_z=float(z),stationary_coordinate_stretch=1.,travel_time=float(travel),
        photon_survival=float(np.exp(-epsilon*travel)),photon_flux_out=float(FR),T_flux_out=float(FT),M_flux_out=float(FM),
        deposit_growth_rate=float(cap),flux_error=flux_error,comparison=compare))
(HERE/'stationary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(rows,indent=2))
