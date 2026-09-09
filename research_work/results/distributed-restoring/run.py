from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
HERE=Path(__file__).resolve().parent
U=.003
records=[]
for m in [0.,.01,.03,.1]:
    end=8*np.pi/(m if m else .01)
    def rhs(t,y):return [y[1],U/y[0]**2-m*m*(y[0]-1)]
    def maximum(t,y):return y[1]
    maximum.direction=-1
    def minimum(t,y):return y[1]
    minimum.direction=1
    runs=[]
    for tol in [1e-9,1e-11]:
        sol=solve_ivp(rhs,[0,end],[1.,0.],method='DOP853',rtol=tol,atol=tol/1e4,max_step=.5,dense_output=True,events=[maximum,minimum])
        assert sol.success
        runs.append(sol)
    times=np.linspace(0,end,5001)
    n,v=runs[-1].sol(times)
    energy=v*v/2+m*m*(n-1)**2/2+U/n
    err=float(np.max(abs(energy-U))/U)
    diff=float(np.max(abs(runs[0].sol(times)-runs[1].sol(times))))
    assert n.min()>0 and err<1e-7 and diff<1e-6
    r=dict(m=m,end=end,relative_energy_error=err,convergence_max_abs=diff,n_min=float(n.min()),n_max_sampled=float(n.max()),rate_min=float(v.min()),rate_max=float(v.max()))
    if m:
        tmax=float(runs[-1].t_events[0][0]);nmax=float(runs[-1].y_events[0][0,0])
        returns=runs[-1].t_events[1];treturn=float(returns[returns>1e-6][0])
        exact=(1+np.sqrt(1+8*U/(m*m)))/2
        nreturn,vreturn=runs[-1].sol(treturn)
        assert abs(nmax/exact-1)<1e-7 and abs(treturn/(2*tmax)-1)<1e-7
        assert abs(nreturn-1)<1e-7 and v.min()<0
        r.update(first_max_time=tmax,first_max_n=nmax,analytic_max_n=float(exact),return_time=treturn,return_photon_energy=float(U/nreturn),return_field_energy=float(vreturn*vreturn/2+m*m*(nreturn-1)**2/2))
    records.append(r)
(HERE/'results.json').write_text(json.dumps(dict(protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),cases=records,checks_passed=True),indent=2)+'\n',newline='\n')
print(json.dumps(records,indent=2))
