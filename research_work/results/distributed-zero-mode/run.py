from pathlib import Path
import json, hashlib
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
records=[]; probes=[]
grid=np.linspace(0,80,801)
for U in [0.,.003,.006]:
    sols=[]
    for tol in [1e-9,1e-11]:
        sol=solve_ivp(lambda t,y:[y[1],U/y[0]**2],[0,80],[1.,0.],method='DOP853',rtol=tol,atol=tol/100,dense_output=True)
        assert sol.success
        sols.append(sol)
    n,v=sols[-1].sol(grid)
    diff=float(np.max(abs(sols[0].sol(grid)-sols[1].sol(grid))))
    energy=U/n+v*v/2
    err=float(np.max(abs(energy-U))/(U if U else 1))
    if U:
        t_exact=(np.sqrt(n*(n-1))+np.arcsinh(np.sqrt(n-1)))/np.sqrt(2*U)
        time_error=float(np.max(abs(t_exact-grid)))
    else:
        assert np.all(n==1) and np.all(v==0)
        time_error=0.
    assert diff<1e-6 and err<1e-7 and time_error<1e-5
    records.append(dict(U0=U,n_final=float(n[-1]),rate_final=float(v[-1]),asymptotic_rate=float(np.sqrt(2*U)),photon_energy_final=float(U/n[-1]),field_energy_final=float(v[-1]**2/2),relative_energy_error=err,convergence_max_abs=diff,implicit_solution_time_error=time_error))
    background=sols[-1]
    def ray(te):
        def rhs(s,y):
            nn,vv=background.sol(y[0]);return [nn,-vv]
        r=solve_ivp(rhs,[0,2],[te,0.],method='DOP853',rtol=1e-11,atol=1e-13)
        assert r.success and r.y[0,-1]<80
        return r.y[:,-1]
    for te in [1.,5.,20.]:
        to,logratio=ray(te)
        predicted=background.sol(to)[0]/background.sol(te)[0]
        for h in [1e-4,5e-5]:
            J=(ray(te+h)[0]-ray(te-h)[0])/(2*h)
            assert abs(J/predicted-1)<1e-7
            for momentum in [1.,2.]:
                wi=momentum/background.sol(te)[0]
                wo=wi*np.exp(logratio)
                assert abs(wi/wo/predicted-1)<1e-8
                probes.append(dict(U0=U,launch=te,arrival=float(to),momentum=momentum,launch_offset=h,carrier_stretch=float(wi/wo),event_stretch=float(J),endpoint_n_ratio=float(predicted)))
out=dict(protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),backgrounds=records,probes=probes,checks_passed=True,scope='Homogeneous finite periodic 3D sector only; reference-clock predictions, no astronomical fit or gravity completion')
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps(records,indent=2))
