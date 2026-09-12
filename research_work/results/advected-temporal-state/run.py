from pathlib import Path
import json,csv
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp, quad
HERE=Path(__file__).resolve().parent
def matrices(K,R,C,eta):
    n=1+C/K;a=R/(K*n);d=1-a
    return n,d,np.array([[1.,a],[0.,d]]),np.array([[1/n,-R/(K*n*n)],[0.,eta/(n*n)]])
rows=[];errors=[]
for K in [.05,.5,2.]:
    for R in [.01,.22,1.]:
        for C in [.001,.01,.22,1.]:
            for eta in [1.,.5]:
                n,d,M,B=matrices(K,R,C,eta);assert abs(d)>1e-12
                A=np.linalg.solve(M,B);expected=np.array([1/n,eta/(n*n*d)])
                err=float(max(abs(np.sort(np.linalg.eigvals(A))-np.sort(expected))));errors.append(err);assert err<1e-10
                local=eta/(n*d)
                rows.append(dict(K=K,R=R,C=C,eta=eta,time_matrix_determinant=d,local_radiation_speed=1.,local_state_characteristic=float(local),state_energy_flux_speed_local=eta,
                    capture_depletes_state=d>0,forward_and_within_local_c=bool(0<local<=1+1e-12),coincident=bool(abs(local-1)<1e-12),eigenvalue_error=err))
with (HERE/'states.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
K=.5;R=C=.22;n,d,M,B=matrices(K,R,C,1.);A=np.linalg.solve(M,B);v=1/n;off=A[0,1]
assert abs(A[1,1]-v)<1e-12 and abs(off)>0
jordan=[]
for k in [1,10,100,1000]:
    predicted=np.exp(-1j*k*v)*np.array([-1j*k*off,1.])
    numerical=expm(-1j*k*A)@np.array([0.,1.]);err=float(max(abs(predicted-numerical)));assert err<1e-8
    jordan.append(dict(wavenumber=k,time=1.,radiation_perturbation_amplitude=float(abs(numerical[0])),total_amplitude=float(np.linalg.norm(numerical)),independent_check_error=err))
_,d0,M0,B0=matrices(.5,.51,.01,1.);assert abs(d0)<1e-12
capture=[]
for K in [.5,2.]:
    R0=.22;C0=1.;Gamma=1.;invariant=R0*(K+C0)
    crossing=(-K+np.sqrt(K*K+4*invariant))/2
    def rhs(t,y):
        R,C,D=y;d=1-R/(K+C);dc=-Gamma*C/d
        return [-R*dc/(K+C),dc,Gamma*C]
    def event(t,y): return y[1]-y[0]
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,[0,20],[R0,C0,0.],events=event,rtol=1e-12,atol=1e-14,max_step=.01)
    assert sol.success and len(sol.t_events[0])==1
    t_exact=quad(lambda C:(1-invariant/(K+C)**2)/(Gamma*C),crossing,C0,epsabs=1e-12)[0]
    time_error=abs(sol.t_events[0][0]-t_exact)
    budget_error=float(max(abs(sol.y.sum(axis=0)-(R0+C0))))
    state_error=float(abs(sol.y_events[0][0][1]-crossing))
    assert time_error<1e-8 and state_error<1e-9 and budget_error<1e-9
    capture.append(dict(K=K,R_initial=R0,C_initial=C0,Gamma=Gamma,
        crossing_C=float(crossing),crossing_time=float(t_exact),time_check_error=float(time_error),
        crossing_state_error=state_error,reference_budget_error=budget_error,
        positive_singular_C=float(np.sqrt(invariant)-K) if invariant>K*K else None,
        implication='Capture reaches the defective equal-speed boundary from the initially regular C>R domain.'))
result=dict(states=len(rows),max_eigenvalue_error=max(errors),capture_domain_tests=capture,by_eta={str(eta):dict(n=sum(r['eta']==eta for r in rows),capture_reversal=sum(r['eta']==eta and not r['capture_depletes_state'] for r in rows),forward_within_c=sum(r['eta']==eta and r['forward_and_within_local_c'] for r in rows),coincident=sum(r['eta']==eta and r['coincident'] for r in rows)) for eta in [1.,.5]},
    representative=[r for r in rows if r['K']==.5 and r['R']==.22 and r['C']==.01],
    singular_state=dict(K=.5,R=.51,C=.01,determinant=float(np.linalg.det(M0))),
    coincident_state=dict(K=.5,R=.22,C=.22,off_diagonal=float(off),fourier=jordan),
    scope='Local characteristic and capture audit of a specified constitutive closure, not observational evidence')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(json.dumps(result,indent=2))
