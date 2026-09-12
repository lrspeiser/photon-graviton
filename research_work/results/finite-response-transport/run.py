from pathlib import Path
import csv
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent

def principal(K,R,T,M,sigma):
    n=1+M/K;v=1/n
    return np.array([[v,0,-R/(K*n*n)-R*sigma/(K*n**3)],
                     [0,v,-T/(K*n*n)+R*sigma/(K*n**3)],
                     [0,0,sigma/(n*n)]])

states=[]
for K in [.05,.5,2.]:
    for R in [.01,.22,1.]:
        for M in [0.,.1,1.]:
            for sigma in [1.,.5]:
                T=.01;n=1+M/K;A=principal(K,R,T,M,sigma)
                error=float(max(abs(np.sort(np.linalg.eigvals(A))-np.sort([1/n,1/n,sigma/n**2]))))
                assert error<1e-12
                defective=bool(M==0 and sigma==1)
                condition=None;diag_error=None
                if not defective:
                    w=sigma/n**2;v=1/n
                    V=np.array([[1,0,-A[0,2]/(v-w)], [0,1,-A[1,2]/(v-w)], [0,0,1.]])
                    V=V/np.linalg.norm(V,axis=0)
                    diag_error=float(np.max(abs(A@V-V@np.diag([v,v,w]))))
                    condition=float(np.linalg.cond(V));assert diag_error<1e-10
                states.append(dict(K=K,R=R,T=T,M=M,sigma=sigma,
                    radiation_and_T_local_speed=1.,M_flux_local_speed=sigma,
                    M_characteristic_local_speed=sigma/n,defective=defective,
                    eigenvalue_error=error,eigenvector_condition=condition,diagonalization_error=diag_error))
with (HERE/'characteristics.csv').open('w',newline='',encoding='utf8') as f:
    writer=csv.DictWriter(f,fieldnames=list(states[0]),lineterminator='\n');writer.writeheader();writer.writerows(states)

A=principal(.5,.22,.01,0,1.);N=A-np.eye(3);assert np.max(abs(N@N))==0
jordan=[]
for k in [1,10,100,1000]:
    initial=np.array([0.,0.,1.])
    exact=np.exp(-1j*k)*(initial-1j*k*N@initial)
    direct=expm(-1j*k*A)@initial
    error=float(max(abs(exact-direct)));assert error<1e-8
    jordan.append(dict(k=k,gain=float(np.linalg.norm(direct)),error=error))

def cell(K,R0,M0,Gamma,gated,T0=.01):
    def rhs(t,y):
        R,T,M,D=y;n=1+M/K;b=R/(K*n)
        availability=max(0.,1-b) if gated else 1.
        Q=T-M*availability
        return [-b*Q,(b-1)*Q-Gamma*T,Q,Gamma*T]
    def negative(t,y):return y[1]+1e-9
    negative.terminal=True;negative.direction=-1
    sol=solve_ivp(rhs,[0,20],[R0,T0,M0,0.],method='DOP853',rtol=1e-11,atol=1e-13,max_step=.025,events=negative)
    assert sol.success
    initial_sum=R0+T0+M0;initial_Rn=R0*(1+M0/K)
    budget=float(max(abs(sol.y.sum(axis=0)-initial_sum)))
    invariant=float(max(abs(sol.y[0]*(1+sol.y[2]/K)-initial_Rn)))/max(1.,initial_Rn)
    minimum=float(sol.y.min());failed=bool(len(sol.t_events[0]))
    assert budget<1e-9 and invariant<1e-9
    if gated:assert not failed and minimum>=-1e-9
    return dict(K=K,R0=R0,T0=T0,M0=M0,Gamma=Gamma,gated=gated,
        negative_T_event=failed,time_reached=float(sol.t[-1]),minimum=minimum,
        budget_error=budget,Rn_scaled_error=invariant,
        R_final=float(sol.y[0,-1]),T_final=float(sol.y[1,-1]),M_final=float(sol.y[2,-1]),D_final=float(sol.y[3,-1]))

cells=[cell(K,R,M,Gamma,gated) for K in [.05,.5,2.] for R in [.01,.22,1.] for M in [0.,.1] for Gamma in [0.,1.] for gated in [False,True]]
controls=[cell(K,.22,0,1,True,T0=0) for K in [.05,.5,2.]]
assert all(r['T_final']==r['M_final']==r['D_final']==0 for r in controls)
with (HERE/'cells.csv').open('w',newline='',encoding='utf8') as f:
    writer=csv.DictWriter(f,fieldnames=list(cells[0]),lineterminator='\n');writer.writeheader();writer.writerows(cells)
boundary=[]
K=.5;R=.22;T=0.;M=.1;sigma=.5;n=1+M/K;b=R/(K*n)
Q=T-M*max(0.,1-b)
for gradient in [0.,.1,1.]:
    M_t=Q-sigma*gradient/n**2
    original=b*M_t-Q
    reduced=-principal(K,R,T,M,sigma)[1,2]*gradient+(b-1)*Q
    error=abs(original-reduced);assert error<1e-12
    boundary.append(dict(K=K,R=R,T=T,M=M,sigma=sigma,M_x=gradient,
        T_t=float(original),independent_error=float(error),nonnegative_boundary_direction=bool(original>=0)))
assert boundary[-1]['T_t']<0
result=dict(characteristic_states=len(states),defective_same_speed=sum(r['defective'] for r in states),
    half_speed_max_eigenvector_condition=max(r['eigenvector_condition'] for r in states if r['sigma']==.5),
    fourier=jordan,cell_runs=len(cells),raw_negative_runs=sum(r['negative_T_event'] for r in cells if not r['gated']),
    gated_negative_runs=sum(r['negative_T_event'] for r in cells if r['gated']),
    max_budget_error=max(r['budget_error'] for r in cells),max_Rn_scaled_error=max(r['Rn_scaled_error'] for r in cells),
    no_seed_controls=controls,spatial_positivity_boundary=boundary,
    scope='Constitutive local-characteristic and homogeneous positivity audit; no spatial or observational redshift claim.')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(result,indent=2))
