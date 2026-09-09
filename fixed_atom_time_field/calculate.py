"""Synthetic fixed-atom propagation-field model; not an observational fit.

Dimensionless time s=gamma*t, distance X=gamma*R/c0, n(today)=1.
eta = radiation energy(today)/(K*gamma**2); photon bath energy is eta/n.
Integrate backward along a ray parametrized by increasing geometric distance.
"""
from pathlib import Path
import json, math
import numpy as np
from scipy.integrate import solve_ivp

P=Path(__file__).resolve().parent
kappa=0.000077315
grid=[0.01,0.1,1.,3.,7.]
rows=[]
for eta in [0.,1e-6,1e-4,0.001,0.1]:
    # y=(n, dn/ds, s). Along ray ds/dX=-n.
    sol=solve_ivp(lambda X,y:[-y[0]*y[1],-eta/y[0],-y[0]],
        [0,7],[1,1,0],method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True,max_step=.02)
    assert sol.success
    check=sol.sol(np.linspace(0,7,2001))
    energy=.5*check[1]**2+eta/check[0]
    values=[]
    for X in grid:
        n,v,s=sol.sol(X)
        S=1/n
        values.append(dict(X=X,distance_Mly=X/kappa,stretch=float(S),
            exponential_stretch=math.exp(X),
            relative_difference_from_exponential=float(S/math.exp(X)-1),
            emission_time_gamma_t=float(s),field_velocity_gamma=float(v)))
    rows.append(dict(eta=eta,max_relative_energy_error=float(np.max(abs(energy/(.5+eta)-1))),
        maximum_possible_stretch=None if eta==0 else 1+.5/eta,values=values))

# Independent finite-pulse arrival test in prescribed affine n=1+s.
# Source at X=0.1. Arrival time sr=(1+se)*exp(X)-1.
X=.1; se=math.exp(-X)-1; ds=1e-6
arrivals=[]
for start in [se,se+ds]:
    def end(t,y):return y[0]-X
    end.terminal=True;end.direction=1
    sol=solve_ivp(lambda t,y:[1/(1+t)],[start,.2],[0],events=end,
        rtol=1e-12,atol=1e-14,max_step=.001)
    arrivals.append(float(sol.t_events[0][0]))
stretch=(arrivals[1]-arrivals[0])/ds

T=2.72548; S=3000/T
result=dict(status='Synthetic derivation checks only; unchanged atomic physics is stipulated, not microscopically proven',
    kappa_per_Mly=kappa,gamma_per_year=kappa/1e6,
    free_field_exponential_max_relative_error=max(abs(v['relative_difference_from_exponential']) for v in rows[0]['values']),
    conservative_bath_cases=rows,
    finite_pulse=dict(X=X,stretch=stretch,expected=math.exp(X),relative_error=stretch/math.exp(X)-1),
    illustrative_3000K_to_CMB=dict(stretch=S,maximum_eta_to_reach_stretch=1/(2*(S-1)),
        minimum_field_to_radiation_energy_today=S-1,
        caveat='Necessary bound in this homogeneous conserved-bath toy model, not sufficient for an exponential law; no 3000 K source is inferred'))
(P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
