from pathlib import Path
import numpy as np,json
from scipy.integrate import solve_ivp
p=.26390611655880836;gamma=7.7315e-11
# Dimensionless tau=gamma*t, V/(K gamma^2)=-n^(2p)/2-beta/n.
# y=(n,dn/dtau,u/(K gamma^2)). Closed homogeneous field+radiation ODE.
rows=[]
for beta in [.001,.1,1.]:
 def rhs(t,y):
  n,v,u=y
  vp=-p*n**(2*p-1)+beta/n**2
  return [v,-vp+u/n,-v*u/n]
 initial=.3;end=1.5
 stop=(end**(1-p)-initial**(1-p))/(1-p)
 t=np.linspace(0,stop,1001)
 sol=solve_ivp(rhs,[0,stop],[initial,initial**p,beta/initial],t_eval=t,rtol=2e-11,atol=2e-12)
 assert sol.success
 n,v,u=sol.y;expected=(initial**(1-p)+(1-p)*t)**(1/(1-p))
 energy=.5*v*v-.5*n**(2*p)-beta/n+u
 err=float(max(abs(n-expected)))
 assert err<1e-7
 rows.append(dict(beta=beta,max_n_error=err,max_energy_error=float(max(abs(energy))),max_photon_invariant_error=float(max(abs(u*n-beta)))))
out=dict(status='Synthetic dynamical reconstruction; no new observations or complete gravity/matter action',p=p,
 potential='V(n)=V0-K gamma^2 n^(2p)/2-U/n; radiation u=U/n; K>0',
 field_equation='K n_ddot+V_prime=u/n',radiation_equation='u_dot=-(n_dot/n)u',
 total_energy='K n_dot^2/2+V+u=V0 along target trajectory',
 restricted_curvature_today_over_gamma2=p*(1-2*p),
 restricted_curvature_assumption='Fixed homogeneous photon invariant U, second derivative of V+U/n; not full coupled photon-field perturbations',
 past_positive_n_branch_endpoint_years_before_today=1/((1-p)*gamma),tests=rows,
 limitations=['Potential reverse-engineered from fitted p,gamma and chosen background radiation U.', 'No atomic protection or dynamical spacetime has been included.', 'Positive restricted potential curvature does not prove full stability.', 'Potential is unbounded below on n>0; needs modification outside finite tested interval for a globally bounded energy.', 'Energy cancellation tunes potential to chosen photon invariant; emission/absorption and density perturbations alter it.', 'Finite past n=0 endpoint is a breakdown of this branch, not a demonstrated origin of the universe.'])
Path(__file__).with_name('dynamics_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
