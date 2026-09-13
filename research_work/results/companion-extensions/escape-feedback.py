"""One-zone stimulated return and escaping-wave energy accounting."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
eps=.5;B=1.01**3;ga=.001

def run(mu,tau,gb,tol=1e-8):
 def rhs(A,reverse):
  def f(t,y):
   p0,pa,pb,n=y[:4];flow=(n+1)*pa-n*pb
   return [-A*p0+(reverse+ga)*pa+gb*pb,A*p0-(reverse+ga)*pa-flow,flow-gb*pb,mu*flow-n/tau,A*p0-reverse*pa,ga*pa+(1-eps)*gb*pb,eps*n/(mu*tau)]
  return f
 def jac(A,reverse):
  def f(t,y):
   _,pa,pb,n=y[:4];d=pa-pb
   return np.array([
    [-A,reverse+ga,gb,0,0,0,0],
    [A,-(reverse+ga)-n-1,n,-d,0,0,0],
    [0,n+1,-n-gb,d,0,0,0],
    [0,mu*(n+1),-mu*n,mu*d-1/tau,0,0,0],
    [A,-reverse,0,0,0,0,0],
    [0,ga,(1-eps)*gb,0,0,0,0],
    [0,0,0,eps/(mu*tau),0,0,0]])
  return f
 initial=[1.,0,0,0,0,0,0]
 pump=solve_ivp(rhs(1.,B),[0,100],initial,method='Radau',jac=jac(1.,B),rtol=tol,atol=tol*.01)
 dark=solve_ivp(rhs(0.,0.),[0,1000],pump.y[:,-1],method='Radau',jac=jac(0.,0.),rtol=tol,atol=tol*.01)
 assert pump.success and dark.success
 def record(y):
  p0,pa,pb,n,inputE,rad,escaped=y;stored=pa+(1-eps)*pb;waves=eps*n/mu
  assert min(p0,pa,pb,n)>-1e-9 and abs(p0+pa+pb-1)<1e-8
  assert abs(inputE-stored-waves-rad-escaped)<1e-8
  return dict(protected_probability=float(pb),bright_probability=float(pa),wave_occupation=float(n),stored_energy=float(stored),wave_energy=float(waves),net_input_energy=float(inputE),other_radiation_energy=float(rad),escaped_relaxation_energy=float(escaped))
 return dict(sites_per_mode=mu,escape_time=tau,protected_decay_rate=gb,pump=record(pump.y[:,-1]),dark=record(dark.y[:,-1]))

out=dict(scope='Homogeneous rate closure with mean escape time; illustrative site/mode ratios and rates, not spatial radiative transfer',rows=[],checks=[])
for mu in [.001,1.,1000.]:
 for tau in [.1,10.,1000.]:
  for gb in [0.,1e-5]:
   row=run(mu,tau,gb);out['rows'].append(row)
   if mu in [.001,1000.] and tau in [.1,1000.]:
    fine=run(mu,tau,gb,1e-10);err=max(abs(row[s][k]-fine[s][k]) for s in ['pump','dark'] for k in row[s]);assert err<1e-5
    out['checks'].append(dict(mu=mu,tau=tau,gb=gb,max_difference=err))
out['stationary_checks']=[]
for row in out['rows']:
 mu,tau,gb=row['sites_per_mode'],row['escape_time'],row['protected_decay_rate']
 def balance(pb):
  n=mu*tau*gb*pb
  pa=pb*(gb+n)/(1+n)
  return (1+B+ga)*pa+(1+gb)*pb-1
 pb=brentq(balance,0,1,xtol=1e-14)
 n=mu*tau*gb*pb;pa=pb*(gb+n)/(1+n);p0=1-pa-pb
 flow=(n+1)*pa-n*pb
 residual=max(abs(-p0+(B+ga)*pa+gb*pb),abs(p0-(B+ga)*pa-flow),abs(flow-gb*pb),abs(mu*flow-n/tau))
 assert residual<1e-9
 out['stationary_checks'].append(dict(mu=mu,tau=tau,gb=gb,protected_probability=pb,residual=residual))
(P/'escape-feedback-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([r for r in out['rows'] if r['protected_decay_rate']==0 and r['escape_time']==1000],indent=2))
