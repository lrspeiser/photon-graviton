"""Self-heating reservoir for the protected-state candidate."""
from pathlib import Path
import json,math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
B=1.01**3;g=.001;theta0=.01

def run(eps,C,cooling,tol=1e-8):
 def rhs(A,back):
  def f(t,y):
   p0,pa,pb,theta=y[:4];assert theta>0
   n=1/math.expm1(eps/theta);down=n+1;up=n
   flow=down*pa-up*pb;escape=C*cooling*(theta-theta0)
   return [-A*p0+(back+g)*pa,A*p0-(back+g)*pa-flow,flow,(eps*flow-escape)/C,A*p0-back*pa,g*pa,escape]
  return f
 y0=[1.,0,0,theta0,0,0,0]
 pump=solve_ivp(rhs(1.,B),[0,20],y0,method='Radau',rtol=tol,atol=tol*.01)
 dark=solve_ivp(rhs(0.,0.),[0,1000],pump.y[:,-1],method='Radau',rtol=tol,atol=tol*.01)
 assert pump.success and dark.success
 def record(y):
  p0,pa,pb,theta,net,radiation,escape=y;U=pa+(1-eps)*pb;bath=C*(theta-theta0)
  assert abs(p0+pa+pb-1)<1e-8 and min(p0,pa,pb)>-1e-10
  assert abs(net-U-bath-radiation-escape)<1e-8
  assert abs(bath+escape-eps*pb)<1e-8
  return dict(protected_probability=float(pb),bath_temperature=float(theta),stored_energy=float(U),bath_energy_gain=float(bath),net_optical_input=float(net),radiated_energy=float(radiation),net_exported_bath_energy=float(escape),thermal_up_down_ratio=math.exp(-eps/theta))
 return dict(relaxation_energy_fraction=eps,dimensionless_heat_capacity=C,cooling_rate=cooling,pump=record(pump.y[:,-1]),dark=record(dark.y[:,-1]))

if __name__=='__main__':
 out=dict(scope='Constant heat capacity per site and linear ambient cooling are hypotheses; no astronomical rate or bath identified',initial_temperature=theta0,rows=[],refinements=[],closed_bath_stationary=[])
 for eps in [.1,.5,.9]:
  for C in [.1,1.,10.,100.]:
   for cooling in [0.,.01,1.]:
    r=run(eps,C,cooling);out['rows'].append(r)
    if eps==.5 and C in [.1,10.]:
     fine=run(eps,C,cooling,1e-10)
     error=max(abs(r[stage][k]-fine[stage][k]) for stage in ['pump','dark'] for k in r[stage])
     assert error<1e-6
     out['refinements'].append(dict(heat_capacity=C,cooling_rate=cooling,max_difference=error))
   pb=brentq(lambda p:p-1/(1+(1+B+g)*math.exp(-eps/(theta0+eps*p/C))),0.,1.)
   # 1+(B+g)/A with A=1 is the normalization coefficient above.
   out['closed_bath_stationary'].append(dict(eps=eps,heat_capacity=C,protected_probability=pb,bath_temperature=theta0+eps*pb/C))
 out['capacity_design_example']=dict(eps=.5,target_protected_probability=.9,target_up_down_ratio=.001,minimum_heat_capacity=.5*.9/(.5/math.log(1000)-theta0))
 (P/'finite-bath-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps([r for r in out['rows'] if r['relaxation_energy_fraction']==.5 and r['cooling_rate']==0],indent=2))
