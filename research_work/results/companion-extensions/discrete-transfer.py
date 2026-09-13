"""Exact killed pure-death energy ladder and joint width/survival requirements."""
from pathlib import Path
import json,math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import binom
P=Path(__file__).resolve().parent

def exact(n,gap,tau,r):
 q=math.exp(-(1+r)*tau)
 logt=math.log1p(-r/(1+r)*(-math.expm1(-(1+r)*tau)))
 logT=n*logt;T=math.exp(logT);p=q/math.exp(logt)
 C=gap/r*(-math.expm1(logT)) if r else n*gap*(-math.expm1(-tau))
 beam=T*n*gap*p;H=n*gap-beam-C
 return dict(survival=T,conditional_energy_eV=n*gap*p,conditional_relative_width=math.sqrt((1-p)/(n*p)),companion_energy_eV=C,removed_energy_eV=H,beam_energy_eV=beam,p=p,zero_energy_conditional_probability=(1-p)**n)

if __name__=='__main__':
 out=dict(scope='Known continuous-time Markov energy ladder applied to proposed fixed-gap transfer and removal; not a derived interaction or arrival-time law',checks=[],joint_designs=[],drift_comparisons=[])
 # Direct master equation, including a separate removed state and energy accounts.
 for n in [5,20,50]:
  for r in [0.,.01,.5]:
   gap=.1;tau=.7
   def rhs(t,y):
    prob=y[:n+1];k=np.arange(n+1);d=-(1+r)*k*prob;d[:-1]+=k[1:]*prob[1:]
    return np.r_[d,r*np.dot(k,prob),gap*np.dot(k,prob),r*gap*np.dot(k*k,prob)]
   y=np.zeros(n+4);y[n]=1
   sol=solve_ivp(rhs,[0,tau],y,rtol=1e-11,atol=1e-13)
   a=exact(n,gap,tau,r);k=np.arange(n+1);expected=a['survival']*binom(n,k)*a['p']**k*(1-a['p'])**(n-k)
   error=max(float(np.max(abs(sol.y[:n+1,-1]-expected))),abs(sol.y[-2,-1]-a['companion_energy_eV'])/(n*gap),abs(sol.y[-1,-1]-a['removed_energy_eV'])/(n*gap))
   assert sol.success and error<1e-9
   assert abs(sol.y[:n+2,-1].sum()-1)<1e-10
   out['checks'].append(dict(initial_quanta=n,removal_ratio=r,max_normalized_error=error))
 znear=math.expm1(.0002488993286382367*30.660139)
 for z in [znear,.1,1.]:
  for w in [1e-3,1e-5]:
   n=math.ceil(z/w**2);p=1/(1+z);r=math.expm1(-math.log(.9)/n)/(1-p)
   tau=math.log((1+r*(1-p))/p)/(1+r)
   a=exact(n,2/n,tau,r)
   assert abs(a['survival']-.9)<1e-10 and a['conditional_relative_width']<=w*(1+1e-10)
   out['joint_designs'].append(dict(z=z,illustrative_width=w,initial_quanta=n,optical_gap_eV=2/n,max_removal_ratio=r,path_integral=tau,energy_result=a))
  n=200000000;gap=1e-8;tau=math.log1p(z);r=-math.log(.9)*gap/(2*z/(1+z));a=exact(n,gap,tau,r)
  out['drift_comparisons'].append(dict(z=z,ratio_from_previous_drift=r,survival_difference_from_point9=a['survival']-.9,conditional_energy_difference_from_drift_eV=a['conditional_energy_eV']-2/(1+z),exact_width=a['conditional_relative_width']))
 (P/'discrete-transfer-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps(out['joint_designs'],indent=2))
