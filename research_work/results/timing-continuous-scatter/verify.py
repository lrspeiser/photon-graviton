from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from integration import averaged_likelihood
HERE=Path(__file__).resolve().parent
x=np.linspace(0,5,81);row=.02+np.exp(-.5*((x-2.7)/.3)**2)
checks=[]
for mu in [0.,.013,2.731,4.99,5.]:
 for sigma in [0.,1e-8,1e-4,.01,.03,.1,.6]:
    value=float(averaged_likelihood(x,row[None,:],np.array([mu]),sigma)[0])
    if sigma==0: reference=float(np.interp(mu,x,row))
    else:
      low=max(-12.,(x[0]-mu)/sigma);high=min(12.,(x[-1]-mu)/sigma)
      points=((x-mu)/sigma);points=points[(points>low)&(points<high)]
      reference=quad(lambda u:np.interp(mu+sigma*u,x,row)*np.exp(-u*u/2)/np.sqrt(2*np.pi),low,high,points=points,epsabs=1e-12,epsrel=1e-11,limit=200)[0]/(ndtr(high)-ndtr(low))
    error=abs(value-reference);assert error<2e-10,(mu,sigma,value,reference)
    checks.append(dict(mu=mu,sigma=sigma,absolute_error=error))
out=dict(scope='Exact integration check for declared interpolant; not accuracy of underlying event curves or coverage',checks=checks,maximum_absolute_error=max(r['absolute_error'] for r in checks))
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(checks=len(checks),maximum_absolute_error=out['maximum_absolute_error'])))
