"""Conditional incoherent oscillator ensemble times Gaussian inelastic rate."""
from pathlib import Path
import importlib.util,json,math
import numpy as np
from scipy.integrate import quad
P=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('spatial',P/'spatial-response.py');sp=importlib.util.module_from_spec(spec);spec.loader.exec_module(sp)

def response(E,p,eta,low,high,tol=1e-10):
 a,b=math.log(low/E),math.log(high/E)
 def integrand(t):
  v=math.exp(2*t)
  return math.exp((p+1)*t)/((v-1)**2+eta*eta*v)
 val=quad(integrand,a,b,points=([0.] if a<0<b else None),epsabs=tol,epsrel=tol,limit=300)[0]
 return E**(p-3)*val

if __name__=='__main__':
 out=dict(scope='Inverse-designed distribution of independent damped resonances; not a complete optical or energy-conserving interaction',size_inverse_eV=1e6,gap_eV=sp.gap,cases=[])
 for p in [1.,2.,3.]:
  for eta in [.1,1.,3.]:
   for low,high in [(1e-5,1e5),(1e-3,1e3)]:
    normalization=sp.rate(1.,1e6)[0]*response(1.,p,eta,low,high)
    rows=[]
    for E in np.geomspace(.01,100.,17):
     r=response(E,p,eta,low,high);ref=response(E,p,eta,low,high,1e-12)
     assert abs(ref/r-1)<1e-8
     rows.append(dict(energy_eV=float(E),relative_fractional_loss=float(sp.rate(E,1e6)[0]*r/normalization),response_refinement_relative_change=ref/r-1))
    out['cases'].append(dict(spectral_power=p,damping_ratio=eta,minimum_resonance_eV=low,maximum_resonance_eV=high,max_fractional_departure_from_1eV=max(abs(r['relative_fractional_loss']-1) for r in rows),rows=rows))
 # For eta=2,p=2, infinite dimensionless integral u^2/(u^2+1)^2=pi/4.
 value=response(1.,2.,2.,1e-8,1e8)
 assert abs(value-math.pi/4)<1.1e-8
 (P/'resonant-spectrum-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps([{k:v for k,v in q.items() if k!='rows'} for q in out['cases'] if q['spectral_power']==2],indent=2))
