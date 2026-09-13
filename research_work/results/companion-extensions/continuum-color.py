"""Positive receiving spectra in the existing point-response interaction."""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
P=Path(__file__).resolve().parent

def response(E,q,low,high,n):
 top=min(E,high)
 if top<=low:return 0.,0.
 x,w=leggauss(n);a,b=np.log(low),np.log(top);t=(a+b)/2+(b-a)*x/2;weights=w*(b-a)/2
 gap=np.exp(t);base=np.exp((q+2)*t)
 alpha=float(np.sum(weights*base*(E-gap)**3))
 derivative=float(np.sum(weights*base*3*(E-gap)**2))
 return alpha,E*derivative/alpha

if __name__=='__main__':
 out=dict(scope='Forward creation into an energy-independent positive gap spectrum, heavy-store point response; no inverse population or spatial overlap',energy_unit='eV; arbitrary common rate normalization cancels',cases=[])
 for q in [-2,-1,0,1]:
  for high in [.001,1.,100.]:
   reference,_=response(1.,q,1e-8,high,128);rows=[]
   for E in [.01,.03,.1,.3,1.,3.,10.]:
    value,slope=response(E,q,1e-8,high,128);fine,fs=response(E,q,1e-8,high,256)
    rows.append(dict(photon_energy_eV=E,relative_fractional_loss=value/reference,logarithmic_slope=slope,refinement_relative_difference=fine/value-1,required_squared_amplitude_relative_to_1eV=reference/value))
    assert slope>=3-1e-12
    assert abs(fine/value-1)<1e-10
   out['cases'].append(dict(spectral_power=q,min_gap_eV=1e-8,max_gap_eV=high,rows=rows))
 # Independent analytic integral for q=0 and a cutoff above E.
 E=.7;lo=1e-8
 primitive=lambda d:E**3*d*d/2-E*E*d**3+3*E*d**4/4-d**5/5
 numerical,_=response(E,0,lo,1.,256)
 assert abs(numerical/(primitive(E)-primitive(lo))-1)<1e-12
 (P/'continuum-color-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps(dict(cases=len(out['cases']),min_slope=min(r['logarithmic_slope'] for c in out['cases'] for r in c['rows']),max_slope=max(r['logarithmic_slope'] for c in out['cases'] for r in c['rows']),ten_eV_relative_loss_range=[min(c['rows'][-1]['relative_fractional_loss'] for c in out['cases']),max(c['rows'][-1]['relative_fractional_loss'] for c in out['cases'])]),indent=2))
