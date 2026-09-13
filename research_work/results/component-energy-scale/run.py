"""Separate shared deposited-source normalization from radial shape mismatch."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import brentq,minimize_scalar
H=Path(__file__).resolve().parent;R=H.parents[2];p=H.parent/'component-loading-cepheids/results.json';d=json.loads(p.read_text(encoding='utf-8'))
for name,value in d['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
outputs=[]
for result in d['results']:
 bins=result['bins'];y=np.array([v['inferred_jeans_speed'] for v in bins]);base=np.array([v['ordinary_speed']**2 for v in bins]);delta=np.array([v['predicted_balance_speed']**2 for v in bins])-base;assert np.all(delta>0)
 def residual(k):return np.sqrt(base+k*delta)-y
 # Mean-square loss is strictly convex in k for positive y,base,delta.
 def derivative(k):return float(np.mean(delta*(1-y/np.sqrt(base+k*delta))))
 hi=1.
 while derivative(hi)<0:hi*=2
 k=brentq(derivative,0,hi) if derivative(0)<0 else 0.
 opt=minimize_scalar(lambda v:float(np.mean(residual(v)**2)),bounds=(0,hi),method='bounded',options={'xatol':1e-11});assert opt.success and abs(opt.x-k)<1e-6
 def balance(k):v=residual(k);return float(v.max()+v.min())
 upper=hi
 while balance(upper)<0:upper*=2
 minimax=brentq(balance,0,upper) if balance(0)<0 else 0.
 needed=(y*y-base)/delta
 np.testing.assert_allclose(np.sqrt(base+needed*delta),y,rtol=1e-13)
 outputs.append(dict(refined=result['refined'],best_rms_scale=k,rms_kms=float(np.sqrt(np.mean(residual(k)**2))),bias_kms=float(np.mean(residual(k))),minimax_scale=minimax,minimum_worst_residual_kms=float(np.max(abs(residual(minimax)))),required_scale_range=[float(needed.min()),float(needed.max())],bins=[dict(bin=v['bin'],radius_mean=v['radius_mean'],required_scale=float(needed[i]),observed_proxy=float(y[i]),rms_fit_prediction=float(y[i]+residual(k)[i]),rms_fit_residual=float(residual(k)[i]),minimax_prediction=float(y[i]+residual(minimax)[i])) for i,v in enumerate(bins)]))
 print(json.dumps({k:v for k,v in outputs[-1].items() if k!='bins'}))
out=dict(scope='One-parameter calibration on exposed Cepheid proxies; not energy supply prediction or blind validation',hashes={str(v.relative_to(R)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [Path(__file__),p]},results=outputs)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
