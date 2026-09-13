"""Paraxial source-to-observer moments of independent angular kicks."""
from pathlib import Path
import json,math
import numpy as np
from scipy.integrate import quad
P=Path(__file__).resolve().parent;PC=3.085677581491367e16;c=299792458.;ARC=180*3600/math.pi
old=json.loads((P/'quantum-bridge-results.json').read_text());angles=json.loads((P/'angular-outlet-results.json').read_text())
d=next(r['step_over_source_radius'] for r in angles['rows'] if r['input_occupation']==1e-12 and r['target_reverse_forward']==.9)
R=6.957e8/PC;alpha=old['alpha_per_Mpc']/1e6;delta=1e-8
out=dict(scope='Independent zero-mean kicks, paraxial endpoint geometry, mean fixed-gap event rate; prescribed kernels, no angular occupation feedback',rows=[],geometry_checks=[])
for D in [old['distance_Mpc']*1e6,1e8,1e9]:
 for s0 in [1/206264.80624709636,1.,1000.]:
  for kernel in ['fixed_at_one_pc','tracks_apparent_source_radius']:
   def f(t,kind):
    s=s0*math.exp(t);E=2*math.exp(-alpha*(s-s0));nu=alpha*E/delta
    theta=d*R/(s if kernel=='tracks_apparent_source_radius' else 1.)
    q=nu*theta*theta
    if kind=='image':return q*(s/D)**2*s
    if kind=='delay':return q*s*(D-s)/D*s*PC/(2*c)
    return q*s
   upper=math.log(D/s0)
   values={k:quad(lambda t:f(t,k),0,upper,epsabs=1e-60,epsrel=1e-9)[0] for k in ['image','delay','raw']}
   t=np.linspace(0,upper,8193)
   for k in values:
    approx=np.trapezoid([f(x,k) for x in t],t)
    assert abs(approx/values[k]-1)<1e-4
   if kernel=='tracks_apparent_source_radius':
    exact=(d*R/D)**2*2*(-math.expm1(-alpha*(D-s0)))/delta
    assert abs(values['image']/exact-1)<1e-10
   out['rows'].append(dict(distance_pc=D,start_pc=s0,kernel=kernel,apparent_direction_rms_arcsec=math.sqrt(values['image'])*ARC,unconditioned_direction_rms_proxy_arcsec=math.sqrt(values['raw'])*ARC,mean_geometric_delay_seconds=values['delay'],mean_net_steps=2*(-math.expm1(-alpha*(D-s0)))/delta))
# Independent geometric derivation for a single kick and endpoint constraint.
for x in [.1,.5,.9]:
 kick=1e-4;initial=-(1-x)*kick;final=x*kick
 displacement=initial*x+final*(1-x)
 delay=.5*(initial**2*x+final**2*(1-x))
 assert abs(displacement)<1e-18 and abs(delay/(.5*kick*kick*x*(1-x))-1)<1e-12
 out['geometry_checks'].append(dict(fraction=x,endpoint_residual=displacement,dimensionless_delay=delay))
(P/'observer-transport-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
for r in out['rows']:
 if r['distance_pc']==old['distance_Mpc']*1e6:print(r)
