"""Analytic central support obstruction evaluated on fitted galaxy profiles."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
fit=json.loads((HERE/'results.json').read_text())['models']['attenuated']
energy=json.loads((HERE/'energy-results.json').read_text())
lenses=json.loads((HERE/'lensing-results.json').read_text())
k0=fit['k0_per_kpc'];s=fit['scale_to_disk']
critical=(-3*np.pi/2+np.sqrt((3*np.pi/2)**2+48))/2
mu,w=leggauss(128)
rows=[]
inputs=[(d['galaxy'],'SPARC',d['split'],d['capture_scale_kpc']) for d in energy['rows']]
inputs+=[('Milky Way','MW','transfer',2.6*s)]
inputs+=[(d['Name'],'SLACS','transfer',d['capture_scale_kpc']) for d in lenses['rows'] if d['model']=='attenuated']
for name,sample,split,a in inputs:
    K=k0*a;coef=K*K/6+np.pi*K/4-2
    x=1e-4;t=x*mu;B2=1+x*x*(1-mu*mu);B=np.sqrt(B2)
    tau=K*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
    # Divide out the exact central exp(-pi K/4), avoiding tiny absolute densities.
    ratio=.5*np.sum(np.exp(-tau+np.pi*K/4)*w)/(1+x*x)**2
    numeric=(ratio-1)/x**2
    assert abs(numeric-coef)<1e-3*(1+abs(coef))
    rows.append(dict(name=name,sample=sample,split=split,capture_scale_kpc=a,K=K,central_quadratic_coefficient=coef,direct_coefficient=numeric,isotropic_central_support_excluded=bool(coef>0)))
summary=[]
for sample in ['SPARC','MW','SLACS']:
    rr=[d for d in rows if d['sample']==sample]
    summary.append(dict(sample=sample,n=len(rr),excluded_by_central_condition=sum(d['isotropic_central_support_excluded'] for d in rr)))
out=dict(critical_K=float(critical),critical_disk_scale_kpc=float(critical/(k0*s)),summary=summary,rows=rows)
(HERE/'support-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
