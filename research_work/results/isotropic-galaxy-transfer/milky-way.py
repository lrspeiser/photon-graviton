"""No-refit transfer of the exact deposited profile to two MW observables."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
source=HERE.parent/'joint-galaxy-audit/milky-way-predictions.json'
old=json.loads(source.read_text()); fit=json.loads((HERE/'results.json').read_text())
G=4.30091727003628e-6; rd=2.6
r=np.linspace(0,30,4001);mu,w=leggauss(96)
out=dict(R_disk_kpc=rd,scale_source='https://arxiv.org/abs/astro-ph/0510520',scores=[],input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,HERE/'results.json']})
rows=[]
for model in ['baryons','power','transparent_control','attenuated']:
    if model in fit['models']:
        p=fit['models'][model];a=rd*p['scale_to_disk'];x=r/a
        shape=(1+x*x)**-2
        if model=='attenuated':
            t=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu**2);B=np.sqrt(B2)
            tau=p['k0_per_kpc']*a*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
            J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1)
        else:J=np.ones_like(x)
        mass=4*np.pi*p['C_Msun_kpc3']*cumulative_trapezoid(shape*J*r*r,r,initial=0)
    for variant in ['I','II']:
        selected=[d for d in old if d['baryons']==variant and d['model']==('sparc_transfer' if model=='power' else 'baryons')]
        assert len(selected)==81
        local=[]
        for d in selected:
            pred=d['predicted'];R=d['R_kpc'];z=d['z_kpc'];rr=float(np.hypot(R,z))
            if model in fit['models']:
                M=float(np.interp(rr,r,mass))
                if d['observable']=='vc_kms':pred=float(np.sqrt(pred**2+G*M/R))
                else:pred+=M*abs(z)/(2*np.pi*rr**3*1e6)
            row=dict(d,model=model,predicted=pred,residual=pred-d['observed'])
            rows.append(row);local.append(row)
        for subset in ['inner','outer','rotation_all','vertical_provisional']:
            sel=[d for d in local if (d['observable']=='vc_kms' if subset=='rotation_all' else d['split']==subset)]
            e=np.array([d['residual'] for d in sel])
            out['scores'].append(dict(model=model,baryons=variant,subset=subset,n=len(sel),RMSE=float(np.sqrt(np.mean(e*e))),bias=float(np.mean(e))))
for fn,obj in [('milky-way-results.json',out),('milky-way-predictions.json',rows)]:
    (HERE/fn).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out['scores'],indent=2))
