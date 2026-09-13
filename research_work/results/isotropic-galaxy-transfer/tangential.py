"""Sign of exact beta=-1/2 inversion for fixed fitted lens deposits."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'slacs-component-refit'))
from model import deproject,G,ARCSEC
fit=json.loads((HERE/'results.json').read_text())['models']['attenuated']
systems=[d for d in json.loads((HERE/'lensing-results.json').read_text())['rows'] if d['model']=='attenuated']
profiles={d['Name']:d for d in json.loads((HERE.parent/'slacs-light-profile-audit/results.json').read_text())['rows']}
geo={d['Name']:d for d in json.loads((HERE.parent/'lensing-data-readiness/conditional-geometry.json').read_text())}
mu,w=np.polynomial.legendre.leggauss(128)
xcheck=np.geomspace(1e-3,1e3,1601);rows=[]
for d in systems:
    name=d['Name'];dl=geo[name]['conditional_Dl_Mpc']*1000;ac=d['capture_scale_kpc']
    comps=[dict(R=c['R_arcsec']*dl/ARCSEC,n=c['n'],amp=c['amp_at_R'],bn=c['bn']) for c in profiles[name]['components']]
    outputs=[]
    for n in [2001,4001]:
        x=np.geomspace(1e-5,1e5,n);r=ac*x;lr=np.log(r)
        nu,ltotal=deproject(r,comps,256)
        slope=np.log(nu[1]/nu[0])/np.log(r[1]/r[0])
        mb=4*np.pi*(nu[0]*r[0]**3/(3+slope)+cumulative_trapezoid(nu*r*r,r,initial=0))
        masscheck=mb[-1]/ltotal;mb=mb/mb[-1]*d['mass_Msun']
        t=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu*mu);B=np.sqrt(B2)
        tau=fit['k0_per_kpc']*ac*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
        J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1)
        rho=fit['C_Msun_kpc3']*J/(1+x*x)**2
        md=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(rho*r*r,r,initial=0))
        g=G*(mb+md)/r**2
        ell=np.gradient(np.log(rho/r),lr,edge_order=2)
        m=np.gradient(np.log(g),lr,edge_order=2)
        Q=ell**2-ell+np.gradient(ell,lr,edge_order=2)-ell*m
        q=np.interp(np.log(xcheck),np.log(x),Q)
        outputs.append(dict(n=n,Q=q,minimum_Q=float(q.min()),minimum_radius_kpc=float(ac*xcheck[q.argmin()]),stellar_mass_integral_ratio=float(masscheck)))
    coarse,fine=outputs;change=float(np.max(abs(fine['Q']-coarse['Q'])))
    bad=fine['Q'] < -max(.001,2*change)
    rows.append(dict(Name=name,deposit_beta=-.5,minimum_Q=fine['minimum_Q'],minimum_radius_kpc=fine['minimum_radius_kpc'],coarse_minimum_Q=coarse['minimum_Q'],max_Q_change=change,negative_region_resolved=bool(bad.any()),negative_radius_min_kpc=float(ac*xcheck[bad][0]) if bad.any() else None,negative_radius_max_kpc=float(ac*xcheck[bad][-1]) if bad.any() else None,stellar_mass_integral_ratios=[v['stellar_mass_integral_ratio'] for v in outputs]))
out=dict(scope='Fixed beta=-1/2 sign test; not formation or global stability',rows=rows)
(HERE/'tangential-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
