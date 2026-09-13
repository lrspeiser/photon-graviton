"""Common s=2 mass redistribution; inner stellar fits, outer/lens predictions."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.integrate import quad
from scipy.linalg import cho_factor,cho_solve
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'slacs-component-refit'))
from model import ComponentModel,G,C,ARCSEC
def read(folder,file='results.json'):return json.loads((HERE.parent/folder/file).read_text())
data=read('slacs-resolved-input-audit');base=read('slacs-motion-lensing-pilot')
geo={r['Name']:r for r in read('lensing-data-readiness','conditional-geometry.json')}
profiles={r['Name']:r for r in read('slacs-light-profile-audit')['rows']}
pilot={r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
cfg=read('slacs-outer-bin-check','protocol.json');par=base['parameters']
A,p=par['A'],par['p'];astar=par['a_star_m_s2']*3.085677581491367e19/1e6
prior={r['Name']:r for r in read('slacs-outer-bin-check')['rows'] if r['model']=='empirical_extra'}
rows=[]
for item in data['systems']:
    name=item['Name']
    if name not in prior:continue
    dl=geo[name]['conditional_Dl_Mpc']*1000;ratio=geo[name]['conditional_Dls_over_Ds'];a=pilot[name]['scale_a_kpc']
    edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*dl/ARCSEC
    psf=item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
    components=[dict(R=q['R_arcsec']*dl/ARCSEC,n=q['n'],amp=q['amp_at_R'],bn=q['bn']) for q in profiles[name]['components']]
    model=ComponentModel(a,edges,psf,A,p,astar,20,components)
    y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);fac=cho_factor(cov[:-1,:-1])
    crossw=cho_solve(fac,cov[:-1,-1]);csd=np.sqrt(cov[-1,-1]-cov[-1,:-1]@crossw)
    rt=model.rt
    def original_extra(r,mass):
        if r<=rt:return A*astar*(G*mass*model.mass_fraction(r)/r**2/astar)**p
        return A*astar*(G*mass*model.mass_fraction(rt)/rt**2/astar)**p*(rt/r)**2
    for s in [1.,2.]:
        model.forces[1]=np.array([original_extra(r/s,1e11)/s**2 for r in model.r])
        def objective(q):
            e=y[:-1]-model.predict(np.exp(q[0]),q[1],True)[:-1]
            return float(e@cho_solve(fac,e))
        fits=[minimize(objective,[np.log(1e11),beta],method='L-BFGS-B',
            bounds=[np.log(cfg['mass_Msun_bounds']),cfg['constant_beta_bounds']],
            options={'ftol':1e-11,'maxiter':1000}) for beta in cfg['starts_beta']]
        good=[f for f in fits if f.success and np.isfinite(f.fun)];assert good
        fit=min(good,key=lambda f:f.fun);mass=float(np.exp(fit.x[0]));beta=float(fit.x[1])
        pred=model.predict(mass,beta,True);outer=pred[-1]+crossw@(y[:-1]-pred[:-1])
        def lens_residual(b):
            def integrand(t):
                r=b/np.cos(t);g=G*mass*model.mass_fraction(r)/r**2+original_extra(r/s,mass)/s**2
                return g*b/np.cos(t)
            edge=np.arccos(b/(s*rt)) if b<s*rt else None
            bend=4/C**2*quad(integrand,0,np.pi/2,points=[edge] if edge is not None else None,epsabs=1e-6,epsrel=1e-7,limit=200)[0]
            return ratio*bend-b/dl
        lo=a*1e-7;hi=a*1e5
        assert lens_residual(lo)*lens_residual(hi)<0
        angle=brentq(lens_residual,lo,hi,xtol=1e-9)/dl*ARCSEC
        if s==1:
            assert abs(fit.fun-prior[name]['chi2'])<.01 and abs(angle-prior[name]['predicted_angle_arcsec'])<.002
        rows.append(dict(Name=name,radial_scale=s,mass_Msun=mass,beta=beta,
            orbit_boundary=bool(np.min(abs(beta-np.array(cfg['constant_beta_bounds'])))<1e-5),
            inner_chi2=float(fit.fun),outer_observed_kms=float(y[-1]),outer_prediction_kms=float(outer),
            outer_conditional_standardized_residual=float((y[-1]-outer)/csd),
            lens_prediction_arcsec=angle,lens_catalog_arcsec=pilot[name]['catalog_SIE_arcsec'],
            lens_fractional_residual=angle/pilot[name]['catalog_SIE_arcsec']-1,
            optimizer_successes=sum(bool(f.success) for f in fits)))
        print(name,s,fit.fun,angle,flush=True)
summary=[]
for s in [1.,2.]:
    rr=[r for r in rows if r['radial_scale']==s]
    summary.append(dict(radial_scale=s,inner_chi2_sum=sum(r['inner_chi2'] for r in rr),
        lens_fractional_rms=float(np.sqrt(np.mean([r['lens_fractional_residual']**2 for r in rr]))),
        outer_conditional_residual_square_sum=sum(r['outer_conditional_standardized_residual']**2 for r in rr),
        orbit_boundary_count=sum(r['orbit_boundary'] for r in rr)))
out=dict(scope='Fixed redistribution and inner-bin nuisance refit on exposed galaxies; no source derivation',summary=summary,rows=rows)
(HERE/'redistribution-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summary,indent=2))
