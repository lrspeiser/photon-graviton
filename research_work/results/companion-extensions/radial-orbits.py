"""Radial anisotropy with fixed lens-required masses and constant stellar M/L."""
from pathlib import Path
import sys,hashlib,json
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from scipy.linalg import cho_factor,cho_solve
P=Path(__file__).resolve().parent;ROOT=P.parents[2];L=P.parent/'isotropic-galaxy-transfer'
sys.path.insert(0,str(P.parent/'slacs-component-refit'))
from model import ComponentModel,G,ARCSEC
def read(folder,name='results.json'):return json.loads((P.parent/folder/name).read_text())
constant=json.loads((L/'capacity-reference-exact-lens-optics-results.json').read_text())
free={(r['Name'],r['model']):r for r in json.loads((L/'capacity-reference-optics-results.json').read_text())['rows']}
data={r['Name']:r for r in read('slacs-resolved-input-audit')['systems']}
profiles={r['Name']:r for r in read('slacs-light-profile-audit')['rows']}
pilot={r['Name']:r for r in read('slacs-motion-lensing-pilot')['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
geo={r['Name']:r for r in read('lensing-data-readiness','conditional-geometry.json')}
cfg=read('slacs-outer-bin-check','protocol.json');fit=json.loads((L/'third-radiation-retention-results.json').read_text())['models']['attenuated']
out=dict(scope='Known radial-anisotropy family, fixed Re transition, fixed exact-lens mass and constant M/L; endpoints fitted to inner stars, no positive distribution function derived',rows=[],summary=[])
cache={};mu,w=np.polynomial.legendre.leggauss(96)
for old in constant['rows']:
    name=old['Name'];item=data[name];dl=old['geometry']['angular_Dl_Mpc']*1000;ac=old['capture_scale_kpc']
    if name not in cache:
        a=pilot[name]['scale_a_kpc']*dl/(geo[name]['conditional_Dl_Mpc']*1000)
        edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*dl/ARCSEC
        psf=item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
        components=[dict(R=q['R_arcsec']*dl/ARCSEC,n=q['n'],amp=q['amp_at_R'],bn=q['bn']) for q in profiles[name]['components']]
        cache[name]=ComponentModel(a,edges,psf,0,.5,1,20,components)
    model=cache[name];r=model.r;nu=model.nu
    Re=profiles[name]['computed_equal_area_half_light_arcsec']*dl/ARCSEC
    x=r/ac;t=x[:,None]*mu;b2=1+x[:,None]**2*(1-mu*mu);b=np.sqrt(b2)
    tau=fit['k0_per_kpc']*ac*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
    J=np.exp(-tau)@w/2
    rho=fit['source_C_before_retention_Msun_kpc3']*constant['C0_multiplier']*old['retention_mapping']['eta']*J/(1+x*x)**2
    cm=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(r*r*rho,r,initial=0))
    stellar_fraction=np.array([model.mass_fraction(v) for v in r])
    model.forces=np.array([G*1e11*stellar_fraction/r**2,G*cm/r**2])
    mass=old['mass_Msun'];force=model.forces[0]*mass/1e11+model.forces[1]
    y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);fac=cho_factor(cov[:-1,:-1]);cw=cho_solve(fac,cov[:-1,-1]);csd=np.sqrt(cov[-1,-1]-cov[-1,:-1]@cw)
    def predict(q):
        b0,bi=q;beta=b0+(bi-b0)*r*r/(r*r+Re*Re)
        factor=np.exp(2*b0*np.log(r/model.a)+(bi-b0)*np.log((r*r+Re*Re)/(model.a*model.a+Re*Re)))
        pressure=-cumulative_trapezoid((nu*force*factor)[::-1],r[::-1],initial=0)[::-1]/factor
        v2=np.trapezoid(r*r*pressure*(model.W-beta*model.T),r,axis=1)/model.den
        assert min(v2)>0 and np.isfinite(v2).all()
        return np.sqrt(v2)
    baseline=predict([old['beta'],old['beta']]);drift=float(max(abs(baseline-np.array(old['predicted_stellar_vrms']))))
    assert drift<1e-6
    for beta in [-1.,0.,.4]:
        coef=model.coefficients(beta)
        assert max(abs(predict([beta,beta])-np.sqrt(coef[0]*mass/1e11+coef[1])))<1e-8
    def loss(q):
        e=y[:-1]-predict(q)[:-1];return float(e@cho_solve(fac,e))
    starts=[[a,b] for a in [-1.5,-.5,0,.4] for b in [-1.5,-.5,0,.4]]+[[old['beta'],old['beta']]]
    fits=[minimize(loss,q,bounds=[cfg['constant_beta_bounds']]*2,method='L-BFGS-B',options={'ftol':1e-11,'maxiter':500}) for q in starts]
    good=[f for f in fits if f.success and np.isfinite(f.fun)];assert good
    opt=min(good,key=lambda f:f.fun);pred=predict(opt.x);outer=pred[-1]+cw@(y[:-1]-pred[:-1])
    assert opt.fun<=old['inner_chi2']+1e-5
    bounds=cfg['constant_beta_bounds'];fr=free[name,old['model']]
    out['rows'].append(dict(Name=name,population=old['retention_mapping']['population'],beta0=float(opt.x[0]),beta_infinity=float(opt.x[1]),transition_Re_kpc=Re,
        endpoint_boundary=bool(any(min(abs(v-bounds[0]),abs(v-bounds[1]))<1e-4 for v in opt.x)),stellar_mass_Msun=mass,
        free_constant_inner_chi2=fr['inner_chi2'],exact_constant_inner_chi2=old['inner_chi2'],radial_inner_chi2=float(opt.fun),
        free_outer_residual=fr['outer_conditional_standardized_residual'],exact_constant_outer_residual=old['outer_conditional_standardized_residual'],radial_outer_residual=float((y[-1]-outer)/csd),
        observed_stellar_vrms=y.tolist(),predicted_stellar_vrms=pred.tolist(),constant_limit_max_velocity_difference_kms=drift,optimizer_successes=len(good)))
    if opt.x[1]>.4499:
        expanded_bounds=[cfg['constant_beta_bounds'],[-2,.95]]
        seeds=[opt.x.tolist()]+[[v,b] for v in [-.2,.2,opt.x[0]] for b in [.7,.94]]
        trials=[minimize(loss,q,bounds=expanded_bounds,method='L-BFGS-B',options={'ftol':1e-11,'maxiter':500}) for q in seeds]
        successful=[f for f in trials if f.success];assert successful
        expanded=min(successful,key=lambda f:f.fun);ev=predict(expanded.x);eo=ev[-1]+cw@(y[:-1]-ev[:-1])
        assert expanded.fun<=opt.fun+1e-5
        out['rows'][-1]['expanded_outer_bound_followup']=dict(beta0=float(expanded.x[0]),beta_infinity=float(expanded.x[1]),
            outer_bound=.95,outer_boundary=bool(abs(expanded.x[1]-.95)<1e-4),inner_chi2=float(expanded.fun),outer_standardized_residual=float((y[-1]-eo)/csd))
    print(name,old['retention_mapping']['population'],opt.x,opt.fun,flush=True)
for pop in ['Chabrier','Salpeter']:
    rr=[r for r in out['rows'] if r['population']==pop]
    out['summary'].append(dict(population=pop,endpoint_boundary_count=sum(r['endpoint_boundary'] for r in rr),
        inner_chi2={k:sum(r[k] for r in rr) for k in ['free_constant_inner_chi2','exact_constant_inner_chi2','radial_inner_chi2']},
        outer_residual_square_sum={k:sum(r[k]**2 for r in rr) for k in ['free_outer_residual','exact_constant_outer_residual','radial_outer_residual']}))
    out['summary'][-1]['expanded_outer_bound_followup_inner_chi2_sum']=sum(r.get('expanded_outer_bound_followup',{}).get('inner_chi2',r['radial_inner_chi2']) for r in rr)
    out['summary'][-1]['expanded_outer_bound_followup_outer_residual_square_sum']=sum(r.get('expanded_outer_bound_followup',{}).get('outer_standardized_residual',r['radial_outer_residual'])**2 for r in rr)
paths=[L/'capacity-reference-optics-results.json',L/'capacity-reference-exact-lens-optics-results.json',L/'third-radiation-retention-results.json',P.parent/'slacs-component-refit/model.py',P.parent/'slacs-resolved-fit/model.py']
paths += [P.parent/f/n for f,n in [('slacs-resolved-input-audit','results.json'),('slacs-light-profile-audit','results.json'),('slacs-motion-lensing-pilot','results.json'),('lensing-data-readiness','conditional-geometry.json'),('slacs-outer-bin-check','protocol.json')]]
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
out['known_profile_reference']='https://arxiv.org/abs/0705.4109'
assert len(out['rows'])==12
(P/'radial-orbits-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2))
