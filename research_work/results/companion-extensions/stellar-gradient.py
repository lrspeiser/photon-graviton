"""Bounded M/L gradients at exact lens angles, fixed reference companions."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.linalg import cho_factor,cho_solve
P=Path(__file__).resolve().parent;ROOT=P.parents[2];L=P.parent/'isotropic-galaxy-transfer'
sys.path.insert(0,str(P.parent/'slacs-component-refit'))
from model import ComponentModel,G,C,ARCSEC
def read(folder,name='results.json'):return json.loads((P.parent/folder/name).read_text())
source=json.loads((L/'capacity-reference-optics-results.json').read_text())
constant=json.loads((L/'capacity-reference-exact-lens-optics-results.json').read_text())
fixed={(r['Name'],r['model']):r for r in constant['rows']}
data={r['Name']:r for r in read('slacs-resolved-input-audit')['systems']}
profiles={r['Name']:r for r in read('slacs-light-profile-audit')['rows']}
pilot={r['Name']:r for r in read('slacs-motion-lensing-pilot')['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
geo={r['Name']:r for r in read('lensing-data-readiness','conditional-geometry.json')}
cfg=read('slacs-outer-bin-check','protocol.json')
fit=json.loads((L/'third-radiation-retention-results.json').read_text())['models']['attenuated']
out=dict(scope='Fixed companion profile; exact catalogue angle calibrates stellar mass, bounded M/L gradient and constant beta fit inner stars; no independent population gradient evidence',rows=[],summary=[])
cache={};mu,w=np.polynomial.legendre.leggauss(96)
for old in source['rows']:
    name=old['Name'];item=data[name];dl=old['geometry']['angular_Dl_Mpc']*1000;dr=old['geometry']['Dls_over_Ds'];ac=old['capture_scale_kpc']
    if name not in cache:
        a=pilot[name]['scale_a_kpc']*dl/(geo[name]['conditional_Dl_Mpc']*1000)
        edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*dl/ARCSEC
        psf=item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
        components=[dict(R=q['R_arcsec']*dl/ARCSEC,n=q['n'],amp=q['amp_at_R'],bn=q['bn']) for q in profiles[name]['components']]
        cache[name]=ComponentModel(a,edges,psf,0,.5,1,20,components)
    model=cache[name];r=model.r;nu=model.nu
    Re=profiles[name]['computed_equal_area_half_light_arcsec']*dl/ARCSEC
    initial=4*np.pi*nu[0]*r[0]**3/model.inner_power
    lum=initial+4*np.pi*cumulative_trapezoid(r*r*nu,r,initial=0)
    H=1/(1+(r/Re)**2)
    core=initial*H[0]+4*np.pi*cumulative_trapezoid(r*r*nu*H,r,initial=0)
    meanH=core[-1]/lum[-1];core/=core[-1]
    cp=PchipInterpolator(np.log(r),core)
    def core_fraction(rr):
        if rr<r[0]:return core[0]*(rr/r[0])**model.inner_power
        if rr>r[-1]:return 1.
        return float(cp(np.log(rr)))
    x=r/ac;t=x[:,None]*mu;b2=1+x[:,None]**2*(1-mu*mu);b=np.sqrt(b2)
    tau=fit['k0_per_kpc']*ac*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
    J=np.exp(-tau)@w/2
    rho=fit['source_C_before_retention_Msun_kpc3']*source['C0_multiplier']*old['retention_mapping']['eta']*J/(1+x*x)**2
    cm=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(r*r*rho,r,initial=0));ci=PchipInterpolator(np.log(r),cm)
    def comp_mass(rr):
        if rr<r[0]:return cm[0]*(rr/r[0])**3
        if rr>r[-1]:return cm[-1]
        return float(ci(np.log(rr)))
    f0=np.array([model.mass_fraction(v) for v in r])
    model.forces=np.array([G*1e11*f0/r**2,G*1e11*core/r**2,G*cm/r**2])
    target=old['lens_catalog_arcsec']/ARCSEC;impact=target*dl
    def bend(fn):return 4*G/C**2*quad(lambda t:fn(impact/np.cos(t))/(impact/np.cos(t)),0,np.pi/2,epsabs=1e-8,epsrel=1e-8,limit=200)[0]
    s0=1e11*bend(model.mass_fraction);sh=1e11*bend(core_fraction);dc=bend(comp_mass)
    needed=target/dr-dc;assert needed>0
    y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);fac=cho_factor(cov[:-1,:-1]);cw=cho_solve(fac,cov[:-1,-1]);csd=np.sqrt(cov[-1,-1]-cov[-1,:-1]@cw)
    grid=np.linspace(*cfg['constant_beta_bounds'],129);interp=PchipInterpolator(grid,np.array([model.coefficients(beta) for beta in grid]),axis=0)
    def evaluate(q,direct=False,pure_core=False):
        h,beta=q;weight=1. if pure_core else h*meanH/(1+h*meanH)
        mass=1e11*needed/((1-weight)*s0+weight*sh)
        coeff=model.coefficients(beta) if direct else interp(beta)
        v2=mass/1e11*((1-weight)*coeff[0]+weight*coeff[1])+coeff[2]
        assert min(v2)>0 and 1e7<mass<1e14
        pred=np.sqrt(v2);err=y[:-1]-pred[:-1];score=float(err@cho_solve(fac,err))
        return score,pred,mass,weight
    baseline=fixed[name,old['model']]
    zero=evaluate([0,baseline['beta']],True)
    assert abs(zero[0]-baseline['inner_chi2'])<1e-4
    assert abs(zero[2]/baseline['mass_Msun']-1)<1e-7
    bounds=[(-.8,9),cfg['constant_beta_bounds']]
    starts=[[h,beta] for h in [-.5,0,1,4,9] for beta in [-1,0,.3]]
    fits=[minimize(lambda q:evaluate(q)[0],q,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-11,'maxiter':500}) for q in starts]
    good=sorted([f for f in fits if f.success],key=lambda f:f.fun);assert good
    final=[minimize(lambda q:evaluate(q,True)[0],f.x,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-11,'maxiter':500}) for f in good[:3]]
    final=[f for f in final if f.success];assert final
    opt=min(final,key=lambda f:f.fun);score,pred,mass,weight=evaluate(opt.x,True)
    assert score<=baseline['inner_chi2']+1e-5
    h,beta=opt.x;outer=pred[-1]+cw@(y[:-1]-pred[:-1])
    lens_error=abs(dr*(mass/1e11*((1-weight)*s0+weight*sh)+dc)/target-1)
    assert lens_error<1e-10
    out['rows'].append(dict(Name=name,population=old['retention_mapping']['population'],h=float(h),central_to_outer_ML_ratio=float(1+h),gradient_scale_Re_kpc=Re,
        gradient_boundary=bool(min(abs(h+.8),abs(h-9))<1e-4),beta=float(beta),beta_boundary=bool(min(abs(beta-cfg['constant_beta_bounds'][0]),abs(beta-cfg['constant_beta_bounds'][1]))<1e-4),
        stellar_mass_Msun=mass,free_constant_ML_inner_chi2=old['inner_chi2'],exact_constant_ML_inner_chi2=baseline['inner_chi2'],gradient_inner_chi2=score,
        free_outer_standardized_residual=old['outer_conditional_standardized_residual'],constant_exact_outer_standardized_residual=baseline['outer_conditional_standardized_residual'],gradient_outer_standardized_residual=float((y[-1]-outer)/csd),
        predicted_stellar_vrms=pred.tolist(),observed_stellar_vrms=y.tolist(),exact_angle_fractional_residual=lens_error,zero_gradient_inner_chi2_reproduction_difference=abs(zero[0]-baseline['inner_chi2']),
        coefficient_interpolation_final_score_difference=abs(evaluate(opt.x)[0]-score)))
    if h>8.9999:
        checks=[]
        for larger in [99.,999.,None]:
            def core_loss(q):return evaluate([0 if larger is None else larger,q[0]],True,larger is None)[0]
            trials=[minimize(core_loss,[v],bounds=[cfg['constant_beta_bounds']],method='L-BFGS-B',options={'ftol':1e-11,'maxiter':500}) for v in cfg['starts_beta']]
            successes=[v for v in trials if v.success];assert successes
            best=min(successes,key=lambda v:v.fun)
            checks.append(dict(h=larger,pure_core_limit=larger is None,beta=float(best.x[0]),inner_chi2=float(best.fun)))
        out['rows'][-1]['large_gradient_followup']=checks
    print(name,old['retention_mapping']['population'],h,beta,score,flush=True)
for pop in ['Chabrier','Salpeter']:
    rr=[r for r in out['rows'] if r['population']==pop]
    out['summary'].append(dict(population=pop,systems=len(rr),gradient_boundary_count=sum(r['gradient_boundary'] for r in rr),beta_boundary_count=sum(r['beta_boundary'] for r in rr),
        inner_chi2={k:sum(r[k] for r in rr) for k in ['free_constant_ML_inner_chi2','exact_constant_ML_inner_chi2','gradient_inner_chi2']},
        outer_residual_square_sum={k:sum(r[k]**2 for r in rr) for k in ['free_outer_standardized_residual','constant_exact_outer_standardized_residual','gradient_outer_standardized_residual']}))
paths=[L/'capacity-reference-optics-results.json',L/'capacity-reference-exact-lens-optics-results.json',L/'third-radiation-retention-results.json',P.parent/'slacs-component-refit/model.py',P.parent/'slacs-resolved-fit/model.py']
paths += [P.parent/f/n for f,n in [('slacs-resolved-input-audit','results.json'),('slacs-light-profile-audit','results.json'),('slacs-motion-lensing-pilot','results.json'),('lensing-data-readiness','conditional-geometry.json'),('slacs-outer-bin-check','protocol.json')]]
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
assert len(out['rows'])==12
(P/'stellar-gradient-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2))
