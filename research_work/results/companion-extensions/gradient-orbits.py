"""Joint stellar M/L gradient and radial anisotropy, calibrated lens angles."""
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.linalg import cho_factor,cho_solve
P=Path(__file__).resolve().parent;ROOT=P.parents[2];L=P.parent/'isotropic-galaxy-transfer'
parser=argparse.ArgumentParser()
parser.add_argument('--slope-constrained',action='store_true')
args=parser.parse_args()
if args.slope_constrained:
    from orbit_density import deproject_slope
    previous={(r['Name'],r['population']):r for r in json.loads((P/'gradient-orbits-results.json').read_text())['rows']}
    slope_cache={}
    constraint_x=np.geomspace(1e-6,100,4097)
    check_x=np.geomspace(1e-6,100,8193)
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
separate_gradient={(r['Name'],r['population']):r for r in json.loads((P/'stellar-gradient-results.json').read_text())['rows']}
separate_radial={(r['Name'],r['population']):r for r in json.loads((P/'radial-orbits-results.json').read_text())['rows']}
out=dict(scope='Fixed companions; lens-calibrated stellar mass; gradient plus radial anisotropy fitted to inner motions; no positive distribution function or stability established',rows=[],summary=[])
out['slope_constrained']=args.slope_constrained
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
    def evaluate(q):
        h,b0,bi=q;weight=h*meanH/(1+h*meanH)
        mass=1e11*needed/((1-weight)*s0+weight*sh)
        force=mass/1e11*((1-weight)*model.forces[0]+weight*model.forces[1])+model.forces[2]
        beta=b0+(bi-b0)*r*r/(r*r+Re*Re)
        factor=np.exp(2*b0*np.log(r/model.a)+(bi-b0)*np.log((r*r+Re*Re)/(model.a*model.a+Re*Re)))
        pressure=-cumulative_trapezoid((nu*force*factor)[::-1],r[::-1],initial=0)[::-1]/factor
        v2=np.trapezoid(r*r*pressure*(model.W-beta*model.T),r,axis=1)/model.den
        assert min(v2)>0 and np.isfinite(v2).all() and 1e7<mass<1e14
        pred=np.sqrt(v2);err=y[:-1]-pred[:-1];score=float(err@cho_solve(fac,err))
        return score,pred,mass,weight
    key=(name,old['retention_mapping']['population'])
    sg=separate_gradient[key];sr=separate_radial[key]
    rb=sr.get('expanded_outer_bound_followup',{})
    radial_seed=[0,rb.get('beta0',sr['beta0']),rb.get('beta_infinity',sr['beta_infinity'])]
    grad_seed=[sg['h'],sg['beta'],sg['beta']]
    radial_score=rb.get('inner_chi2',sr['radial_inner_chi2'])
    radial_outer=rb.get('outer_standardized_residual',sr['radial_outer_residual'])
    rg=evaluate(radial_seed);gg=evaluate(grad_seed)
    radial_drift=abs(rg[0]-radial_score);gradient_drift=abs(gg[0]-sg['gradient_inner_chi2'])
    assert radial_drift<1e-4 and gradient_drift<1e-4
    bounds=[(-.8,9),(-2,.45),(-2,.95)]
    starts=[radial_seed,grad_seed]+[[h,b0,bi] for h in [-.7,1,8] for b0 in [-.5,.3] for bi in [-1,.4,.94]]
    if args.slope_constrained:
        if name not in slope_cache:
            gamma=deproject_slope(profiles[name],constraint_x,512)[1]
            fine_gamma=deproject_slope(profiles[name],check_x,512)[1]
            nmax=max(c['n'] for c in profiles[name]['components'] if c['amp_at_R']>0)
            assert nmax>1
            slope_cache[name]=(gamma,fine_gamma,1-1/nmax)
        gamma,fine_gamma,central_gamma=slope_cache[name]
        fraction=constraint_x**2/(1+constraint_x**2)
        def central_bound(bi):
            return min(.45,central_gamma/2,float(np.min((gamma/2-bi*fraction)/(1-fraction))))
        def physical(q):
            h,t,bi=q;cap=central_bound(bi);assert cap>-2
            return np.array([h,-2+t*(cap+2),bi])
        def mapped(q):
            h,b0,bi=q
            return [h,float(np.clip((b0+2)/(central_bound(bi)+2),0,1)),bi]
        prev=previous[key];starts.append([prev['h'],prev['beta0'],prev['beta_infinity']])
        fits=[minimize(lambda q:evaluate(physical(q))[0],mapped(q),bounds=[bounds[0],(0,1),bounds[2]],method='L-BFGS-B',options={'ftol':1e-11,'maxiter':600}) for q in starts]
        for f in fits:f.x=physical(f.x)
    else:
        fits=[minimize(lambda q:evaluate(q)[0],q,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-11,'maxiter':600}) for q in starts]
    good=[f for f in fits if f.success and np.isfinite(f.fun)];assert good
    opt=min(good,key=lambda f:f.fun);score,pred,mass,weight=evaluate(opt.x)
    assert score<=(radial_score if args.slope_constrained else min(radial_score,sg['gradient_inner_chi2']))+1e-4
    h,b0,bi=opt.x;outer=pred[-1]+cw@(y[:-1]-pred[:-1])
    lens_error=abs(dr*(mass/1e11*((1-weight)*s0+weight*sh)+dc)/target-1)
    assert lens_error<1e-10
    out['rows'].append(dict(Name=name,population=key[1],h=float(h),beta0=float(b0),beta_infinity=float(bi),Re_kpc=Re,
        stellar_mass_Msun=mass,central_to_outer_ML_ratio=float(1+h),
        boundaries={key:bool(min(abs(v-lo),abs(v-hi))<1e-4) for key,v,(lo,hi) in zip(['h','beta0','beta_infinity'],opt.x,bounds)},
        combined_inner_chi2=score,combined_outer_residual=float((y[-1]-outer)/csd),
        radial_inner_chi2=radial_score,radial_outer_residual=radial_outer,
        gradient_inner_chi2=sg['gradient_inner_chi2'],gradient_outer_residual=sg['gradient_outer_standardized_residual'],
        free_inner_chi2=old['inner_chi2'],free_outer_residual=old['outer_conditional_standardized_residual'],
        observed_stellar_vrms=y.tolist(),predicted_stellar_vrms=pred.tolist(),exact_angle_fractional_residual=lens_error,
        radial_limit_chi2_difference=radial_drift,gradient_limit_chi2_difference=gradient_drift,
        optimizer_successes=len(good),optimizer_attempts=len(fits),successful_objectives=[float(f.fun) for f in good],
        failed_optimizer_messages=[str(f.message) for f in fits if not f.success]))
    if args.slope_constrained:
        fine_beta=b0+(bi-b0)*check_x**2/(1+check_x**2)
        margin=fine_gamma-2*fine_beta
        assert min(margin)>-1e-6 and central_gamma-2*b0>=-1e-10
        out['rows'][-1].update(unconstrained_inner_chi2=prev['combined_inner_chi2'],unconstrained_outer_residual=prev['combined_outer_residual'],central_asymptotic_gamma=central_gamma,central_constraint_margin=float(central_gamma-2*b0),slope_constraint_boundary=bool(abs(b0-central_bound(bi))<1e-5),refined_minimum_slope_margin=float(min(margin)),refined_minimum_radius_Re=float(check_x[np.argmin(margin)]))
    print(name,key[1],opt.x,score,out['rows'][-1]['combined_outer_residual'],flush=True)
for pop in ['Chabrier','Salpeter']:
    rr=[r for r in out['rows'] if r['population']==pop]
    out['summary'].append(dict(population=pop,systems=len(rr),boundary_counts={k:sum(r['boundaries'][k] for r in rr) for k in ['h','beta0','beta_infinity']},
        inner_chi2={k:sum(r[k+'_inner_chi2'] for r in rr) for k in ['free','gradient','radial','combined']},
        outer_residual_square_sum={k:sum(r[k+'_outer_residual']**2 for r in rr) for k in ['free','gradient','radial','combined']}))
paths=[L/'capacity-reference-optics-results.json',L/'capacity-reference-exact-lens-optics-results.json',L/'third-radiation-retention-results.json',P.parent/'slacs-component-refit/model.py',P.parent/'slacs-resolved-fit/model.py']
paths += [P.parent/f/n for f,n in [('slacs-resolved-input-audit','results.json'),('slacs-light-profile-audit','results.json'),('slacs-motion-lensing-pilot','results.json'),('lensing-data-readiness','conditional-geometry.json'),('slacs-outer-bin-check','protocol.json')]]
paths += [P/'stellar-gradient-results.json',P/'radial-orbits-results.json',P/'gradient-orbits-protocol.md']
if args.slope_constrained:paths += [P/'gradient-orbits-results.json',P/'orbit_density.py',P/'constrained-gradient-orbits-protocol.md']
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
assert len(out['rows'])==12
(P/('constrained-gradient-orbits-results.json' if args.slope_constrained else 'gradient-orbits-results.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2))
