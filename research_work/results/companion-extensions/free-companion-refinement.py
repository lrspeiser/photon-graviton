"""Free companion capture scale and density strength, matched to free NFW freedoms."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from functools import lru_cache
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
from orbit_density import deproject_slope
nfw=json.loads((P/"free-nfw-results.json").read_text())
stars=json.loads((P/"stellar-only-lens-control-results.json").read_text())
reference=json.loads((P/"orbit-transition-results.json").read_text())
fit_results=json.loads((P/"free-companion-results.json").read_text())
out=dict(scope="Bounded companion capture scale and density strength plus four stellar parameters fitted to all motion bins; lens calibration consumed; no independent prediction or stability proof",rows=[])
cache={};mu,w=np.polynomial.legendre.leggauss(96)
for old in source['rows']:
    if old['retention_mapping']['population']!='Chabrier':continue
    counterpart=next(v for v in source['rows'] if v['Name']==old['Name'] and v['retention_mapping']['population']=='Salpeter')
    assert counterpart['geometry']==old['geometry'] and counterpart['lens_catalog_arcsec']==old['lens_catalog_arcsec']
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

    required=target/dr
    y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared'])
    full_fac=cho_factor(cov);fac=cho_factor(cov[:-1,:-1]);cw=cho_solve(fac,cov[:-1,-1]);csd=np.sqrt(cov[-1,-1]-cov[-1,:-1]@cw)
    constraint_x=np.geomspace(1e-6,100,4097);check_x=np.geomspace(1e-6,100,8193)
    gamma=deproject_slope(profiles[name],constraint_x,512)[1];fine_gamma=deproject_slope(profiles[name],check_x,512)[1]
    nmax=max(c['n'] for c in profiles[name]['components'] if c['amp_at_R']>0);assert nmax>1
    central_gamma=1-1/nmax
    @lru_cache(maxsize=128)
    def unit_profile(log_rs):
        ac_trial=Re*np.exp(log_rs)
        xx=r/ac_trial;tt=xx[:,None]*mu;bb2=1+xx[:,None]**2*(1-mu*mu);bb=np.sqrt(bb2)
        optical=fit['k0_per_kpc']*ac_trial*(tt/(2*bb2*(bb2+tt*tt))+(np.arctan(tt/bb)+np.pi/2)/(2*bb**3))
        intensity=np.exp(-optical)@w/2
        density=intensity/(1+xx*xx)**2
        mass_unit=4*np.pi*(density[0]*r[0]**3/3+cumulative_trapezoid(r*r*density,r,initial=0))
        interpolator=PchipInterpolator(np.log(r),mass_unit)
        def mass_at(rr):
            if rr<r[0]:return mass_unit[0]*(rr/r[0])**3
            if rr>r[-1]:return mass_unit[-1]
            return float(interpolator(np.log(rr)))
        integral=quad(lambda t:mass_at(impact/np.cos(t))/(impact/np.cos(t)),0,np.pi/2,epsabs=1e-8,epsrel=1e-9,limit=250)[0]
        unit=4*G/C**2*integral
        assert unit>0 and np.all(mass_unit>0)
        return mass_unit,unit
    def unit_deflection(log_rs):return unit_profile(float(log_rs))[1]
    def cap(bi,log_ra):
        ff=constraint_x**2/(np.exp(2*log_ra)+constraint_x**2)
        return min(.45,central_gamma/2,float(np.min((gamma/2-bi*ff)/(1-ff))))
    def physical(q):
        h,t,bi,log_ra,log_rs,f=q
        return np.array([h,-2+t*(cap(bi,log_ra)+2),bi,log_ra,log_rs,f])
    def mapped(q):
        h,b0,bi,log_ra,log_rs,f=q
        return np.array([h,np.clip((b0+2)/(cap(bi,log_ra)+2),0,1),bi,log_ra,log_rs,f])
    def evaluate(q,direct=False):
        h,b0,bi,log_ra,log_rs,f=q;ra=Re*np.exp(log_ra);rs=Re*np.exp(log_rs)
        weight=h*meanH/(1+h*meanH)
        mass=1e11*required*(1-f)/((1-weight)*s0+weight*sh)
        if direct:mass_unit,unit=unit_profile(float(log_rs))
        else:mass_unit=np.exp(mass_table(log_rs));unit=float(np.exp(deflection_table(log_rs)))
        amplitude=required*f/unit
        force=mass/1e11*((1-weight)*model.forces[0]+weight*model.forces[1])+G*amplitude*mass_unit/r**2
        beta=b0+(bi-b0)*r*r/(r*r+ra*ra)
        factor=np.exp(2*b0*np.log(r/model.a)+(bi-b0)*np.log((r*r+ra*ra)/(model.a*model.a+ra*ra)))
        pressure=-cumulative_trapezoid((nu*force*factor)[::-1],r[::-1],initial=0)[::-1]/factor
        v2=np.trapezoid(r*r*pressure*(model.W-beta*model.T),r,axis=1)/model.den
        assert np.all(v2>0) and 1e7<mass<1e14
        pred=np.sqrt(v2);err=y-pred
        return float(err@cho_solve(full_fac,err)),pred,mass,amplitude

    row=next(v for v in fit_results['rows'] if v['Name']==name)
    q=[row['h'],row['beta0'],row['beta_infinity'],np.log(row['orbit_transition_radius_Re']),np.log(row['capture_scale_Re']),row['halo_bending_fraction']]
    mu,w=np.polynomial.legendre.leggauss(96);unit_profile.cache_clear()
    initial=evaluate(q,True)
    assert abs(initial[0]-row['all_motion_chi2'])<1e-5
    checks=[]
    for order in [192,384]:
        mu,w=np.polynomial.legendre.leggauss(order);unit_profile.cache_clear()
        current=evaluate(q,True)
        checks.append(dict(angular_order=order,all_motion_chi2=current[0],max_velocity_difference_from_96_kms=float(max(abs(current[1]-initial[1]))),density_normalization_ratio_to_96=float(current[3]/initial[3])))
    masses,unit=unit_profile(float(q[4]));mass_i=PchipInterpolator(np.log(r),masses)
    quad_checks=[]
    for order in [512,1024]:
        zz,ww=np.polynomial.legendre.leggauss(order);theta=(zz+1)*np.pi/4;rr=impact/np.cos(theta)
        mass_values=mass_i(np.log(np.clip(rr,r[0],r[-1])))
        mass_values=np.where(rr<r[0],masses[0]*(rr/r[0])**3,np.where(rr>r[-1],masses[-1],mass_values))
        fixed_unit=4*G/C**2*np.sum(ww*mass_values/rr)*np.pi/4
        quad_checks.append(dict(order=order,unit_deflection_relative_difference=float(fixed_unit/unit-1)))
    ac_fit=Re*np.exp(q[4])
    breaks=sorted(set([0.,np.pi/2]+[float(np.arccos(impact/rr)) for rr in ac_fit*np.array([.01,.1,1,10,100]) if rr>impact]))
    split_checks=[]
    for order in [256,512]:
        zz,ww=np.polynomial.legendre.leggauss(order);integral=0.
        for lo,hi in zip(breaks[:-1],breaks[1:]):
            theta=lo+(zz+1)*(hi-lo)/2;rr=impact/np.cos(theta)
            mass_values=mass_i(np.log(np.clip(rr,r[0],r[-1])))
            mass_values=np.where(rr<r[0],masses[0]*(rr/r[0])**3,np.where(rr>r[-1],masses[-1],mass_values))
            integral+=np.sum(ww*mass_values/rr)*(hi-lo)/2
        split_checks.append(dict(order_per_interval=order,unit_deflection_relative_difference=float(4*G/C**2*integral/unit-1)))
    out['rows'].append(dict(Name=name,split_deflection_quadrature=split_checks,angular_refinement=checks,deflection_fixed_quadrature=quad_checks,base_score_reproduction_difference=abs(initial[0]-row['all_motion_chi2'])))
    print(name,checks,quad_checks,flush=True)
out['scope']='No refit: fixed fitted parameters, 96/192/384 angle orders and independent fixed deflection quadrature; inherited radial and stellar projection grid unchanged'
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in [P/'free-companion.py',P/'free-companion-results.json']}
(P/'free-companion-refinement-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
