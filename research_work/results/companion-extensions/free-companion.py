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
    rs_grid=np.linspace(np.log(.1),np.log(100),257)
    grid_profiles=[unit_profile(v) for v in rs_grid]
    mass_table=PchipInterpolator(rs_grid,np.log([v[0] for v in grid_profiles]),axis=0)
    deflection_table=PchipInterpolator(rs_grid,np.log([v[1] for v in grid_profiles]))
    del grid_profiles
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
    seeds=[];reproduction=[]
    for row in reference['rows']:
        if row['Name']!=name:continue
        old_case=next(v for v in source['rows'] if v['Name']==name and v['retention_mapping']['population']==row['population'])
        log_rs=np.log(old_case['capture_scale_kpc']/Re)
        ref_density=fit['source_C_before_retention_Msun_kpc3']*source['C0_multiplier']*old_case['retention_mapping']['eta']
        f=ref_density*unit_deflection(log_rs)/required
        q=[row['h'],row['beta0'],row['beta_infinity'],np.log(row['orbit_transition_radius_Re']),log_rs,f]
        err=abs(evaluate(q,True)[0]-row['all_motion_chi2']);assert err<1e-4
        reproduction.append(err)
        seeds.append(q)
    st=next(v for v in stars['rows'] if v['Name']==name and v['population']=='Chabrier')
    star_seed=[st['h'],st['beta0'],st['beta_infinity'],np.log(st['orbit_transition_radius_Re']),np.log(3),0.]
    star_drift=abs(evaluate(star_seed,True)[0]-st['all_motion_chi2']);assert star_drift<1e-5
    seeds.append(star_seed)
    nf=next(v for v in nfw['rows'] if v['Name']==name)
    seeds.append([nf['h'],nf['beta0'],nf['beta_infinity'],np.log(nf['orbit_transition_radius_Re']),np.log(nf['nfw_rs_Re']),nf['halo_bending_fraction']])
    for rs_start in [.3,3,30]:
        for f in [.1,.4,.8]:
            for h,b0,bi in [(0,0,.5),(2,-.5,-.5)]:seeds.append([h,b0,bi,0,np.log(rs_start),f])
    bounds=[(-.8,9),(0,1),(-2,.95),(-np.log(10),np.log(10)),(np.log(.1),np.log(100)),(0,.95)]
    trials=[minimize(lambda q:evaluate(physical(q))[0],mapped(q),method='L-BFGS-B',bounds=bounds,options={'ftol':1e-10,'maxiter':800}) for q in seeds]
    good=sorted([v for v in trials if v.success],key=lambda v:v.fun);assert good
    polished=[minimize(lambda q:evaluate(physical(q),True)[0],v.x,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-11,'maxiter':800}) for v in good[:3]]
    successful=[v for v in polished if v.success];assert successful
    best=min(successful,key=lambda v:v.fun);q=physical(best.x);score,pred,mass,amplitude=evaluate(q,True)
    h,b0,bi,log_ra,log_rs,f=q
    seed_scores=[evaluate(v,True)[0] for v in seeds]
    assert score<=min(seed_scores)+1e-4
    inner_err=y[:-1]-pred[:-1];inner=float(inner_err@cho_solve(fac,inner_err))
    outer=(y[-1]-pred[-1]-cw@inner_err)/csd
    decomposition=abs(score-inner-outer**2);assert decomposition<1e-5
    beta_fine=b0+(bi-b0)*check_x**2/(np.exp(2*log_ra)+check_x**2)
    margin=fine_gamma-2*beta_fine;assert min(margin)>-1e-6 and central_gamma-2*b0>=-1e-10
    interp_difference=abs(float(np.exp(deflection_table(log_rs)))/unit_deflection(log_rs)-1)
    lens_error=abs((mass/1e11*((1-h*meanH/(1+h*meanH))*s0+(h*meanH/(1+h*meanH))*sh)+amplitude*unit_deflection(log_rs))/required-1)
    assert lens_error<1e-10
    out['rows'].append(dict(Name=name,motion_bins=len(y),h=float(h),beta0=float(b0),beta_infinity=float(bi),orbit_transition_radius_Re=float(np.exp(log_ra)),capture_scale_Re=float(np.exp(log_rs)),capture_scale_kpc=float(Re*np.exp(log_rs)),halo_bending_fraction=float(f),stellar_mass_Msun=mass,stored_density_normalization_Msun_kpc3=amplitude,all_motion_chi2=score,inner_chi2=inner,conditional_outer_residual=float(outer),observed_stellar_vrms=y.tolist(),predicted_stellar_vrms=pred.tolist(),stellar_only_chi2=st['all_motion_chi2'],companion_chi2={v['population']:v['all_motion_chi2'] for v in reference['rows'] if v['Name']==name},boundaries={k:bool(min(abs(v-lo),abs(v-hi))<1e-5) for k,v,(lo,hi) in zip(['h','central_slope_coordinate','beta_infinity','log_orbit_radius','log_capture_radius','halo_fraction'],best.x,bounds)},minimum_refined_slope_margin=float(min(margin)),central_slope_margin=float(central_gamma-2*b0),reference_companion_reproduction_max_difference=max(reproduction),zero_halo_reproduction_difference=star_drift,covariance_decomposition_difference=decomposition,lens_fractional_residual=lens_error,deflection_interpolation_fractional_difference=interp_difference,optimizer_attempts=len(trials),optimizer_successes=len(good),successful_objectives=[float(v.fun) for v in good],direct_polish_successes=len(successful),direct_polish_objectives=[float(v.fun) for v in successful],failed_optimizer_messages=[str(v.message) for v in trials if not v.success]))
    out['rows'][-1].update(free_nfw_chi2=nf['all_motion_chi2'],reference_parameter_ratios={v['retention_mapping']['population']:dict(capture_scale_ratio=float(Re*np.exp(log_rs)/v['capture_scale_kpc']),density_normalization_ratio=float(amplitude/(fit['source_C_before_retention_Msun_kpc3']*source['C0_multiplier']*v['retention_mapping']['eta']))) for v in source['rows'] if v['Name']==name},outer_grid_radius_in_capture_scales=float(r[-1]/(Re*np.exp(log_rs))),stored_mass_within_grid_Msun=float(amplitude*unit_profile(float(log_rs))[0][-1]),force_interpolation_max_relative_difference=float(np.max(np.abs(np.exp(mass_table(log_rs))/unit_profile(float(log_rs))[0]-1))))
    print(name,q,score,outer,flush=True)
paths=[L/'capacity-reference-optics-results.json',L/'capacity-reference-exact-lens-optics-results.json',L/'third-radiation-retention-results.json',P.parent/'slacs-component-refit/model.py',P.parent/'slacs-resolved-fit/model.py',P/'orbit_density.py',P/'orbit-transition-results.json',P/'stellar-only-lens-control-results.json',P/'free-companion-protocol.md',P/'free-nfw-results.json']
paths += [P.parent/f/n for f,n in [('slacs-resolved-input-audit','results.json'),('slacs-light-profile-audit','results.json'),('slacs-motion-lensing-pilot','results.json'),('lensing-data-readiness','conditional-geometry.json'),('slacs-outer-bin-check','protocol.json')]]
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
assert len(out['rows'])==6
out['summary']=dict(total_motion_bins=sum(v['motion_bins'] for v in out['rows']),local_parameters_per_galaxy=6,all_motion_chi2=sum(v['all_motion_chi2'] for v in out['rows']),inner_chi2=sum(v['inner_chi2'] for v in out['rows']),outer_residual_square_sum=sum(v['conditional_outer_residual']**2 for v in out['rows']))
(P/'free-companion-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2))
