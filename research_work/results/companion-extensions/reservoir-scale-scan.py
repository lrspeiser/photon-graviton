"""Fixed-scale J1621 diagnostic; inherited setup from free-companion.py."""
from pathlib import Path
import sys,json,hashlib,argparse
parser=argparse.ArgumentParser();parser.add_argument("--angles",type=int,default=96);args=parser.parse_args()
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
cache={};mu,w=np.polynomial.legendre.leggauss(args.angles)
for old in source['rows']:
    if old['retention_mapping']['population']!='Chabrier' or old['Name']!='J1621+3931':continue
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
    def cap(bi,log_ra):
        ff=constraint_x**2/(np.exp(2*log_ra)+constraint_x**2)
        return min(.45,central_gamma/2,float(np.min((gamma/2-bi*ff)/(1-ff))))
    def physical(q):
        h,t,bi,log_ra,log_rs,f=q
        return np.array([h,-2+t*(cap(bi,log_ra)+2),bi,log_ra,log_rs,f])
    def mapped(q):
        h,b0,bi,log_ra,log_rs,f=q
        return np.array([h,np.clip((b0+2)/(cap(bi,log_ra)+2),0,1),bi,log_ra,log_rs,f])
    def evaluate(q,direct=True):
        h,b0,bi,log_ra,log_rs,f=q;ra=Re*np.exp(log_ra);rs=Re*np.exp(log_rs)
        weight=h*meanH/(1+h*meanH)
        mass=1e11*required*(1-f)/((1-weight)*s0+weight*sh)
        mass_unit,unit=unit_profile(float(log_rs))
        amplitude=required*f/unit
        force=mass/1e11*((1-weight)*model.forces[0]+weight*model.forces[1])+G*amplitude*mass_unit/r**2
        beta=b0+(bi-b0)*r*r/(r*r+ra*ra)
        factor=np.exp(2*b0*np.log(r/model.a)+(bi-b0)*np.log((r*r+ra*ra)/(model.a*model.a+ra*ra)))
        pressure=-cumulative_trapezoid((nu*force*factor)[::-1],r[::-1],initial=0)[::-1]/factor
        v2=np.trapezoid(r*r*pressure*(model.W-beta*model.T),r,axis=1)/model.den
        assert np.all(v2>0) and 1e7<mass<1e14
        pred=np.sqrt(v2);err=y-pred
        return float(err@cho_solve(full_fac,err)),pred,mass,amplitude

    frozen=json.loads((P/'free-companion-results.json').read_text())
    fr=next(v for v in frozen['rows'] if v['Name']==name)
    st=next(v for v in stars['rows'] if v['Name']==name and v['population']=='Chabrier')
    out['scope']='Fixed-scale conditional profile diagnostic, five locally fitted parameters; all motion bins and lens calibration consumed'
    out['stellar_only_chi2']=st['all_motion_chi2'];out['free_companion_chi2']=fr['all_motion_chi2']
    previous=None
    for scale in [.1,.3,1,3,10,30,100]:
        logac=float(np.log(scale))
        def expand(z):return np.array([*z[:4],logac,z[4]])
        def objective(z):return evaluate(physical(expand(z)))[0]
        bounds=[(-.8,9),(0,1),(-2,.95),(-np.log(10),np.log(10)),(0,.95)]
        seeds=[]
        for row,fraction in [(st,0),(fr,fr['halo_bending_fraction'])]:
            q=mapped([row['h'],row['beta0'],row['beta_infinity'],np.log(row['orbit_transition_radius_Re']),logac,fraction]);seeds.append(q[[0,1,2,3,5]])
        zero=objective(seeds[0]);assert abs(zero-st['all_motion_chi2'])<1e-5
        for fraction in [.02,.15,.5]:
            for h,b0,bi,ra in [(0,0,.5,1),(2,-.5,-.5,1),(-.6,.375,.95,.3),(-.6,.375,.95,3)]:
                q=mapped([h,b0,bi,np.log(ra),logac,fraction]);seeds.append(q[[0,1,2,3,5]])
        if previous is not None:seeds.append(previous)
        trials=[minimize(objective,z,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-12,'maxiter':1200,'maxls':40}) for z in seeds]
        good=[t for t in trials if t.success];assert good
        best=min(good,key=lambda t:t.fun);previous=best.x
        q=physical(expand(best.x));score,pred,ms,D=evaluate(q);h,b0,bi,logra,_,fraction=q
        assert score<=zero+1e-5
        innererr=y[:-1]-pred[:-1];inner=float(innererr@cho_solve(fac,innererr));outer=float((y[-1]-pred[-1]-cw@innererr)/csd)
        beta=b0+(bi-b0)*check_x**2/(np.exp(2*logra)+check_x**2);margin=float(np.min(fine_gamma-2*beta))
        assert margin>-1e-6 and central_gamma-2*b0>=-1e-10
        weight=h*meanH/(1+h*meanH);lens=abs((ms/1e11*((1-weight)*s0+weight*sh)+D*unit_profile(logac)[1])/required-1)
        assert lens<1e-10 and abs(score-inner-outer**2)<1e-8
        mass=float(D*unit_profile(logac)[0][-1]);tail=float(4*np.pi*D*(Re*scale)**4/r[-1])
        row=dict(capture_scale_Re=scale,capture_scale_kpc=Re*scale,all_motion_chi2=score,improvement_over_stellar_only=zero-score,excess_over_free_companion=score-fr['all_motion_chi2'],h=float(h),beta0=float(b0),beta_infinity=float(bi),orbit_transition_radius_Re=float(np.exp(logra)),bending_fraction=float(fraction),density_Msun_kpc3=D,stellar_mass_Msun=ms,stored_mass_Msun=mass,exterior_mass_bound_Msun=tail,energy_equivalent_J=mass*1.98847e30*299792458**2,observed_vrms=y.tolist(),predicted_vrms=pred.tolist(),inner_chi2=inner,conditional_outer_residual=outer,slope_margin=margin,lens_residual=lens,zero_reproduction=abs(zero-st['all_motion_chi2']),optimizer_attempts=len(trials),optimizer_successes=len(good),trials=[dict(success=bool(t.success),score=float(t.fun),message=str(t.message)) for t in trials])
        out['rows'].append(row);print(scale,score,mass,fraction,flush=True)
    assert abs(out['rows'][-1]['all_motion_chi2']-fr['all_motion_chi2'])<1e-3
paths=[P/'reservoir-scale-scan-protocol.md',P/'free-companion.py',P/'free-companion-results.json',P/'stellar-only-lens-control-results.json',Path(__file__)]
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
out['inherited_input_sha256']=frozen['input_sha256']
out['incoming_angles']=args.angles
(P/('reservoir-scale-scan-results.json' if args.angles==96 else f'reservoir-scale-scan-refined-{args.angles}.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
