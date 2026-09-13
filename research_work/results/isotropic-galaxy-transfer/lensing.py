"""Frozen capture source, inner-star nuisance fits, outer/lens predictions."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid,quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize,brentq
from scipy.linalg import cho_factor,cho_solve
HERE=Path(__file__).resolve().parent
RELAXING='--relaxing-optics' in sys.argv
RETENTION='--retention-optics' in sys.argv
REGULAR='--regular-optics' in sys.argv or RELAXING or RETENTION
OPTICAL_FILE='relaxing-area-results.json' if RELAXING else 'regular-area-results.json'
sys.path.insert(0,str(HERE.parent/'slacs-component-refit'))
from model import ComponentModel,G,C,ARCSEC
def read(folder,file='results.json'):return json.loads((HERE.parent/folder/file).read_text())
data=read('slacs-resolved-input-audit');base=read('slacs-motion-lensing-pilot')
geo={r['Name']:r for r in read('lensing-data-readiness','conditional-geometry.json')}
profiles={r['Name']:r for r in read('slacs-light-profile-audit')['rows']}
pilot={r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
cfg=read('slacs-outer-bin-check','protocol.json')
allowed={r['Name'] for r in read('slacs-outer-bin-check')['rows'] if r['model']=='empirical_extra'}
capture=json.loads((HERE/'results.json').read_text())['models']
capture_file='results.json'
if RETENTION:
    capture_file='bounded-radiation-retention-results.json'
    retained=json.loads((HERE/capture_file).read_text())['models']['attenuated']
    capture={f'attenuated_{imf}':dict(retained,population=imf) for imf in ['Chabrier','Salpeter']}
    photorows=read('lens-photometric-audit','normalization-sensitivity.json')
    photo={(r['Name'],r['imf']):r for r in photorows if r['propagation_branch']=='energy_loss_and_event_stretch' and r['model']=='baryons'}
optical=read('brightness-distance-consistency',OPTICAL_FILE) if REGULAR else None
observations={r['Name']:r for r in read('lensing-data-readiness','lens-observations-and-image-models.json')} if REGULAR else None
mu,w=np.polynomial.legendre.leggauss(96)
rows=[]
for item in data['systems']:
    name=item['Name']
    if name not in allowed:continue
    dl=geo[name]['conditional_Dl_Mpc']*1000;ratio=geo[name]['conditional_Dls_over_Ds'];a=pilot[name]['scale_a_kpc']
    geometry_info=None
    if REGULAR:
        alpha=optical['alpha_per_mpc'];q=optical['q'];zl=observations[name]['zFG'];zs=observations[name]['zBG']
        fl=zl/(1+zl);fs=zs/(1+zs)
        def H(f):
            area=1+f*np.exp(-q*f) if RELAXING else 1+f/(1+q*f)
            return -(1-f)*np.log1p(-f)*np.sqrt(area)
        angular_l=H(fl)/alpha
        newratio=(1+zl)*H(fl)*quad(lambda f:1/H(f)**2,fl,fs,epsabs=1e-9,epsrel=1e-10)[0]
        geometry_info=dict(z_lens=zl,z_source=zs,path_Dl_Mpc=np.log1p(zl)/alpha,angular_Dl_Mpc=angular_l,angular_Ds_Mpc=H(fs)/alpha,Dls_over_Ds=newratio,previous_Dls_over_Ds=ratio)
        a*=angular_l*1000/dl;dl=angular_l*1000;ratio=newratio
    Re=profiles[name]['computed_equal_area_half_light_arcsec']*dl/ARCSEC
    equiv=Re/1.67834699
    edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*dl/ARCSEC
    psf=item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
    components=[dict(R=q['R_arcsec']*dl/ARCSEC,n=q['n'],amp=q['amp_at_R'],bn=q['bn']) for q in profiles[name]['components']]
    model=ComponentModel(a,edges,psf,0,.5,1,20,components)
    r=model.r;y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);fac=cho_factor(cov[:-1,:-1])
    crossw=cho_solve(fac,cov[:-1,-1]);csd=np.sqrt(cov[-1,-1]-cov[-1,:-1]@crossw)
    for branch,cp in capture.items():
        ac=equiv*cp['scale_to_disk'];x=r/ac;shape=(1+x*x)**-2
        if branch.startswith('attenuated'):
            t=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu**2);B=np.sqrt(B2)
            tau=cp['k0_per_kpc']*ac*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
            J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1)
        else:J=np.ones_like(x)
        rho=cp['C_Msun_kpc3']*shape*J
        retention_info=None
        if RETENTION:
            Barea=1+fl/(1+optical['q']*fl)
            photomass=10**photo[name,cp['population']]['conditional_log10_stellar_mass']*Barea
            Lproxy=photomass/.5
            X=Lproxy/1e9/equiv**2;eta=X**cp['q']/(1+X**cp['q'])
            rho*=2*eta
            retention_info=dict(population=cp['population'],population_mass_Msun=photomass,L3_6_proxy_Lsun=Lproxy,proxy_X=X,eta=eta,luminosity_mapping='Population mass divided by 0.5; not observed rest-frame 3.6-micron luminosity')
        mc=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(rho*r*r,r,initial=0))
        mi=PchipInterpolator(np.log(r),mc,extrapolate=False)
        def Mextra(rr):
            if rr<r[0]:return mc[0]*(rr/r[0])**3
            if rr>r[-1]:return mc[-1]
            return float(mi(np.log(rr)))
        model.forces[1]=G*mc/r**2
        def predict(mass,beta):
            cb,cc=model.coefficients(beta)
            return np.sqrt(cb*(mass/1e11)+cc)
        def objective(q):
            e=y[:-1]-predict(np.exp(q[0]),q[1])[:-1]
            return float(e@cho_solve(fac,e))
        fits=[minimize(objective,[np.log(1e11),beta],method='L-BFGS-B',bounds=[np.log(cfg['mass_Msun_bounds']),cfg['constant_beta_bounds']],options={'ftol':1e-11,'maxiter':1000}) for beta in cfg['starts_beta']]
        good=[f for f in fits if f.success and np.isfinite(f.fun)];assert good
        fit=min(good,key=lambda f:f.fun);mass=float(np.exp(fit.x[0]));beta=float(fit.x[1])
        pred=predict(mass,beta);outer=pred[-1]+crossw@(y[:-1]-pred[:-1])
        def lens_residual(b):
            def integrand(t):
                rr=b/np.cos(t)
                return G*(mass*model.mass_fraction(rr)+Mextra(rr))/rr
            bend=4/C**2*quad(integrand,0,np.pi/2,epsabs=1e-6,epsrel=1e-7,limit=200)[0]
            return ratio*bend-b/dl
        lo=a*1e-7;hi=a*1e5
        assert lens_residual(lo)*lens_residual(hi)<0
        angle=brentq(lens_residual,lo,hi,xtol=1e-9)/dl*ARCSEC
        rows.append(dict(Name=name,model=branch,equivalent_disk_scale_kpc=equiv,capture_scale_kpc=ac,mass_Msun=mass,beta=beta,orbit_boundary=bool(np.min(abs(beta-np.array(cfg['constant_beta_bounds'])))<1e-5),inner_chi2=float(fit.fun),outer_observed_kms=float(y[-1]),outer_prediction_kms=float(outer),outer_conditional_standardized_residual=float((y[-1]-outer)/csd),lens_prediction_arcsec=angle,lens_catalog_arcsec=pilot[name]['catalog_SIE_arcsec'],lens_fractional_residual=angle/pilot[name]['catalog_SIE_arcsec']-1,optimizer_successes=sum(bool(f.success) for f in fits),observed_stellar_vrms=y.tolist(),predicted_stellar_vrms=pred.tolist()))
        if REGULAR:rows[-1]['geometry']=geometry_info
        if RETENTION:rows[-1]['retention_mapping']=retention_info
        print(name,branch,fit.fun,angle,flush=True)
summary=[]
for branch in capture:
    rr=[d for d in rows if d['model']==branch]
    summary.append(dict(model=branch,n=len(rr),inner_chi2_sum=sum(d['inner_chi2'] for d in rr),outer_conditional_residual_square_sum=sum(d['outer_conditional_standardized_residual']**2 for d in rr),lens_fractional_rms=float(np.sqrt(np.mean([d['lens_fractional_residual']**2 for d in rr]))),orbit_boundary_count=sum(d['orbit_boundary'] for d in rr)))
out=dict(summary=summary,rows=rows,capture_input_sha256=hashlib.sha256((HERE/capture_file).read_bytes()).hexdigest())
if REGULAR:out['optical_input_sha256']=hashlib.sha256((HERE.parent/'brightness-distance-consistency'/OPTICAL_FILE).read_bytes()).hexdigest()
output='relaxing-optics-results.json' if RELAXING else 'regular-optics-results.json' if REGULAR else 'lensing-results.json'
if RETENTION:
    output='retention-optics-results.json'
    out['photometric_input_sha256']=hashlib.sha256((HERE.parent/'lens-photometric-audit/normalization-sensitivity.json').read_bytes()).hexdigest()
(HERE/output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
