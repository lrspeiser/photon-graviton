"""Free companion capture scale and density strength, matched to free NFW freedoms."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize,brentq
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
fit_results=json.loads((P/"free-companion-results.json").read_text())
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


    fitted=next(v for v in fit_results['rows'] if v['Name']==name)
    ac=fitted['capture_scale_kpc'];D=fitted['stored_density_normalization_Msun_kpc3'];Mstar=fitted['stellar_mass_Msun'];h=fitted['h']
    Rmax=ac*fitted['outer_grid_radius_in_capture_scales']
    mu,w=np.polynomial.legendre.leggauss(192)
    def density_comp(rr):
        rr=np.asarray(rr);shape=rr.shape;xx=rr.reshape(-1)/ac
        tt=xx[:,None]*mu;bb2=1+xx[:,None]**2*(1-mu*mu);bb=np.sqrt(bb2)
        optical=fit['k0_per_kpc']*ac*(tt/(2*bb2*(bb2+tt*tt))+(np.arctan(tt/bb)+np.pi/2)/(2*bb**3))
        jj=np.exp(-optical)@w/2
        return (D*jj/(1+xx*xx)**2).reshape(shape)
    mass_r=np.geomspace(min(r[0],ac*1e-6),Rmax,8001)
    dd=density_comp(mass_r)
    mc=4*np.pi*(dd[0]*mass_r[0]**3/3+cumulative_trapezoid(mass_r**2*dd,mass_r,initial=0))
    mass_interp=PchipInterpolator(np.log(mass_r),mc)
    def comp_enclosed(rr):
        rr=np.asarray(rr)
        val=mass_interp(np.log(np.clip(rr,mass_r[0],mass_r[-1])))
        return np.where(rr<mass_r[0],mc[0]*(rr/mass_r[0])**3,np.where(rr>Rmax,mc[-1],val))
    weight=h*meanH/(1+h*meanH)
    star_frac=(1-weight)*f0+weight*core
    star_mass_interp=PchipInterpolator(np.log(r),star_frac)
    nu_mass=Mstar*nu*(1+h*H)/(lum[-1]*(1+h*meanH))
    log_nu=PchipInterpolator(np.log(r),np.log(np.maximum(nu_mass,1e-300)))
    def star_enclosed(rr):
        rr=np.asarray(rr);val=star_mass_interp(np.log(np.clip(rr,r[0],r[-1])))
        return Mstar*np.where(rr<r[0],star_frac[0]*(rr/r[0])**model.inner_power,np.where(rr>r[-1],1.,val))
    def density_star(rr):
        rr=np.asarray(rr)
        return np.where(rr>r[-1],0,np.exp(log_nu(np.log(np.clip(rr,r[0],r[-1])))))
    sigma_crit=C*C/(4*np.pi*G*dl*dr)
    def projection(b,order=192):
        b=np.atleast_1d(b);zz,ww=np.polynomial.legendre.leggauss(order)
        breaks=np.column_stack([np.zeros_like(b)]+[np.arccos(np.clip(b/(ac*v),0,1)) for v in [.01,.1,1,10,100]]+[np.full_like(b,np.pi/2)])
        mp_comp=np.zeros_like(b);mp_star=np.zeros_like(b)
        for j in range(breaks.shape[1]-1):
            lo=breaks[:,j];span=breaks[:,j+1]-lo
            theta=lo[:,None]+(zz+1)*span[:,None]/2;tr=b[:,None]/np.cos(theta)
            mp_comp+=np.sum(comp_enclosed(tr)*np.cos(theta)*ww,axis=1)*span/2
            mp_star+=np.sum(star_enclosed(tr)*np.cos(theta)*ww,axis=1)*span/2
        umax=np.arccosh(Rmax/b);u=umax[:,None]*(zz+1)/2;rr=b[:,None]*np.cosh(u)
        sigma_comp=np.sum(density_comp(rr)*np.cosh(u)*ww,axis=1)*b*umax
        sigma_star=np.sum(density_star(rr)*np.cosh(u)*ww,axis=1)*b*umax
        bar=(mp_comp+mp_star)/(np.pi*b*b*sigma_crit)
        kap=(sigma_comp+sigma_star)/sigma_crit;shear=bar-kap
        return dict(kappa=kap,mean_kappa=bar,tangential_shear=shear,reduced_shear=shear/(1-kap),lambda_t=1-bar,lambda_r=1-2*kap+bar,companion_kappa=sigma_comp/sigma_crit,companion_mean_kappa=mp_comp/(np.pi*b*b*sigma_crit),stellar_kappa=sigma_star/sigma_crit)
    scan=np.geomspace(.01,10000,801)
    collected=[]
    for start in range(0,len(scan),24):collected.append(projection(scan[start:start+24]))
    curves={key:np.concatenate([v[key] for v in collected]) for key in collected[0]}
    roots={}
    for key in ['lambda_t','lambda_r']:
        ii=np.where(curves[key][:-1]*curves[key][1:]<0)[0]
        roots[key]=[float(brentq(lambda x:projection([x],384)[key][0],scan[i],scan[i+1],xtol=1e-7)) for i in ii]
    root_checks={key:[float(abs(projection([value],768)[key][0])) for value in values] for key,values in roots.items()}
    assert all(v<1e-4 for values in root_checks.values() for v in values), (name,root_checks)
    potential_comp=G*(mc[0]/mass_r[0]+4*np.pi*np.trapezoid(dd*mass_r,mass_r))/C**2
    potential_star=G*(Mstar*star_frac[0]/r[0]+4*np.pi*np.trapezoid(nu_mass*r,r))/C**2
    potential_tail_bound=2*np.pi*G*D*ac**4/(Rmax**2*C**2)
    radii=np.array([10.,30.,100.,300.,1000.,3000.,10000.]);pr=projection(radii,384);coarse=projection(radii,192)
    cat=projection([impact],384);lens_error=abs(cat['mean_kappa'][0]-1)
    assert lens_error<5e-4
    speed=np.sqrt(G*(comp_enclosed(radii)+star_enclosed(radii))/radii)
    vscan=np.sqrt(G*(comp_enclosed(scan)+star_enclosed(scan))/scan);peak=int(np.argmax(vscan))
    # M_2d derivative identity is checked independently with symmetric perturbations.
    eps=1e-3
    plus=projection(radii*(1+eps),384)['mean_kappa']*(radii*(1+eps))**2
    minus=projection(radii*(1-eps),384)['mean_kappa']*(radii*(1-eps))**2
    deriv=(plus-minus)/(4*eps*radii*radii)
    consistency=float(max(abs(deriv-pr['kappa'])))
    out['rows'].append(dict(Name=name,capture_scale_kpc=ac,sigma_critical_Msun_kpc2=sigma_crit,source_geometry=old['geometry'],catalogue_lens_radius_kpc=impact,catalogue_lens_mean_kappa_error=float(lens_error),mass_grid_upper_radius_kpc=Rmax,mass_inside_grid_Msun=float(mc[-1]),relative_mass_difference_from_fit_grid=float(mc[-1]/fitted['stored_mass_within_grid_Msun']-1),outer_predictions=[dict(radius_kpc=float(b),angle_arcsec=float(b/dl*ARCSEC),circular_speed_kms=float(v),stellar_only_circular_speed_kms=float(np.sqrt(G*star_enclosed(b)/b)),companion_enclosed_Msun=float(comp_enclosed(b)),**{k:float(vv[i]) for k,vv in pr.items()}) for i,(b,v) in enumerate(zip(radii,speed))],critical_radii_kpc=roots,refined_critical_curve_eigenvalue_residuals=root_checks,central_potential_magnitude_over_c2_within_grid=float(potential_comp+potential_star),exterior_companion_potential_bound_over_c2=float(potential_tail_bound),peak_scanned_circular_speed_kms=float(vscan[peak]),peak_speed_radius_kpc=float(scan[peak]),peak_scanned_compactness_2GM_rc2=float(2*max(vscan*vscan)/C**2),maximum_scanned_companion_kappa=float(max(curves['companion_kappa'])),surface_projection_derivative_max_kappa_difference=consistency,quadrature_max_absolute_kappa_difference=float(max(abs(pr['kappa']-coarse['kappa']))),quadrature_max_absolute_mean_kappa_difference=float(max(abs(pr['mean_kappa']-coarse['mean_kappa'])))))
    print(name,'peak vc',vscan[peak],'roots',roots,'lens error',lens_error,'projection check',consistency,flush=True)
out['scope']='Frozen extrapolations, not observed outer data or refits; spherical isolated weak-field thin-lens model at the adopted source geometry'
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in [P/'free-companion-results.json',P/'outer-companion-predictions-protocol.md',L/'capacity-reference-optics-results.json',P.parent/'slacs-light-profile-audit/results.json',P.parent/'slacs-component-refit/model.py',P.parent/'slacs-resolved-fit/model.py',P.parent/'slacs-resolved-input-audit/results.json',P.parent/'slacs-motion-lensing-pilot/results.json',P.parent/'lensing-data-readiness/conditional-geometry.json',L/'third-radiation-retention-results.json']}
(P/'outer-companion-predictions-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
