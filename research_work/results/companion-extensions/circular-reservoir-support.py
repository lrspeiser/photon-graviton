"""Frozen circular-orbit particle support diagnostic."""
from pathlib import Path
import sys,json,hashlib,argparse
parser=argparse.ArgumentParser();parser.add_argument("--refine",action="store_true");args=parser.parse_args()
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
scan_results=json.loads((P/'reservoir-scale-scan-results.json').read_text())
selected=[dict(v,label=v['Name']+' free') for v in fit_results['rows']]
extended=next(v for v in fit_results['rows'] if v['Name']=='J1621+3931')
for v in scan_results['rows']:
    if v['density_Msun_kpc3']>0:
        selected.append(dict(v,Name='J1621+3931',label='J1621 scale '+str(v['capture_scale_Re']),stored_density_normalization_Msun_kpc3=v['density_Msun_kpc3'],outer_grid_radius_in_capture_scales=extended['capture_scale_kpc']*extended['outer_grid_radius_in_capture_scales']/v['capture_scale_kpc']))
for fitted in selected:
    old=next(v for v in source['rows'] if v['Name']==fitted['Name'] and v['retention_mapping']['population']=='Chabrier')
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


    ac=fitted['capture_scale_kpc'];D=fitted['stored_density_normalization_Msun_kpc3'];Mstar=fitted['stellar_mass_Msun'];h=fitted['h']
    Rmax=ac*fitted['outer_grid_radius_in_capture_scales']
    mu,w=np.polynomial.legendre.leggauss(384 if args.refine else 192)
    def density_comp(rr):
        rr=np.asarray(rr);shape=rr.shape;xx=rr.reshape(-1)/ac
        tt=xx[:,None]*mu;bb2=1+xx[:,None]**2*(1-mu*mu);bb=np.sqrt(bb2)
        optical=fit['k0_per_kpc']*ac*(tt/(2*bb2*(bb2+tt*tt))+(np.arctan(tt/bb)+np.pi/2)/(2*bb**3))
        jj=np.exp(-optical)@w/2
        return (D*jj/(1+xx*xx)**2).reshape(shape)
    mass_r=np.geomspace(min(r[0],ac*1e-6),Rmax,16001 if args.refine else 8001)
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

    total_mass=mc+star_enclosed(mass_r)
    density_total=dd+density_star(mass_r)
    vc2=G*total_mass/mass_r
    omega2=vc2/mass_r**2
    radial2=G*(total_mass+4*np.pi*mass_r**3*density_total)/mass_r**3
    assert np.all(radial2>0) and np.all(vc2>0)
    dm_dr=4*np.pi*mass_r**2*dd
    integrated_mass=np.trapezoid(dm_dr,mass_r)
    K=.5*np.trapezoid(vc2*dm_dr,mass_r)
    Labs=np.trapezoid(mass_r*np.sqrt(vc2)*dm_dr,mass_r)
    peak=int(np.argmax(dd));rp=float(mass_r[peak]);vpeak=float(np.sqrt(vc2[peak]))
    # Finite-difference angular-momentum derivative checks the epicyclic formula.
    radial_fd=np.gradient(G*total_mass*mass_r,mass_r,edge_order=2)/mass_r**3
    use=(mass_r>max(r[0]*100,ac*1e-3))&(mass_r<Rmax/10)
    fd_error=float(max(abs(radial_fd[use]/radial2[use]-1)))
    assert fd_error<.01
    row=dict(label=fitted['label'],capture_scale_kpc=ac,stored_mass_Msun=float(mc[-1]),kinetic_energy_J=float(K*1.98847e36),kinetic_to_rest_energy=float(K/(mc[-1]*C*C)),mass_weighted_absolute_specific_angular_momentum_kpc_kms=float(Labs/integrated_mass),net_vector_angular_momentum='zero by paired orientations; not zero individual angular momentum',max_circular_speed_kms=float(np.sqrt(max(vc2))),max_tangential_pressure_over_rho_c2=float(max(vc2)/(2*C*C)),minimum_radial_frequency_squared_over_omega_squared=float(min(radial2/omega2)),finite_difference_radial_frequency_max_relative_error=fd_error,density_peak_radius_kpc=rp,circular_speed_at_density_peak_kms=vpeak,orbital_period_at_density_peak_years=float(2*np.pi*rp*3.0856775814913673e16/vpeak/(365.25*86400)),isotropic_incident_low_angular_momentum_fraction_at_peak=float(vc2[peak]/C**2),mass_quadrature_relative_error=float(abs(integrated_mass/mc[-1]-1)),radial_points=len(mass_r),incoming_angles=len(mu))
    mass_tail=4*np.pi*D*ac**4/Rmax
    mass_upper=Mstar+mc[-1]+mass_tail
    row['kinetic_exterior_tail_bound_J']=float(np.pi*G*mass_upper*D*ac**4/Rmax**2*1.98847e36)
    row['absolute_angular_momentum_exterior_tail_bound_Msun_kpc_kms']=float(8*np.pi*D*ac**4*np.sqrt(G*mass_upper/Rmax))
    row['absolute_angular_momentum_inside_grid_Msun_kpc_kms']=float(Labs)
    row['mass_grid_outer_radius_kpc']=float(Rmax)
    out['rows'].append(row);print(row['label'],row['kinetic_to_rest_energy'],row['kinetic_exterior_tail_bound_J']/row['kinetic_energy_J'],flush=True)
out['scope']='Singular circular-orbit completion and kinetic inventory; individual frozen-potential radial stability only; no formation or collective stability proof'
out['input_sha256']={f.as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),P/'circular-reservoir-support-protocol.md',P/'free-companion-results.json',P/'reservoir-scale-scan-results.json']}
out['inherited_inputs']=json.loads((P/'outer-companion-predictions-results.json').read_text())['input_sha256']
(P/('circular-reservoir-support-refined.json' if args.refine else 'circular-reservoir-support-results.json')).write_text(json.dumps(out,indent=2)+chr(10),encoding='utf-8',newline=chr(10))
