"""Decompose bending at the observed lens angle with existing fitted stellar models."""
from pathlib import Path
import sys,hashlib,json
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
P=Path(__file__).resolve().parent;ROOT=P.parents[2];L=P.parent/'isotropic-galaxy-transfer'
sys.path.insert(0,str(P.parent/'slacs-component-refit'))
from model import ComponentModel,G,C,ARCSEC
def read(folder,name='results.json'):return json.loads((P.parent/folder/name).read_text())
data={r['Name']:r for r in read('slacs-resolved-input-audit')['systems']}
profiles={r['Name']:r for r in read('slacs-light-profile-audit')['rows']}
pilot={r['Name']:r for r in read('slacs-motion-lensing-pilot')['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
geo={r['Name']:r for r in read('lensing-data-readiness','conditional-geometry.json')}
fit=json.loads((L/'third-radiation-retention-results.json').read_text())['models']['attenuated']
assert abs(fit['source_C_before_retention_Msun_kpc3']/(2*fit['C_Msun_kpc3'])-1)<1e-12
out=dict(scope='Bending decomposition at fixed inner-star nuisance fits and conditional geometry; no new fitted correction or raw-image analysis',rows=[],summary=[])
cache={}
branches=['original','reference-refined','local','recycling-refined']
for branch in branches:
    source=json.loads((L/('capacity-'+branch+'-optics-results.json')).read_text())
    refined=branch.endswith('refined')
    mu,w=np.polynomial.legendre.leggauss(192 if refined else 96)
    for row in source['rows']:
        name=row['Name'];item=data[name];dl=row['geometry']['angular_Dl_Mpc']*1000
        ratio=row['geometry']['Dls_over_Ds'];ac=row['capture_scale_kpc']
        key=(name,refined)
        if key not in cache:
            olddl=geo[name]['conditional_Dl_Mpc']*1000
            a=pilot[name]['scale_a_kpc']*dl/olddl
            edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*dl/ARCSEC
            psf=item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
            components=[dict(R=q['R_arcsec']*dl/ARCSEC,n=q['n'],amp=q['amp_at_R'],bn=q['bn']) for q in profiles[name]['components']]
            cache[key]=ComponentModel(a,edges,psf,0,.5,1,20,components,**(dict(n=8001,order=256,deproj_order=512) if refined else {}))
        model=cache[key];r=model.r;x=r/ac
        t=x[:,None]*mu;b2=1+x[:,None]**2*(1-mu*mu);b=np.sqrt(b2)
        tau=fit['k0_per_kpc']*ac*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
        J=np.exp(-tau)@w/2
        X=row['retention_mapping']['proxy_X'];e=row['retention_mapping']['eta']
        A=fit['source_C_before_retention_Msun_kpc3']*source['C0_multiplier']
        if branch=='local':
            y=np.cbrt(X*J);rho=A/(1+x*x)**2*y/(1+y)
        else:rho=A*e/(1+x*x)**2*(1 if branch.startswith('recycling') else J)
        mass=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(rho*r*r,r,initial=0));mi=PchipInterpolator(np.log(r),mass)
        def extra(rr):
            if rr<r[0]:return mass[0]*(rr/r[0])**3
            if rr>r[-1]:return mass[-1]
            return float(mi(np.log(rr)))
        def bends(angle):
            b=angle/ARCSEC*dl
            star=4*G*row['mass_Msun']/C**2*quad(lambda t:model.mass_fraction(b/np.cos(t))/(b/np.cos(t)),0,np.pi/2,epsabs=1e-8,epsrel=1e-8)[0]
            comp=4*G/C**2*quad(lambda t:extra(b/np.cos(t))/(b/np.cos(t)),0,np.pi/2,epsabs=1e-6,epsrel=1e-8,limit=200)[0]
            return star*ratio*ARCSEC,comp*ratio*ARCSEC
        s0,c0=bends(row['lens_prediction_arcsec']);root_error=abs((s0+c0)/row['lens_prediction_arcsec']-1)
        assert root_error<1e-6,(name,branch,root_error)
        target=row['lens_catalog_arcsec'];star,comp=bends(target)
        sf=star/target;cf=comp/target
        out['rows'].append(dict(branch=branch,Name=name,population=row['retention_mapping']['population'],catalog_angle_arcsec=target,
            star_bending_arcsec=star,companion_bending_arcsec=comp,star_fraction_of_required_bending=sf,companion_fraction_of_required_bending=cf,
            required_companion_multiplier_at_fixed_stars=(1-sf)/cf,
            stellar_mass_ceiling_if_companions_removed_Msun=row['mass_Msun']/sf,
            minimum_stellar_mass_reduction_fraction_if_companions_removed=max(0.,1-1/sf),
            stellar_mass_multiplier_at_fixed_companions=(1-cf)/sf,
            root_relative_reproduction_error=root_error))
    for pop in ['Chabrier','Salpeter']:
        rr=[r for r in out['rows'] if r['branch']==branch and r['population']==pop]
        out['summary'].append(dict(branch=branch,population=pop,systems=len(rr),
            stars_alone_exceed_required_bending=sum(r['star_fraction_of_required_bending']>1 for r in rr),
            star_fraction_range=[min(r['star_fraction_of_required_bending'] for r in rr),max(r['star_fraction_of_required_bending'] for r in rr)],
            required_companion_multiplier_range=[min(r['required_companion_multiplier_at_fixed_stars'] for r in rr),max(r['required_companion_multiplier_at_fixed_stars'] for r in rr)]))
paths=[L/('capacity-'+b+'-optics-results.json') for b in branches]+[L/'third-radiation-retention-results.json',P.parent/'slacs-component-refit/model.py',P.parent/'slacs-resolved-fit/model.py']
paths += [P.parent/f/n for f,n in [('slacs-resolved-input-audit','results.json'),('slacs-light-profile-audit','results.json'),('slacs-motion-lensing-pilot','results.json'),('lensing-data-readiness','conditional-geometry.json')]]
out['input_sha256']={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
(P/'lens-bending-budget-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2))
