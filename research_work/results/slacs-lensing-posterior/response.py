"""Required extra-sector lens response at published lens angle, fixed stellar inference."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'slacs-component-refit'))
from model import ComponentModel,ARCSEC,G,C
def read(folder,file='results.json'):return json.loads((HERE.parent/folder/file).read_text())
data=read('slacs-resolved-input-audit');photo=read('slacs-light-profile-audit')
base=read('slacs-motion-lensing-pilot');geometry=read('lensing-data-readiness','conditional-geometry.json')
post=json.loads((HERE/'refined-results.json').read_text())
profiles={r['Name']:r for r in photo['rows']};geo={r['Name']:r for r in geometry}
posterior={r['Name']:r for r in post['rows']}
pilot={r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
par=base['parameters'];amp,p=par['A'],par['p'];astar=par['a_star_m_s2']*3.085677581491367e19/1e6
rows=[]
for item in data['systems']:
    name=item['Name']
    if name not in posterior:continue
    dl=geo[name]['conditional_Dl_Mpc']*1000;ratio=geo[name]['conditional_Dls_over_Ds']
    a=pilot[name]['scale_a_kpc'];edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*dl/ARCSEC
    psf=item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
    comp=[dict(R=q['R_arcsec']*dl/ARCSEC,n=q['n'],amp=q['amp_at_R'],bn=q['bn']) for q in profiles[name]['components']]
    model=ComponentModel(a,edges,psf,amp,p,astar,20,comp)
    selected=next(q for q in posterior[name]['predictions'] if q['model']=='empirical_extra' and q['prior']=='uniform_log_mass_and_beta')
    masses=selected['resolutions'][-1]['mass_quantiles_Msun']
    theta=posterior[name]['catalog_SIE_arcsec'];b=theta*dl/ARCSEC;rt=model.rt
    def bending(mass):
        gc_edge=amp*astar*(G*mass*model.mass_fraction(rt)/rt**2/astar)**p
        def fun(t,extra):
            r=b/np.cos(t);gb=G*mass*model.mass_fraction(r)/r**2
            g=(amp*astar*(gb/astar)**p if r<=rt else gc_edge*(rt/r)**2) if extra else gb
            return g*b/np.cos(t)
        edge=np.arccos(b/rt) if b<rt else None
        return [ratio*ARCSEC*4/C**2*quad(lambda t:fun(t,extra),0,np.pi/2,
            points=[edge] if edge is not None else None,epsabs=1e-6,epsrel=1e-7,limit=200)[0] for extra in [False,True]]
    values=[]
    for mass in [masses[0],masses[2],masses[-1]]:
        bb,bc=bending(mass);f=(theta-bb)/bc
        values.append(dict(stellar_mass=mass,baryon_deflection_arcsec=bb,extra_deflection_at_f1_arcsec=bc,
                           required_f=f,required_gamma_extra=2*f-1))
    bounds=sorted([values[0]['required_f'],values[-1]['required_f']])
    mid=values[1]
    rows.append(dict(Name=name,catalog_SIE_arcsec=theta,mass_conditioned_values=values,
        required_f_range_from_central95_mass=bounds,
        fixed_gamma0_lens_equation_residual_arcsec=mid['baryon_deflection_arcsec']+.5*mid['extra_deflection_at_f1_arcsec']-theta,
        fixed_gamma1_lens_equation_residual_arcsec=mid['baryon_deflection_arcsec']+mid['extra_deflection_at_f1_arcsec']-theta))
lower=max(r['required_f_range_from_central95_mass'][0] for r in rows)
upper=min(r['required_f_range_from_central95_mass'][1] for r in rows)
out=dict(scope='Exposed six-system spherical lens-summary diagnostic, no joint image likelihood',
    ordinary_matter_lensing_factor_fixed=1,extra_factor_definition='f=(1+gamma_extra)/2',
    stellar_inference_unchanged=True,lens_angles_used_to_infer_required_response_not_blind_predictions=True,
    intersection_lower=lower,intersection_upper=upper,common_range_exists=lower<=upper,rows=rows)
(HERE/'response-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
