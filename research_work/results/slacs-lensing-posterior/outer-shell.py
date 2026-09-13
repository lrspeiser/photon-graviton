"""Required outer-shell addition from existing stellar-conditioned lens residuals."""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
response=json.loads((HERE/'response-results.json').read_text())
geo={x['Name']:x for x in json.loads((HERE.parent/'lensing-data-readiness/conditional-geometry.json').read_text())}
pilot=json.loads((HERE.parent/'slacs-motion-lensing-pilot/results.json').read_text())
scales={x['Name']:x['scale_a_kpc'] for x in pilot['rows'] if x['model']=='empirical_extra' and x['cutoff_in_a']==20}
G=4.300917270036279e-6;c=299792.458;arcsec=206264.80624709636
rows=[]
for r in response['rows']:
    name=r['Name'];dl=1000*geo[name]['conditional_Dl_Mpc'];ratio=geo[name]['conditional_Dls_over_Ds']
    theta=r['catalog_SIE_arcsec'];b=theta*dl/arcsec;radius=20*scales[name]
    assert b<radius
    x=b/radius;projected_fraction=x*x/(1+np.sqrt(1-x*x))
    factor=ratio*arcsec*4*G*projected_fraction/(c*c*b)
    required=[]
    for item in r['mass_conditioned_values']:
        missing=theta-item['baryon_deflection_arcsec']-item['extra_deflection_at_f1_arcsec']
        m=missing/factor
        required.append(dict(stellar_mass=item['stellar_mass'],missing_deflection_arcsec=missing,
            formal_signed_shell_mass_Msun=m,positive_addition_possible=bool(m>=0),
            signed_capture_column_Msun_per_kpc2=m/(np.pi*radius**2)))
    median=required[1]
    recon=factor*median['formal_signed_shell_mass_Msun']
    assert abs(recon-median['missing_deflection_arcsec'])<1e-12
    rows.append(dict(Name=name,impact_radius_kpc=b,shell_radius_kpc=radius,
        projected_shell_fraction=projected_fraction,
        required_at_mass_quantiles=required,
        median_required_external_convergence=median['missing_deflection_arcsec']/theta))
out=dict(scope='Conditional addition to fixed equal-potential lens model; no new stellar refit or physical shell support',
    shell_radius_rule='20 times retained scale a',shell_lensing_gamma=1,
    signed_negative_masses_are_impossible_additions_not_negative_energy_proposals=True,rows=rows)
(HERE/'outer-shell-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
for r in rows:
    q=r['required_at_mass_quantiles'][1]
    print(r['Name'],r['shell_radius_kpc'],q['formal_signed_shell_mass_Msun'],r['median_required_external_convergence'])
