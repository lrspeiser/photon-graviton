"""Rebuild all six physical models under a consistent common distance rescaling."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/'slacs-component-refit'))
from model import ComponentModel, ARCSEC
paths = [HERE.parent/'slacs-resolved-input-audit/results.json',
         HERE.parent/'slacs-light-profile-audit/results.json',
         HERE.parent/'slacs-motion-lensing-pilot/results.json',
         HERE.parent/'lensing-data-readiness/conditional-geometry.json',
         HERE.parent/'slacs-lensing-posterior/refined-results.json']
data, photo, base, geometry, posterior = [json.loads(f.read_text(encoding='utf-8')) for f in paths]
profiles = {r['Name']:r for r in photo['rows']}
geo = {r['Name']:r for r in geometry}
post = {r['Name']:r for r in posterior['rows']}
pilot = {r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
par = base['parameters']
A, p, astar = par['A'], par['p'], par['a_star_m_s2']*3.085677581491367e19/1e6
rows = []
for item in data['systems']:
    name = item['Name']
    if name not in post:
        continue
    Dl = geo[name]['conditional_Dl_Mpc']*1000
    ratio = geo[name]['conditional_Dls_over_Ds']
    a = pilot[name]['scale_a_kpc']
    edges = np.r_[item['inner_arcsec'], item['outer_arcsec'][-1]]*Dl/ARCSEC
    psf = item['psf_fwhm_arcsec']*Dl/ARCSEC/np.sqrt(8*np.log(2))
    comp = [dict(R=c['R_arcsec']*Dl/ARCSEC, n=c['n'], amp=c['amp_at_R'], bn=c['bn']) for c in profiles[name]['components']]
    reference = ComponentModel(a, edges, psf, A, p, astar, 20, comp)
    masses = {r['model']:r['resolutions'][-1]['mass_quantiles_Msun'][2] for r in post[name]['predictions'] if r['prior']=='uniform_log_mass_and_beta'}
    for scale in [.1, 10.]:
        scaled_comp = [dict(c, R=c['R']*scale) for c in comp]
        changed = ComponentModel(a*scale, edges*scale, psf*scale, A, p, astar/scale, 20, scaled_comp)
        for label, extra in [('baryons', False), ('empirical_extra', True)]:
            mass = masses[label]
            ref_angle = reference.angle(mass, Dl, ratio, extra)
            angle = changed.angle(mass*scale, Dl*scale, ratio, extra)
            force_i = 1 if extra else 0
            # Per-bin velocity predictions independently rebuilt at three orbits.
            errors = []
            for beta in [-.5, 0., .3]:
                original = reference.predict(mass, beta, extra)
                transformed = changed.predict(mass*scale, beta, extra)
                errors.append(float(np.max(np.abs(transformed/original-1))))
            rows.append(dict(Name=name, scale=scale, model=label,
                             max_relative_velocity_error=max(errors),
                             reference_angle_arcsec=ref_angle, transformed_angle_arcsec=angle,
                             angle_error_arcsec=abs(angle-ref_angle)))
out = dict(scope='Conditional common-distance scaling invariance; no new fit', scales=[.1,10.],
           beta_checks=[-.5,0,.3], rows=rows,
           gates=dict(max_relative_velocity_error=1e-6, max_angle_error_arcsec=1e-5),
           hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths+[HERE/'run.py', HERE.parent/'slacs-component-refit/model.py', HERE.parent/'slacs-resolved-fit/model.py']})
out['passed'] = all(r['max_relative_velocity_error'] < 1e-6 and r['angle_error_arcsec'] < 1e-5 for r in rows)
(HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
assert out['passed']
print(json.dumps(dict(passed=out['passed'], cases=len(rows), max_velocity_error=max(r['max_relative_velocity_error'] for r in rows), max_angle_error=max(r['angle_error_arcsec'] for r in rows))))
