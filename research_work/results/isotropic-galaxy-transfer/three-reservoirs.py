"""Calibrate common storage normalization on Sun; predict Earth and Moon."""
from pathlib import Path
import json,numpy as np
HERE=Path(__file__).resolve().parent
s=json.loads((HERE/'solar-photon-budget-results.json').read_text())
p=json.loads((HERE/'planet-well-energy-results.json').read_text())
data=[dict(body='Sun',mass=1.9884e30,radius=6.957e8,power=s['nominal_solar_luminosity_W'])]
for name in ['Earth','Moon']:
    r=next(x for x in p['rows'] if x['body']==name)
    data.append(dict(body=name,mass=r['mass_kg'],radius=r['radius_m'],power=r['intercepted_solar_power_W']))
c=299792458.;year=31557600.;alpha=s['alpha_per_Mpc']/3.085677581491367e22
out=dict(status='Restricted geometries/source proxies; shared Sun calibration with exposed Earth/Moon transfer; self-bound independent branch unresolved',cases={})
for mode in ['fitted_alpha_one_radius','full_conversion_ceiling']:
    rows=[]
    for d in data:
        f=-np.expm1(-alpha*d['radius']) if mode.startswith('fitted') else 1.
        pc=d['power']*f;target=d['mass']*c*c
        rows.append(dict(body=d['body'],source_power_W=d['power'],D_m=d['radius'],converted_fraction=f,
                         companion_power_W=pc,target_mass_kg=d['mass'],target_energy_J=target,
                         required_K_seconds=target/pc,required_K_years=target/pc/year,
                         free_escape_time_seconds=d['radius']/c,
                         free_escape_reservoir_mass_upper_kg=pc*d['radius']/c**3))
    K=rows[0]['required_K_seconds']
    for r in rows:
        predicted=K*r['companion_power_W']/c**2
        r.update(shared_K_predicted_mass_kg=predicted,predicted_over_target_mass=predicted/r['target_mass_kg'],
                 independent_required_K_over_solar_K=r['required_K_seconds']/K,
                 free_escape_upper_over_target_mass=r['free_escape_reservoir_mass_upper_kg']/r['target_mass_kg'])
    assert abs(rows[0]['predicted_over_target_mass']-1)<1e-12
    out['cases'][mode]=dict(K_calibrated_on_Sun_seconds=K,rows=rows)
(HERE/'three-reservoirs-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
for mode,case in out['cases'].items():
    print(mode)
    for r in case['rows']:print(r['body'], 'power',r['companion_power_W'],'K_years',r['required_K_years'],'pred/actual',r['predicted_over_target_mass'],'free_mass',r['free_escape_reservoir_mass_upper_kg'])
