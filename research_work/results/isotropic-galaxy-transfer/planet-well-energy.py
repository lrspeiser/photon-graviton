"""Common spherical potential and solar-fed increment for planets and Moon."""
from pathlib import Path
import json,numpy as np
HERE=Path(__file__).resolve().parent
solar=json.loads((HERE/'solar-photon-budget-results.json').read_text())
G=6.67430e-11;c=299792458.;year=31557600.;Mpc=3.085677581491367e22;Rsun=6.957e8
# NASA summary masses (1e24 kg), diameters (km), solar distances (1e6 km).
# Moon solar distance uses Earth's; its 0.384 table entry is geocentric.
inputs=[('Mercury',.330,4879,57.9),('Venus',4.87,12104,108.2),
        ('Earth',5.97,12756,149.6),('Moon',.073,3475,149.6),
        ('Mars',.642,6792,228.),('Jupiter',1898,142984,778.5),
        ('Saturn',568,120536,1432.),('Uranus',86.8,51118,2867.),
        ('Neptune',102,49528,4515.),('Pluto',.0130,2376,5906.4)]
rows=[]
for name,m,d,a in inputs:
    mass=m*1e24;r=d*500;distance=a*1e9
    depth=G*mass/r;power=solar['nominal_solar_luminosity_W']*r*r/(4*distance**2)
    f=-np.expm1(-solar['alpha_per_Mpc']*(distance-Rsun)/Mpc)
    delta=G*power*year/(c*c*r)
    rows.append(dict(body=name,mass_kg=mass,radius_m=r,solar_distance_m=distance,
        self_well_depth_J_kg=depth,spherical_surface_g_m_s2=depth/r,
        spherical_escape_speed_km_s=np.sqrt(2*depth)/1000,
        intercepted_solar_power_W=power,converted_fraction=f,
        all_incident_energy_retained_depth_gain_per_year_J_kg=delta,
        fitted_conversion_perfect_capture_depth_gain_per_year_J_kg=delta*f,
        all_incident_retained_mass_equivalent_kg_per_year=power*year/c**2,
        fitted_conversion_retained_mass_equivalent_kg_per_year=power*year/c**2*f,
        effective_full_illumination_years_for_total_mass=mass*c*c/power/year))
earth=next(r for r in rows if r['body']=='Earth');moon=next(r for r in rows if r['body']=='Moon')
out=dict(status='Known potential consistency plus hypothetical retained-energy increments; not an independent gravity validation',
    sources=['https://nssdc.gsfc.nasa.gov/planetary/factsheet/'],rows=rows,
    Earth_to_Moon_existing_depth_ratio=earth['self_well_depth_J_kg']/moon['self_well_depth_J_kg'],
    Earth_to_Moon_solar_fed_depth_ratio=earth['all_incident_energy_retained_depth_gain_per_year_J_kg']/moon['all_incident_energy_retained_depth_gain_per_year_J_kg'],
    Earth_to_Moon_required_efficiency_times_duration_ratio=earth['effective_full_illumination_years_for_total_mass']/moon['effective_full_illumination_years_for_total_mass'])
assert all(abs(r['spherical_escape_speed_km_s']**2*1e6/(2*r['self_well_depth_J_kg'])-1)<1e-12 for r in rows)
(HERE/'planet-well-energy-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
for r in rows: print(r['body'],*[f'{r[k]:.6g}' for k in ['self_well_depth_J_kg','spherical_surface_g_m_s2','intercepted_solar_power_W','all_incident_energy_retained_depth_gain_per_year_J_kg','fitted_conversion_perfect_capture_depth_gain_per_year_J_kg']])
print({k:v for k,v in out.items() if k.startswith('Earth')})
