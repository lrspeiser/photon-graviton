"""Finite stored-mass and energy-equivalent bounds for the fitted companion profiles."""
from pathlib import Path
import json,hashlib,math
P=Path(__file__).resolve().parent
source=P/'free-companion-results.json';fit=json.loads(source.read_text())
rows=[]
for v in fit['rows']:
    a=v['capture_scale_kpc'];D=v['stored_density_normalization_Msun_kpc3']
    R=a*v['outer_grid_radius_in_capture_scales'];M=v['stored_mass_within_grid_Msun']
    # J<=1 and [1+(r/a)^2]^-2 <= a^4/r^4 imply this exterior-mass bound.
    tail=4*math.pi*D*a**4/R
    unit=1.98847e30*299792458.**2
    rows.append(dict(Name=v['Name'],integration_radius_kpc=R,mass_inside_grid_Msun=M,exterior_mass_upper_bound_Msun=tail,total_mass_upper_bound_Msun=M+tail,exterior_bound_fraction_of_grid_mass=tail/M,energy_equivalent_inside_grid_J=M*unit,energy_equivalent_total_upper_bound_J=(M+tail)*unit))
out=dict(scope='Reference effective mass equals deposited energy divided by c squared; not a unique local gravitational-field energy; no universe-age assumption',solar_mass_kg=1.98847e30,speed_of_light_m_s=299792458.,rows=rows,input_sha256={source.name:hashlib.sha256(source.read_bytes()).hexdigest()})
(P/'free-companion-inventory-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(rows,indent=2))
