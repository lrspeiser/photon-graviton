"""One conditional mass-energy response for both motion and deflection."""
from pathlib import Path
import ast,json,hashlib
import numpy as np
from scipy.integrate import quad

OUT=Path(__file__).resolve().parent
projection_path=OUT/'projection.json'; data=json.loads(projection_path.read_text())
# Load projector definitions and inputs, stopping before its numerical driver.
path=OUT/'project.py'; tree=ast.parse(path.read_text()); nodes=[]
for node in tree.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='coarse' for t in node.targets): break
    nodes.append(node)
ns={'__file__':str(path)}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)
radial=ns['calculate'](600); q=radial[5]; total=data['independent_ray_power_Lsun']

def enclosed(r):
    if r>=1: return total
    return 4*np.pi*quad(lambda x:x*x*q(x),0,r,epsabs=1e-8,epsrel=2e-8,limit=200)[0]

G=6.67430e-11;c=299792458.;Lsun=3.828e26;Mpc=3.085677581491367e22
Gyr=31557600.*1e9;Msun=1.98847e30
mass_per_power=Lsun*Gyr/c**2
rows=[]
for p in data['profile']:
    b=p['projected_radius_kpc']/1000
    if b<.001: continue
    p3=p['sphere_enclosed_power_Lsun'];p2=p['projected_aperture_power_Lsun']
    # Integrate the gravitational field to infinity, including outside receiver.
    # z=b*tan(theta) => projected-equivalent power = integral M(b sec(theta))*cos(theta) dtheta.
    split=np.arccos(b) if b<1 else 0.
    line=quad(lambda theta:enclosed(b/np.cos(theta))*np.cos(theta),0,split,epsabs=.001,epsrel=2e-7,limit=150)[0]+total*(1-np.sin(split))
    err=abs(line/p2-1)
    assert err<2e-5
    v2=G*mass_per_power*p3/(b*Mpc)
    bend=4*G*mass_per_power*p2/(c*c*b*Mpc)
    ratio=4*p2/p3
    assert abs(bend/(v2/c**2)/ratio-1)<1e-14
    rows.append(dict(radius_kpc=b*1000,deposited_mass_equivalent_inside_sphere_Msun_per_Gyr=mass_per_power*p3/Msun,
        circular_speed_squared_km2_s2_per_Gyr=v2/1e6,deflection_arcsec_per_Gyr=bend*180*3600/np.pi,
        deflection_over_v_squared_over_c_squared=ratio,field_integral_vs_aperture_relative=err))
result=dict(scope='Conditional static, spherical, weak-field GR-like equal-potential response to stationary stored energy. Unit Gyr is a normalization, not an age. M87 deposits only; no total measured stellar speed or lensing comparison.',
 projection_sha256=hashlib.sha256(projection_path.read_bytes()).hexdigest(),source_total_mass_equivalent_Msun_per_Gyr=total*mass_per_power/Msun,
 max_independent_field_integral_relative=max(r['field_integral_vs_aperture_relative'] for r in rows),rows=rows)
(OUT/'response.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
