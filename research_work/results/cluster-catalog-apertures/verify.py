"""Compare volume-sampled spherical apertures to separate ray integration."""
from pathlib import Path
import json,hashlib
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
result=json.loads((OUT/'results.json').read_text())
ray_path=ROOT/'research_work/results/cluster-catalog-deposit-profile/results.json'
ray=json.loads(ray_path.read_text())['runs'][0]
m87=json.loads((ROOT/'research_work/results/cluster-m87-finite-emitter/projection.json').read_text())
L=1.504e11;a=.0002488993265191759;k=10.
rows=[]
for row in result['rows']:
    r=row['radius_Mpc'];i=int(np.argmin(abs(np.array(ray['radii_Mpc'])-r)))
    assert abs(ray['radii_Mpc'][i]-r)<1e-12
    base=ray['internal_enclosed_power_Lsun'][i]+ray['external_enclosed_power_Lsun'][i]
    old=L*(-np.expm1(-a*r)-a*np.exp(-a*r)*(-np.expm1(-(k-a)*r))/(k-a))
    new=next(v['sphere_enclosed_power_Lsun'] for v in m87['profile'] if abs(v['projected_radius_kpc']/1000-r)<1e-12)
    independent=base-old+new
    difference=row['mean_sphere_power_Lsun']/independent-1
    assert abs(difference)<.03
    rows.append(dict(radius_Mpc=r,independent_ray_sphere_power_Lsun=independent,volume_mean_relative_difference=difference,
       conditional_ratio_using_ray_denominator=4*row['mean_projected_power_Lsun']/independent))
output=dict(ray_input_sha256=hashlib.sha256(ray_path.read_bytes()).hexdigest(),rows=rows,
    caveat='Four scramble ranges are not confidence intervals; smaller apertures converge more slowly. Ray and source extent approximations retain their original limitations.')
(OUT/'verification.json').write_bytes((json.dumps(output,indent=2)+'\n').encode())
print(json.dumps(output,indent=2))
