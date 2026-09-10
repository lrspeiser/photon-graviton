"""Verify the orbit ordinary field matches the completion's source field."""
from pathlib import Path
import json
import numpy as np
import run as orbit

HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(92167)
R=np.exp(rng.uniform(np.log(.01),np.log(20.),60))
z=rng.uniform(-2,2,60);phi=rng.uniform(-np.pi,np.pi,60)
points=np.c_[R*np.cos(phi),R*np.sin(phi),z]
rad=np.hypot(R,z);mu=z/rad
disk,_=orbit.full.disk_model()
nuc=orbit.full.old.CachedAxisymmetric(orbit.ROOT/'research_work/data-cache/bar-field/nuclei-L16.npz')
out=np.array(disk.evaluate(rad,mu))+np.array(nuc.evaluate(rad,mu))
potential,gr,gt=out
aR=-gr*R/rad-gt*z/rad
az=-gr*z/rad+gt*R/rad
acc=np.c_[aR*np.cos(phi),aR*np.sin(phi),az]
bar=orbit.full.load_field(orbit.ROOT/'research_work/data-cache/bar-field/bar-L64.npz')
for extra in [bar.evaluate(points),orbit.point_potential(points)]:
    potential+=extra[0];acc+=extra[1]
p,a=orbit.Field('ordinary').evaluate(points)
force_error=float(np.max(np.linalg.norm(acc-a,axis=1)/np.linalg.norm(acc,axis=1)))
potential_error=float(np.max(abs(p-potential))/np.max(abs(potential)))
assert force_error<1e-10 and potential_error<1e-10
result=dict(scope='Independent Legendre-to-real-harmonic source-field consistency at 60 synthetic points.',
            acceleration_fractional_error=force_error,potential_scaled_error=potential_error,
            no_observed_data_used=True)
(HERE/'source-field-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
