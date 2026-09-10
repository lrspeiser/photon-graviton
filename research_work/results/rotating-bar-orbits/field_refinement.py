"""Resolve the prior failed high-order bar calculation; report field uncertainty."""
import json
import time
import numpy as np
from run import HERE,FOUNDATION,FIELD_CACHE,save,digest
from field import Multipole,bar,disks,disk_evaluate

started=time.monotonic()
path=FIELD_CACHE/'bar-L64.npz'
meta=FIELD_CACHE/'bar-L64-provenance.json'
source_hash=digest(FOUNDATION/'field.py')
density_hash=digest(FOUNDATION/'published_bar.py')
settings=dict(lmax=64,nr=1600,ntheta=224,nphi=288)
expected=dict(field_source_sha256=source_hash,density_source_sha256=density_hash,settings=settings)
if path.exists():
    assert meta.exists() and json.loads(meta.read_text())==expected,'Unverified high-order cache'
    high=Multipole.load(path)
else:
    print('Building bar L64 with stable radial integration',flush=True)
    high=Multipole.build(bar,**settings)
    assert np.isfinite(high.coeff).all()
    high.save(path)
    meta.write_text(json.dumps(expected,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Bar L64 built after',round(time.monotonic()-started,1),'seconds',flush=True)
records=json.loads((FOUNDATION/'field-predictions.json').read_text())
baseline=[r for r in records if r['model']=='ordinary_matter']
xyz=np.array([[r['x_kpc'],r['y_kpc'],r['z_kpc']] for r in baseline])
base_a=np.array([r['acceleration_kms2_per_kpc'] for r in baseline])
low=Multipole.load(FIELD_CACHE/'bar-L40.npz')
_,aa=low.evaluate(xyz);_,bb=high.evaluate(xyz)
diff=np.linalg.norm(bb-aa,axis=1)
relative=diff/np.linalg.norm(bb,axis=1)
total=diff/np.linalg.norm(base_a,axis=1)
result=dict(bar_L40_to_L64=dict(max_fraction_of_bar_force=float(relative.max()),
    median_fraction_of_bar_force=float(np.median(relative)),max_fraction_of_previous_total_force=float(total.max()),
    worst_position_kpc=xyz[relative.argmax()].tolist(),mass_Msun=high.integrated_mass),
    provenance=expected,bar_cache_sha256=digest(path))
save('field-refinement.json',result)
forces={}
for order in [24,40,64]:
    print('Disk order',order,flush=True)
    potential=disks(order)
    _,forces[order]=disk_evaluate(potential,xyz)
for a,b in [(24,40),(40,64)]:
    frac=np.linalg.norm(forces[b]-forces[a],axis=1)/np.linalg.norm(forces[b],axis=1)
    result[f'disk_N{a}_to_N{b}']=dict(max_fraction_of_disk_force=float(frac.max()),
         median_fraction_of_disk_force=float(np.median(frac)),worst_position_kpc=xyz[frac.argmax()].tolist(),
         max_fraction_of_previous_total_force=float(np.max(np.linalg.norm(forces[b]-forces[a],axis=1)/np.linalg.norm(base_a,axis=1))))
result['elapsed_seconds']=time.monotonic()-started
result['status']='Field-resolution sensitivity on existing 72-probe grid; no orbital likelihood or observational fit'
save('field-refinement.json',result)
print(json.dumps(result,indent=2),flush=True)
