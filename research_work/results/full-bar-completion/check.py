"""Independent boundary/reference checks for the cached full-bar prototype."""
from pathlib import Path
import json
import hashlib
import numpy as np
import run as model

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'results.json').read_text())
points=np.array([r['position_kpc'] for r in data['rows']])
disk,dp=model.disk_model()
base_path=model.ROOT/next(k for k in data['input_hashes'] if 'axisymmetric-reference-80' in k)
base=model.load_field(base_path)
delta=model.load_field(model.CACHE/'finer.npz')
full=model.SumField(base,delta)
pot,acc=full.evaluate(points)
assert np.max(abs(acc-np.array([r['fine_extra_acceleration'] for r in data['rows']])))<1e-10

print('Refining axisymmetric reference only',flush=True)
high,hp=model.axis_field(disk,ellmax=256,angular_nodes=1024)
reference_change=np.linalg.norm(high.evaluate(points)[1]-base.evaluate(points)[1],axis=1)/np.linalg.norm(acc,axis=1)

print('Moving source boundary to 100 kpc',flush=True)
outer,op=model.axis_field(disk,rmax=100.)
outer_delta,odp=model.build('outer100',768,256,192,64,disk,rmax=100.)
outer_acc=model.SumField(outer,outer_delta).evaluate(points)[1]
boundary_change=np.linalg.norm(outer_acc-acc,axis=1)/np.linalg.norm(acc,axis=1)

axis_acc=base.evaluate(points)[1]
bar_delta=acc-axis_acc
bar_fraction=np.linalg.norm(bar_delta,axis=1)/np.linalg.norm(axis_acc,axis=1)
# Reflection-symmetric source should generate a reflection-symmetric force.
mirror=points.copy();mirror[:,2]*=-1
ma=full.evaluate(mirror)[1];ma[:,2]*=-1
reflection_error=float(np.max(np.linalg.norm(ma-acc,axis=1)/np.linalg.norm(acc,axis=1)))
assert reflection_error<1e-10

out=dict(scope='Reference-order and finite-boundary sensitivities, not continuum proof or observational fit.',
         base_reference_order=128,refined_reference_order=256,
         reference_force_change_max=float(reference_change.max()),
         boundary_radii_kpc=[80,100],boundary_force_change_max=float(boundary_change.max()),
         reflection_error=reflection_error,
         all_checked_force_changes_below_one_percent=bool(max(reference_change.max(),boundary_change.max(),data['force_refinement_maximum'])<.01),
         disk_field_refined=False,observer_orbits_fitted=False,
         rows=[dict(position_kpc=x.tolist(),reference_change_fraction=float(a),boundary_change_fraction=float(b),
                    full_bar_extra_force_change_relative_axis=float(c)) for x,a,b,c in zip(points,reference_change,boundary_change,bar_fraction)],
         input_hashes={str(p.relative_to(model.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       [Path(__file__),HERE/'run.py',HERE/'results.json',base_path,model.CACHE/'finer.npz',hp,op,odp,dp]})
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Reference change max',reference_change.max(),flush=True)
print('Boundary change max',boundary_change.max(),flush=True)
print('Bar correction max',bar_fraction.max(),flush=True)
