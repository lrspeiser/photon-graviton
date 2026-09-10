"""Independent disk-mass integrals and metric/source-integrity checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.special import kv

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
s=json.loads((HERE/'results.json').read_text());rows=json.loads((HERE/'predictions.json').read_text())
for name,digest in s['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
params=[(1.332e9,2.,2.7),(8.97e8,2.8,2.7),(5.81e7,7.,4.),(2.68e9,1.5,12.)]
# Integrating the normalized vertical profile analytically leaves the holed exponential surface density.
mass=[4*np.pi*sigma*rd*hole*kv(2,2*np.sqrt(hole/rd)) for sigma,rd,hole in params]
exact={'stellar':sum(mass[:2]),'gas':sum(mass[2:])}
mass_error={name:abs(s['baseline_disk_masses_Msun'][name]/value-1) for name,value in exact.items()}
assert max(mass_error.values())<.001
errors=[]
for case,key in [('refined_original','original'),('refined_component_candidate','candidate'),('refined_balanced_candidate','balanced')]:
    trace=next(t for t in s['trace'] if t['label']==case)
    for obs,metric in [('rotation_training','rotation_rms_kms'),('vertical_exposed','vertical_rms_surface_equivalent')]:
        selected=[r for r in rows if r['observable']==obs]
        value=np.sqrt(np.mean([(r[key]-r['observed'])**2 for r in selected]))
        errors.append(abs(value-trace[metric]));assert abs(value-trace[metric])<1e-10
assert all(.7-1e-12<=x<=1.3+1e-12 for x in s['candidate_scales'])
assert all(.7-1e-12<=x<=1.3+1e-12 for x in s['balanced_candidate_scales'])
assert s['checks']['max_candidate_coarse_refined_rotation_kms']<.1
assert s['checks']['max_candidate_coarse_refined_vertical_units']<1.
out=dict(independent_analytic_disk_masses_Msun=exact,relative_mass_errors=mass_error,
    max_stored_metric_error=max(errors),source_hashes_unchanged=True,
    candidate_in_declared_sensitivity_box=True,
    qualification='Numerical decomposition and mass integration checks; does not validate the illustrative mass bounds or establish a global fit.')
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
