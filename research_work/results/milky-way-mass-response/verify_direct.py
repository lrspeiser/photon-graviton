"""Re-solve a rescaled field to check the analytic homogeneity shortcut."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'conservative-field-completion'
expected=json.loads((BASE/'verification.json').read_text())['input_sha256']
relative='research_work\\results\\conservative-field-completion\\run.py'
assert hashlib.sha256((BASE/'run.py').read_bytes()).hexdigest()==expected[relative]
spec=importlib.util.spec_from_file_location('conservative_field',BASE/'run.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
LAMBDA=.8
source=mod.Baryons(True)
class Scaled:
    def evaluate(self,r,mu):return LAMBDA*source.evaluate(r,mu)
print('Independent field solve for lambda=0.8',flush=True)
extra=mod.Completion(Scaled(),True,200)
original=json.loads((BASE/'predictions-refined.json').read_text())
checks={}
for obs in ['rotation','vertical']:
    b=[r for r in original if r['observable']==obs and r['model']=='ordinary']
    t=[r for r in original if r['observable']==obs and r['model']=='conservative_completion']
    R=np.array([r['R_kpc'] for r in b]);z=np.array([r['z_kpc'] for r in b])
    ordinary=mod.force(Scaled(),R,z);additional=mod.force(extra,R,z)
    full=ordinary+additional
    component=0 if obs=='rotation' else 1
    assert np.all(ordinary[:,component]<0) and np.all(full[:,component]<0)
    assert np.all(additional[:,component]<0)
    power=2 if obs=='rotation' else 1
    bb=np.array([r['predicted']**power for r in b]);cc=np.array([r['predicted']**power for r in t])-bb
    shortcut=(LAMBDA*bb+LAMBDA**mod.P*cc)**(1/power)
    direct=np.sqrt(-R*full[:,0]) if obs=='rotation' else abs(full[:,1])/(2*np.pi*mod.G*1e6)
    maximum=float(np.max(abs(shortcut-direct)))
    assert maximum<1e-7,maximum
    checks[obs]=dict(maximum_absolute_direct_shortcut_difference=maximum,
        all_ordinary_extra_total_accelerations_have_expected_inward_sign=True)
result=dict(lambda_mass=LAMBDA,checks=checks,
    parent_solver_hash_unchanged=True,new_observational_fit=False)
(HERE/'direct-verification.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
